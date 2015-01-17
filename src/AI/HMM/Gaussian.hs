{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE OverloadedLists #-}
module AI.HMM.Gaussian where

import Algorithms.GLasso (glasso)
import Numeric.LinearAlgebra.HMatrix (Matrix, Vector, (<>), invlndet, (!), tr, asRow, vector, matrix, reshape, size, flatten)
import Control.Monad.Primitive (PrimMonad, PrimState)
import Control.Monad (forM_, foldM)
import Control.Monad.ST (ST,runST)
import Data.Function (on)
import Data.List (groupBy, sortBy)
import qualified Data.Matrix.Unboxed as M
import qualified Data.Matrix.Unboxed.Mutable as MM
import Data.Ord (comparing)
import Data.STRef
import qualified Data.Vector as V
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Generic.Mutable as GM
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Unboxed.Mutable as UM
import Statistics.Sample (meanWeighted, mean)
import System.Random.MWC

import AI.Clustering.KMeans
import AI.Function
import AI.HMM.Class

import Debug.Trace


data GaussianHMM = GaussianHMM
    { _startProb :: !(U.Vector Double)
    , _transitonProb :: !(M.Matrix Double)
    , _emission :: !(V.Vector MVN)
    } deriving (Show)

instance HMM GaussianHMM (Vector Double) where
    nSt (GaussianHMM x _ _) = U.length x
    {-# INLINE nSt #-}

    π (GaussianHMM x _ _) s = x `U.unsafeIndex` s
    {-# INLINE π #-}

    a (GaussianHMM _ x _) i j = x `M.unsafeIndex` (i,j)
    {-# INLINE a #-}

    b h i o = exp (b' h i o)

    b' (GaussianHMM _ _ x) i = logPDF (x `G.unsafeIndex` i)
    {-# INLINE b #-}

    baumWelch ob h@(GaussianHMM _ trans _) = GaussianHMM ini' trans' em'
      where
        ini' = G.zipWith (\x y -> exp $ x + y - G.head scales) (fw `M.takeColumn` 0) (bw `M.takeColumn` 0)
        trans' = MM.create $ do
            mat <- MM.new (n,n)
            forM_ [0..n-1] $ \i ->
                forM_ [0..n-1] $ \j -> do
                    let a_ij = trans `M.unsafeIndex` (i,j)
                    temp <- UM.new $ G.length ob
                    forM_ [1 .. G.length ob - 1] $ \t -> do
                        let o = ob `G.unsafeIndex` t
                            b_jo = b' h j o
                            α_it' = fw `M.unsafeIndex` (i,t-1)
                            β_jt = bw `M.unsafeIndex` (j,t)
                        GM.unsafeWrite temp t $ b_jo + α_it' + β_jt
                    x <- logSumExpM temp
                    MM.unsafeWrite mat (i,j) $ log a_ij + x
            normalizeByRow n n mat
            return mat

        em' = G.generate n $ \i -> let ws = G.map f $ G.enumFromN 0 (G.length ob)
                                       f t = let α_it = fw `M.unsafeIndex` (i,t)
                                                 β_it = bw `M.unsafeIndex` (i,t)
                                                 sc = scales `G.unsafeIndex` t
                                              in exp $ α_it + β_it - sc
                                       (m, cov) = weightedMeanCovMatrix ws ob'
                                   in mvn m (convert cov) -- $ fst $ glasso cov 0.01)

        (fw, scales) = forward' h ob
        bw = backward' h ob scales
        n = nSt h
        ob' = M.fromRows $ G.toList ob

convert mat = reshape c $ M.flatten mat
  where
    c = M.cols mat

-- normalize log probabilities
normalizeByRow :: Int -> Int -> MM.MMatrix s Double -> ST s ()
normalizeByRow r c mat = forM_ [0..r-1] $ \i -> do
    sc <- logSumExpM $ MM.takeRow mat i
    forM_ [0..c-1] $ \j -> MM.unsafeRead mat (i,j) >>= MM.unsafeWrite mat (i,j) . (exp . subtract sc)

-- | construct inital HMM model by kmeans clustering
kMeansInitial :: (PrimMonad m, G.Vector v (Vector Double))
              => Gen (PrimState m)
              -> v (Vector Double)
              -> Int
              -> m GaussianHMM
kMeansInitial g ob k = do 
    (membership, _) <- kmeans g k . M.fromRows . G.toList $ ob
    let pi = G.map (/ fromIntegral n) $ G.create $ do
            vec <- GM.replicate k 0
            V.forM_ membership $ \s -> GM.unsafeRead vec s >>= GM.unsafeWrite vec s . (+1)
            return vec

        trans = MM.create $ do
            mat <- MM.replicate (k,k) 0
            V.sequence_ . V.zipWith ( \i j -> MM.unsafeRead mat (i,j) >>=
                MM.unsafeWrite mat (i,j) . (+1) ) membership . V.tail $
                membership
            normalizeByRow k k mat
            return mat
        
        emisson = G.fromList $ map (uncurry mvn . meanCov) clusters
          where
            clusters = map (M.fromRows . snd . unzip) . groupBy ((==) `on` fst) . sortBy (comparing fst) $ zip (G.toList membership) $ G.toList ob
    return $ GaussianHMM pi trans emisson
  where
    n = G.length ob

meanCov dat = (meanVec, reshape p $ M.flatten $ fst $ glasso covs 0.01)
  where
    covs = MM.create $ do
        mat <- MM.new (p,p)
        forM_ [0..p-1] $ \i -> 
            forM_ [i..p-1] $ \j -> do
                let cov = covWithMean (meanVec G.! i, dat `M.takeColumn` i) (meanVec G.! j, dat `M.takeColumn` j) 
                MM.unsafeWrite mat (i,j) cov
                MM.unsafeWrite mat (j,i) cov
        return mat

    meanVec = G.fromList . map mean . M.toColumns $ dat
    p = G.length meanVec

covWithMean :: (Double, Vector Double) -> (Double, Vector Double) -> Double
covWithMean (mx, xs) (my, ys) | n == 1 = 0
                              | otherwise = G.sum (G.zipWith f xs ys) / (n - 1)
  where
    f x y = (x - mx) * (y - my)
    n = fromIntegral $ G.length xs

train :: V.Vector (Vector Double) -> Int -> Int -> GaussianHMM
train ob k n = run init 0
  where
    run !hmm !i | i >= n = hmm
                | otherwise = run (baumWelch ob hmm) (i+1)
    init = runST $ do
        g <- create
        kMeansInitial g ob k

{-
kmeansHMM :: V.Vector (Vector Double) -> Int -> GaussianHMM
kmeansHMM ob k = 
  where
    clusters = kme
    -}

-- test:
-- Gaussian Hidden Markov Model with 2 States
--
-- Transition matrix:
--    0.1 0.9
--    0.5 0.5
--
-- Emission parameters:
-- [N(1.0, 1.0), N(-1.0, 1.0)]
--
-- Initial probabilities: [0.5000, 0.5000]
--
-- viterbi path
-- ([1, 0, 1, 0, 1, 1, 0, 1, 0, 1], -16.67738270170788)
--
-- P(O) = -15.226179570606487
test = do
   let hmm = GaussianHMM (U.fromList [0.5,0.5]) (M.fromLists [[0.1,0.9],[0.5,0.5]]) (V.fromList [m1,m2])
       m1 = mvn (vector [1]) (matrix 1 [1])
       m2 = mvn (vector [-1]) (matrix 1 [1])
       obs = V.fromList $ map (vector.return) [-1.6835, 0.0635, -2.1688, 0.3043, -0.3188, -0.7835, 1.0398, -1.3558, 1.0882, 0.4050]
--   print $ viterbi hmm obs
--   print $ forward' hmm obs
       (f, sc) = forward hmm obs
       b = backward hmm obs sc
   loop obs hmm 0
 where
   loop o h i | i > 400 = 1
              | otherwise = let h' = baumWelch o h
                                (_, sc) = forward h o
                            in traceShow (loglikFromScales sc) $ loop o h' (i+1)

