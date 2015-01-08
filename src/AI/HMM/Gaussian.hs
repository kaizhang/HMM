{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE OverloadedLists #-}
module AI.HMM.Gaussian where

import Numeric.LinearAlgebra.HMatrix (Matrix, Vector, (<>), invlndet, (!), tr, asRow, vector, matrix, reshape, size, flatten)
import AI.Clustering.KMeans
import AI.HMM.Class
import Algorithms.GLasso (glasso)
import Control.Monad.Primitive (PrimMonad, PrimState)
import Control.Monad (forM_)
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
import Statistics.Sample (meanWeighted, mean)
import System.Random.MWC

import Debug.Trace

-- | multivariate normal distribution
data MVN = MVN
    { _mean :: !(Vector Double)
    , _cov :: !(Matrix Double)
    , _invcov :: !(Matrix Double)
    , _logdet :: !Double  -- ^ log determinant of covariance matrix
    } deriving (Show)

mvn :: Vector Double -> Matrix Double -> MVN
mvn m cov = MVN m cov invcov logdet
  where
    (invcov, (logdet, _)) = invlndet cov

-- | log probability of MVN
logPDF :: MVN -> Vector Double -> Double
logPDF (MVN m _ invcov logdet) x = -0.5 * ( d * log (2*pi) + logdet
                                        + (x' <> invcov <> tr x') ! 0 ! 0
                                        )
  where
    x' = asRow $ x - m
    d = fromIntegral . G.length $ m
{-# INLINE logPDF #-}

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
        ini' = G.map (/(G.head scales)) $ G.zipWith (*) (fw `M.takeColumn` 0) (bw `M.takeColumn` 0)
        trans' = MM.create $ do
            mat <- MM.replicate (n,n) 0
            forM_ [1 .. G.length ob - 1] $ \t -> do
                let o = ob `G.unsafeIndex` t
                forM_ [0..n-1] $ \i -> do
                    let α_it' = fw `M.unsafeIndex` (i,t-1)
                    forM_ [0..n-1] $ \j -> do
                        let a_ij = trans `M.unsafeIndex` (i,j)
                            b_jo = b h j o
                            β_jt = bw `M.unsafeIndex` (j,t)
                        MM.unsafeRead mat (i,j) >>= MM.unsafeWrite mat (i,j) .
                            (+) (α_it' * a_ij * b_jo * β_jt)
            normalizeByRow n n mat
            return mat

        em' = G.generate n $ \i -> let ws = G.map f $ G.enumFromN 0 (G.length ob)
                                       f t = let α_it = fw `M.unsafeIndex` (i,t)
                                                 β_it = bw `M.unsafeIndex` (i,t)
                                                 sc = scales `G.unsafeIndex` t
                                              in α_it * β_it / sc
                                       (m, cov) = weightedMeanCovMatrix ws ob'
                                   in mvn m (convert $ fst $ glasso cov 0.01)

        (fw, scales) = forward h ob
        bw = backward h ob scales
        n = nSt h
        ob' = M.fromRows $ G.toList ob

convert mat = reshape c $ M.flatten mat
  where
    c = M.cols mat

normalizeByRow :: Int -> Int -> MM.MMatrix s Double -> ST s ()
normalizeByRow r c mat = forM_ [0..r-1] $ \i -> do
    temp <- newSTRef 0
    forM_ [0..c-1] $ \j -> MM.unsafeRead mat (i,j) >>= modifySTRef' temp . (+) 
    s <- readSTRef temp
    forM_ [0..c-1] $ \j -> MM.unsafeRead mat (i,j) >>= MM.unsafeWrite mat (i,j) . (/s)

covWeighted :: U.Vector Double -> (Vector Double, Double) -> (Vector Double, Double) -> Double
covWeighted ws (xs, mx) (ys, my) = uncurry (/) . G.foldl f (0,0) $ U.enumFromN 0 $ G.length ws
  where
    f (!a,!b) i = let w = ws G.! i
                      x = xs G.! i
                      y = ys G.! i
                  in (a + w * (x - mx) * (y - my), b+w)
{-# INLINE covWeighted #-}

--weightedMeanCovMatrix :: G.Vector v Double => U.Vector Double -> M.Matrix Double -> (Vector Double, Matrix Double)
weightedMeanCovMatrix ws xs = (means, covMat)
  where
    means = G.generate n $ \i -> meanWeighted . G.zip ws . G.convert $ (xs `M.takeColumn` i)
    covMat = MM.create $ do
        mat <- MM.new (n,n)
        forM_ [0..n-1] $ \i ->
            forM_ [i..n-1] $ \j -> do
                let cov = covWeighted ws (xs `M.takeColumn` i, means G.! i) (xs `M.takeColumn` j, means G.! j) 
                MM.unsafeWrite mat (i,j) cov
                MM.unsafeWrite mat (j,i) cov
        return mat

    n = M.cols xs

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
        g <- initialize $ U.fromList [2,2,1,234,95,293,29992,232]
        kMeansInitial g ob k

test1 = do
    print $ train [[0.2,3.1,2],[2,1,3],[0.5,1.5,10],[0.6,2,0.2],[2,3,2],[2,3,6],[2,3,5]] 2 3

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
-- ([1, 0, 1, 0, 1, 1, 0, 1, 0, 1], -16.67738270170788)
test = do
   let hmm = GaussianHMM (U.fromList [0.5,0.5]) (M.fromLists [[0.1,0.9],[0.5,0.5]]) (V.fromList [m1,m2])
       m1 = mvn (vector [1]) (matrix 1 [1])
       m2 = mvn (vector [-1]) (matrix 1 [1])
       obs = V.fromList $ map (vector.return) [-1.6835, 0.0635, -2.1688, 0.3043, -0.3188, -0.7835, 1.0398, -1.3558, 1.0882, 0.4050]
   print $ viterbi hmm obs
   print $ forward hmm obs
