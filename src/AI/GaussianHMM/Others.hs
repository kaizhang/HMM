{-# LANGUAGE BangPatterns #-}

module AI.GaussianHMM.Algorithms where

import Control.Monad.Primitive 
import Control.Monad (forM_, foldM)
import Control.Monad.ST (runST)
import Numeric.LinearAlgebra.HMatrix (Matrix, Vector, (<>), invlndet, (!), tr, asRow, vector, matrix, reshape)
import qualified Data.Vector as V
import qualified Data.Vector.Mutable as VM
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Generic.Mutable as GM
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Unboxed.Mutable as UM
import qualified Data.Matrix.Unboxed as MU
import qualified Data.Matrix.Generic as M
import qualified Data.Matrix.Generic.Mutable as MM
import qualified Data.Matrix.Storable as S
import Algorithms.GLasso (glasso)
import Statistics.Sample (mean)
import Data.List (groupBy, sortBy)
import Data.Ord
import Data.Function
import System.Random.MWC

import AI.Function
import AI.GaussianHMM.Types
import AI.Clustering.KMeans

import Debug.Trace

-------------------------------------------------------------------------------

baumWelch :: Observation -> GaussianHMM -> (GaussianHMM, U.Vector Double)
baumWelch ob h = (GaussianHMM ini' trans' em', scales)
  where
    ini' = G.zipWith (\x y -> x + y - G.head scales) (fw `M.takeColumn` 0) (bw `M.takeColumn` 0)

    trans' = MM.create $ do
        mat <- MM.new (n,n)
        forM_ [0..n-1] $ \i ->
            forM_ [0..n-1] $ \j -> do
                let a_ij = a h (i,j)
                temp <- UM.new m
                forM_ [1 .. m - 1] $ \t -> do
                    let o = ob `M.takeRow` t
                        b_jo = b h j o
                        α_it' = fw `M.unsafeIndex` (i,t-1)
                        β_jt = bw `M.unsafeIndex` (j,t)
                    GM.unsafeWrite temp t $ b_jo + α_it' + β_jt
                x <- logSumExpM temp
                MM.unsafeWrite mat (i,j) $ a_ij + x
        normalizeRow mat
        return mat

    em' = G.generate n $ \i -> let ws = G.map f $ G.enumFromN 0 m
                                   f t = let α_it = fw `M.unsafeIndex` (i,t)
                                             β_it = bw `M.unsafeIndex` (i,t)
                                             sc = scales `G.unsafeIndex` t
                                          in exp $ α_it + β_it - sc
                                   (mean, cov) = traceShow (G.sum ws) $ weightedMeanCovMatrix ws ob
                               in mvn mean (convert $ fst $ glasso cov 0.1)

    (fw, scales) = forward h ob
    bw = backward h ob scales
    n = nSt h
    m = M.rows ob

forward :: GaussianHMM -> Observation -> (MU.Matrix Double, U.Vector Double)
forward h ob = runST $ do
    mat <- MM.new (r,c)
    scales <- GM.new c

    -- update first column
    temp0 <- UM.new r
    forM_ [0..r-1] $ \i -> do
        let x = π h i + b h i (M.takeRow ob 0)
        GM.unsafeWrite temp0 i x
    s0 <- fmap negate . logSumExpM $ temp0
    GM.unsafeWrite scales 0 s0
    -- normalize
    forM_ [0..r-1] $ \i -> GM.unsafeRead temp0 i >>= MM.unsafeWrite mat (i,0) . (+s0)

    -- update the rest of columns
    forM_ [1..c-1] $ \t -> do
        temp <- UM.new r
        forM_ [0..r-1] $ \j -> do
            sum_α_a <- foldM ( \acc i -> do
                α_it' <- MM.unsafeRead mat (i,t-1)
                return $! acc + exp (α_it' + a h (i,j)) ) 0 [0..r-1]
            GM.unsafeWrite temp j $ log sum_α_a + b h j (ob `M.takeRow` t)

        s <- fmap negate . logSumExpM $ temp
        GM.unsafeWrite scales t s
        -- normalize
        forM_ [0..r-1] $ \i -> GM.unsafeRead temp i >>= MM.unsafeWrite mat (i,t) . (+s)
    
    mat' <- MM.unsafeFreeze mat
    scales' <- G.unsafeFreeze scales
    return (mat', scales')
  where
    r = nSt h
    c = M.rows ob
{-# INLINE forward #-}

-- | backward in log scale
backward :: GaussianHMM -> Observation -> U.Vector Double -> MU.Matrix Double
backward h ob scales = MM.create $ do
    mat <- MM.new (r,c)
    -- fill in last column
    forM_ [0..r-1] $ \i -> MM.unsafeWrite mat (i,c-1) $ G.last scales
    
    forM_ [c-2,c-3..0] $ \t -> do
        let sc = scales `G.unsafeIndex` t
        forM_ [0..r-1] $ \i -> do
            temp <- UM.new r
            forM_ [0..r-1] $ \j -> do
                let b_jo = b h j $ ob `M.takeRow` (t+1)
                    a_ij = a h (i,j)
                β_jt' <- MM.unsafeRead mat (j,t+1)
                UM.unsafeWrite temp j $! b_jo + a_ij + β_jt'
            x <- logSumExpM temp
            MM.unsafeWrite mat (i,t) $! x + sc
    return mat
  where
    r = nSt h
    c = M.rows ob
{-# INLINE backward #-}


-------------------------------------------------------------------------------

update :: Observation -> GaussianHMM -> GaussianHMM
update ob h = GaussianHMM ini' trans' em'
  where
    ini' = p_St `M.takeColumn` 0
    trans' = MM.create $ do
        mat <- MM.new (m,m)
        forM_ [0..m-1] $ \i ->
            forM_ [0..m-1] $ \j -> MM.unsafeWrite mat (i,j) $ 
                logSumExp $ U.generate n $ \t -> (q G.! t) M.! (i,j)

        normalizeRow mat
        return mat

    em' = G.generate m $ \i -> let ws = G.map exp $ p_St `M.takeRow` i
                                   (mean, cov) = weightedMeanCovMatrix ws ob
                               in mvn mean (convert $ fst $ glasso cov 0.1)

    q = recursion ob h
    p_St = forwardRecursion q
    m = nSt h
    n = M.rows ob

convert mat = reshape c $ M.flatten mat
  where
    c = M.cols mat

-- | compute q(s_t | s_{t-1}, O)
recursion :: Observation -> GaussianHMM -> V.Vector (MU.Matrix Double)
recursion ob hmm = runST $ do
    mats <- V.replicateM n $ MM.new (m,m)

    -- fill in last matrix m(i,j) = p(j|i)
    let lst = V.last mats
    forM_ [0..m-1] $ \i ->
        forM_ [0..m-1] $ \j -> do
            MM.unsafeWrite lst (i,j) $ a hmm (i,j) + b hmm j (ob `M.takeRow` (n-1))
    -- normalize
    normalizeRow lst

    -- fill in rest
    forM_ [n-2,n-3..1] $ \t -> do
        let o = ob `M.takeRow` t
            q_t = V.unsafeIndex mats t
        forM_ [0..m-1] $ \u_t' ->
            forM_ [0..m-1] $ \u_t -> do
                tmp <- UM.new m
                forM_ [0..m-1] $ \u_t'' -> do
                    q_t' <- V.unsafeIndex mats (t+1) `MM.unsafeRead` (u_t,u_t'') -- q(S_t+1|S_t,O)
                    
                    -- q(s_t|s_t-1,s_t+1,O_t)
                    tmp' <- UM.new m
                    forM_ [0..m-1] $ \i -> GM.unsafeWrite tmp' i $ b hmm i o + a hmm (u_t',i) + a hmm (i, u_t'')
                    s' <- logSumExpM tmp'
                    qq <- fmap (subtract s') $ GM.unsafeRead tmp' u_t

                    GM.unsafeWrite tmp u_t'' $ q_t' - qq
                s <- logSumExpM tmp
                MM.unsafeWrite q_t (u_t',u_t) (-s)

    -- fill in first matrix
    let first = V.head mats
    forM_ [0..m-1] $ \u_1 -> do
        tmp <- UM.new m
        forM_ [0..m-1] $ \u_2 -> do
            q_t' <- V.unsafeIndex mats 1 `MM.unsafeRead` (u_1,u_2) -- q(S_2|S_1,O)

            -- q(s_1|s_2,O)
            tmp' <- UM.new m
            forM_ [0..m-1] $ \i -> GM.unsafeWrite tmp' i $ b hmm i (ob `M.takeRow` 0) + a hmm (i,u_2)
            s' <- logSumExpM tmp'
            qq <- fmap (subtract s') $ GM.unsafeRead tmp' u_1
            GM.unsafeWrite tmp u_2 $ q_t' - qq
        s <- logSumExpM tmp
        MM.unsafeWrite first (0,u_1) (-s)

    G.mapM MM.unsafeFreeze mats
  where
    m = nSt hmm
    n = M.rows ob

forwardRecursion :: V.Vector (MU.Matrix Double) -> MU.Matrix Double
forwardRecursion mats = MM.create $ do
    mat <- MM.new (m, n)
    -- fill in first column
    let first = G.head mats
    forM_ [0..m-1] $ \i -> MM.unsafeWrite mat (i,0) $ first M.! (0,i)

    -- fill in rest
    forM_ [1..n-1] $ \t -> do
        let q_t = mats G.! t
        forM_ [0..m-1] $ \u_t -> do
            tmp <- UM.new m
            forM_ [0..m-1] $ \u_t' -> do
                p_u_t' <- MM.unsafeRead mat (u_t', t-1)
                UM.unsafeWrite tmp u_t' $ (q_t M.! (u_t', u_t)) + p_u_t'
            s <- logSumExpM tmp
            MM.unsafeWrite mat (u_t,t) s
    return mat
  where
    m = M.rows $ G.head mats
    n = G.length mats

-------------------------------------------------------------------------------

-------------------------------------------------------------------------------
--
--
{-
randomInitial :: PrimMonad m
              => Gen (PrimState m)
              -> Observation
              -> Int
              -> m GaussianHMM
randomInitial g k = do
    uniform
    -}

-- | construct inital HMM model by kmeans clustering
kMeansInitial :: PrimMonad m
              => Gen (PrimState m)
              -> Observation
              -> Int
              -> m GaussianHMM
kMeansInitial g ob k = do 
    (membership, _) <- kmeans g k ob
    let pi = G.map (log . (/ fromIntegral n)) $ G.create $ do
            vec <- GM.replicate k 0.0001
            V.forM_ membership $ \s -> GM.unsafeRead vec s >>= GM.unsafeWrite vec s . (+1)
            return vec

        trans = M.map log $ MM.create $ do
            mat <- MM.replicate (k,k) 0.0001
            V.sequence_ . V.zipWith ( \i j -> MM.unsafeRead mat (i,j) >>=
                MM.unsafeWrite mat (i,j) . (+1) ) membership . V.tail $
                membership
            normalizeByRow k mat
            return mat
        
        emisson = G.fromList $ map (uncurry mvn . meanCov) clusters
          where
            clusters = map (M.fromRows . snd . unzip) . groupBy ((==) `on` fst) . sortBy (comparing fst) $ zip (G.toList membership) $ M.toRows ob
    return $ GaussianHMM pi trans emisson
  where
    n = M.rows ob
    normalizeByRow x mat = forM_ [0..x-1] $ \i -> do
        sc <- G.foldM' (\acc j -> fmap (+acc) $ MM.unsafeRead mat (i,j)) 0 $ U.enumFromN 0 x
        forM_ [0..x-1] $ \j -> MM.unsafeRead mat (i,j) >>= MM.unsafeWrite mat (i,j) . (/sc)

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
