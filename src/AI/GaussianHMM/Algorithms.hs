{-# LANGUAGE BangPatterns #-}

module AI.GaussianHMM.Algorithms
   ( baumWelch
   , forward
   , backward
   , viterbi
   , kMeansInitial
   ) where

import Control.Monad.Primitive 
import Control.Monad (forM_)
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
import Data.List (groupBy, sortBy, maximumBy, foldl')
import Data.Ord
import Data.Function
import System.Random.MWC

import AI.Function
import AI.MVN
import AI.GaussianHMM.Types
import AI.Clustering.KMeans

import Debug.Trace

baumWelch :: Observation -> GaussianHMM -> (GaussianHMM, U.Vector Double)
baumWelch ob h = traceShow em' $ (GaussianHMM ini' trans' em', scales)
  where
    ini' = γ `M.takeColumn` 0

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

    em' = G.generate n $ \i -> let ws = G.map exp $ γ `M.takeRow` i
                                   (mean, cov) = weightedMeanCovMatrix ws ob
                               in diagCov mean cov

    (fw, scales) = forward h ob
    bw = backward h ob scales
    γ :: MU.Matrix Double
    γ = MU.generate (n,m) $ \(s,t) -> fw `M.unsafeIndex` (s,t) + bw `M.unsafeIndex` (s,t) - scales `G.unsafeIndex` t
    n = nSt h
    m = M.rows ob
    f :: MU.Matrix Double -> MU.Matrix Double
    f x = runST $ do
        x' <- MM.unsafeThaw $ MU.tr $ MU.map (\i -> if i < -1e200 then -1e200 else i) x
        normalizeRow x'
        x'' <- MM.unsafeFreeze x'
        return $ MU.tr x''

diagCov m cov = mvn m $ convert $ M.diag (M.takeDiag cov :: V.Vector Double)

mvnEstimate m cov = mvn m $ reshape (M.rows cov) $ fst $ glasso (M.rows cov) (M.flatten cov) 0.01

fullCov m cov = MVN m (convert cov) invcov logdet
  where
    (invcov, (logdet, _)) = invlndet $ convert cov


{-
-- | the E step in EM algorithm
eStep :: GaussianHMM -> Observation -> (
eStep h ob = 
  where
    (fw, sc) = forward h ob
    bw = backward h ob sc
    γ = G.zipWith (-) (G.zipWith (+) fw bw) sc
    ξ = 
-}

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
            temp' <- UM.new r
            forM_ [0..r-1] $ \i -> do
                α_it' <- MM.read mat (i,t-1)
                GM.unsafeWrite temp' i $ α_it' + a h (i,j)
                
            sum_α_a <- logSumExpM temp'

            GM.unsafeWrite temp j $ sum_α_a + b h j (ob `M.takeRow` t)

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

viterbi :: GaussianHMM -> Observation -> ([Int], Double)
viterbi h ob = foldl track ([],undefined) . init . foldl' f [] . M.toRows $ ob
  where
    f [] o = [U.generate n $ \i -> (π h i + b h i o, -1)]
    f acc@(l_t:_) o =
        let score i j = fst (G.unsafeIndex l_t i) + a h (i,j)
            vec = U.generate n $ \j -> let (i, x) = maximumBy (comparing snd) . map (\i' -> (i', score i' j)) $ [0..n-1]
                                       in (x + b h j o, i)
        in vec : acc
    track ([], _) v = let ((p, prev), i) = G.maximumBy (comparing (fst.fst)) $ G.zip v $ G.enumFromN 0 n
                      in ([prev, i], p)
    track (path@(i:_), p) v = (snd (G.unsafeIndex v i) : path, p)
    n = nSt h

{-
randomInitial :: PrimMonad m
              => Gen (PrimState m)
              -> Observation
              -> Int
              -> m GaussianHMM
randomInitial g ob k = do
    vec <- uniformVector g $ n * k
  where
    n = M.rows ob
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

meanCov dat = (meanVec, reshape p $ fst $ glasso p (M.flatten covs) 0.01)
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
covWithMean (mx, xs) (my, ys) | n == 1 = 1e-200
                              | otherwise = G.sum (G.zipWith f xs ys) / (n - 1)
  where
    f x y = (x - mx) * (y - my)
    n = fromIntegral $ G.length xs

convert mat = reshape c $ M.flatten mat
  where
    c = M.cols mat

hmmExample :: (GaussianHMM, Observation)
hmmExample = (hmm, obs2)
  where
    hmm = GaussianHMM (U.fromList $ map log [0.5,0.5])
                      (M.fromLists $ ((map.map) log) [[0.1,0.9],[0.5,0.5]])
                      (V.fromList [m1,m2])
    m1 = mvn (vector [1]) (matrix 1 [1])
    m2 = mvn (vector [-1]) (matrix 1 [1])
    obs :: Observation
    obs = M.fromLists $ map return [ -1.6835, 0.0635, -2.1688, 0.3043, -0.3188
                                   , -0.7835, 1.0398, -1.3558, 1.0882, 0.4050 ]

    obs2 = M.fromLists $ map return [ -0.2370699812741054, -0.7650421248531205, -1.2449364749699783
                                    , -0.4940448810044036, -1.0751810551549668, -0.573496578580446
                                    , -0.7244465394993085, -0.44320733687218794, -0.5137119768758072]

test = do
    let (hmm, obs) = hmmExample
    loop obs hmm 0
  where
    loop o h i | i > 100 = print h
               | otherwise = let h' = fst $ baumWelch o h
                                 (f, sc) = forward h o
                             in traceShow (G.sum sc) $ loop o h' (i+1)
