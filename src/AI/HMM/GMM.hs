{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE TypeFamilies #-}
module AI.HMM.GMM where

import Control.Arrow (first)
import Control.Monad (liftM, replicateM, forM_)
import Control.Monad.ST (runST)
import Data.Default.Class
import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Unboxed.Mutable as UM
import qualified Data.Matrix.Unboxed as MU
import System.Random.MWC

import AI.Function
import AI.HMM.Type
import AI.HMM.Class

import Debug.Trace

data GMMHMM = GMMHMM
    { _startProb :: !(U.Vector Double)
    , _transitionProb :: !(MU.Matrix Double)
    , _emission :: !(V.Vector MixMVN)
    , _mixtures :: !Int
    } deriving (Show)

bM :: GMMHMM -> MU.Matrix Double -> Int -> Int -> Int -> Double
bM h ob t i m = w + logPDF x o
  where
    MixMVN xs = _emission h `V.unsafeIndex` i
    (w, x) = xs `V.unsafeIndex` m
    o = ob `MU.takeRow` t

instance HMMLike GMMHMM (MU.Matrix Double) where
    data MStepOpt GMMHMM = MStepOpt
        { _covEstimator :: !CovEstimator
        }

    nSt = U.length . _startProb

    len _ = MU.rows

    isLogProb = const True

    π' h i = _startProb h `U.unsafeIndex` i

    a' h i j = _transitionProb h `MU.unsafeIndex` (i,j)

    b' h ob t i = logPDFMix (_emission h `V.unsafeIndex` i) o
      where
        o = ob `MU.takeRow` t

    setInitProb x h = h {_startProb=x}

    setTransProb x h = h {_transitionProb=MU.fromVector (nSt h,nSt h) x}

    updateEmission ob ps h opt = h {_emission=em'}
      where
        em' = case _covEstimator opt of
            FullCov -> V.generate r $ \i -> f . V.generate m $ \k ->
                let ws = U.slice ((i*m+k)*c) c ps
                    (mean, cov) = weightedMeanCovMatrix ws ob
                in (U.sum ws, mvn mean cov)
            DiagonalCov -> V.generate r $ \i -> f . V.generate m $ \k ->
                let ws = U.slice ((i*m+k)*c) c ps
                    (mean, cov) = weightedMeanDiagCov ws ob
                in (U.sum ws, mvnDiag mean cov)
            _ -> undefined
        f :: V.Vector (Double, MVN) -> MixMVN
        f xs = MixMVN $ V.map (first (log . (/s))) xs
          where s = V.foldl' (\acc x -> fst x + acc) 0 xs
        r = nSt h
        c = len h ob
        m = _mixtures h

    -- | forward probability
    forward' h ob = runST $ do
        mat <- UM.new (r*c*m)
        scales <- UM.new c

        -- update first column
        temp0 <- UM.new (r*m)
        forM_ [0..r*m-1] $ \i -> do
            let x = π' h (i `div` m) + bM h ob 0 (i `div` m) (i `mod` m)
            UM.unsafeWrite temp0 i x
        s0 <- fmap negate . logSumExpM $ temp0
        UM.unsafeWrite scales 0 s0
        -- normalize
        forM_ [0..r*m-1] $ \i -> UM.unsafeRead temp0 i >>= UM.unsafeWrite mat (i*c) . (+s0)

        -- update the rest of columns
        forM_ [1..c-1] $ \t -> do
            temp <- UM.new (r*m)
            forM_ [0..r*m-1] $ \j -> do
                temp' <- UM.new (r*m)
                forM_ [0..r*m-1] $ \i -> do
                    α_it' <- UM.unsafeRead mat $ i*c+t-1
                    UM.unsafeWrite temp' i $ α_it' + a' h (i `div` m) (j `div` m)
                    
                sum_α_a <- logSumExpM temp'

                UM.unsafeWrite temp j $ sum_α_a + bM h ob t (j `div` m) (j `mod` m)

            s <- fmap negate . logSumExpM $ temp
            UM.unsafeWrite scales t s
            -- normalize
            forM_ [0..r*m-1] $ \i -> UM.unsafeRead temp i >>= UM.unsafeWrite mat (i*c+t) . (+s)
        
        mat' <- U.unsafeFreeze mat
        scales' <- U.unsafeFreeze scales
        return (mat', scales')
      where
        r = nSt h
        c = len h ob
        m = _mixtures h
    {-# INLINE forward' #-}

{-
    backward' h ob scales = U.create $ do
        mat <- UM.new (r*c*m)
        -- fill in last column
        forM_ [0..r*m-1] $ \i -> UM.unsafeWrite mat (i*c+c-1) $ U.last scales + getWeight (_emission h `V.unsafeIndex` (i `div` m)) (i `mod` m)
        
        forM_ [c-2,c-3..0] $ \t -> do
            let sc = scales `U.unsafeIndex` t
            forM_ [0..r*m-1] $ \i -> do
                temp <- UM.new (r*m)
                forM_ [0..r*m-1] $ \j -> do
                    let b_jmo = bM h ob (t+1) (j `div` m) (j `mod` m)
                        a_ij = a' h (i `div` m) (j `div` m)
                    β_jmt' <- UM.unsafeRead mat (j*c+t+1)
                    UM.unsafeWrite temp j $! b_jmo + a_ij + β_jmt'
                x <- logSumExpM temp
                UM.unsafeWrite mat (i*c+t) $! x + sc
        return mat
      where
        r = nSt h
        c = len h ob
        m = _mixtures h
    {-# INLINE backward' #-}
    -}

    baumWelch' opt ob h = (setInitProb ini' . setTransProb trans' . updateEmission ob (U.map exp γ') h $ opt, scales)
      where
        ini' = U.generate r $ \i -> γ `U.unsafeIndex` (i*c)

        trans' = U.create $ do
            mat <- UM.new (r*r)
            forM_ [0..r-1] $ \i ->
                forM_ [0..r-1] $ \j -> do
                    let a_ij = a' h i j
                    temp <- UM.new $ c - 1
                    forM_ [1 .. c-1] $ \t -> do
                        let b_jo = b' h ob t j
                            α_it' = fw `U.unsafeIndex` (i*c+t-1)
                            β_jt = bw `U.unsafeIndex` (j*c+t)
                        UM.unsafeWrite temp (t-1) $ b_jo + α_it' + β_jt
                    x <- logSumExpM temp
                    UM.unsafeWrite mat (i*r+j) $ a_ij + x
            forM_ [0..r-1] $ \i -> do
                s <- logSumExpM $ UM.slice (i*r) r mat
                forM_ [0..r-1] $ \j -> UM.unsafeRead mat (i*r+j) >>= UM.unsafeWrite mat (i*r+j) . subtract s
            return mat

        (fw', scales) = forward' h ob
        fw = U.generate (r*c) $ \i -> logSumExp $ U.generate m $ \x -> fw' `U.unsafeIndex` (((i `div` c) * m + x) * c +  (i `mod` c))
        bw = backward' h ob scales
        γ' = U.generate (r*c*m) $ \i -> fw' `U.unsafeIndex` i + bw `U.unsafeIndex` (i `div` c `div` m * c + i `mod` c) - scales `U.unsafeIndex` (i `mod` c)
        γ = U.generate (r*c) $ \i -> logSumExp $ U.generate m $ \x -> γ' `U.unsafeIndex` (((i `div` c) * m + x) * c +  (i `mod` c))
        r = nSt h
        c = len h ob
        m = _mixtures h
    {-# INLINE baumWelch' #-}

    randHMM g opt = do
        let randProbVector = liftM normalize . uniformVector g
            normalize xs = U.map log $ U.map (/ U.sum xs) xs
            randMVN = do
                mean <- uniformVector g p
                d <- uniformVector g p
                return $ mvn mean $ MU.flatten $ MU.diag (d :: V.Vector Double)
            randMVNMix = do mvns <- replicateM m randMVN
                            ws <- randProbVector m
                            return $ MixMVN $ V.zip (V.convert ws) $ V.fromList mvns
        startProb <- randProbVector s
        transProb <- liftM MU.fromRows $ replicateM s $ randProbVector s
        em <- liftM V.fromList $ replicateM s randMVNMix
        return $ GMMHMM startProb transProb em m
      where
        s = _nStates opt
        p = _nObs opt
        m = _nMix opt
        

instance Default (HMMOpt GMMHMM) where
    def = HMMOpt
        { _seed = toSeed $ U.fromList [22]
        , _initialization = Random
        , _nIter = 50
        , _nStates = 2
        , _nObs = 1
        , _nMix = 2
        , _mStepOpt = MStepOpt DiagonalCov
        }

test :: IO ()
test = do
    let obs = MU.fromLists $ map return [ -1.6835, 0.0635, -2.1688, 0.3043, -0.3188
                                    , -0.7835, 1.0398, -1.3558, 1.0882, 0.4050 ]
        h = fitHMM obs def :: GMMHMM
    print h
