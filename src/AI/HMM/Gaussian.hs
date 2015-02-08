{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE TypeFamilies #-}
module AI.HMM.Gaussian where

import Control.Monad (liftM, replicateM)
import Data.Default.Class
import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as U
import qualified Data.Matrix.Unboxed as MU
import System.Random.MWC

import AI.Function
import AI.HMM.Type
import AI.HMM.Class

data GaussianHMM = GaussianHMM
    { _startProb :: !(U.Vector Double)
    , _transitionProb :: !(MU.Matrix Double)
    , _emission :: !(V.Vector MVN)
    } deriving (Show)

instance HMMLike GaussianHMM (MU.Matrix Double) where
    data MStepOpt GaussianHMM = MStepOpt
        { _covEstimator :: !CovEstimator
        }

    nSt = U.length . _startProb

    len _ = MU.rows

    isLogProb = const True

    π' h i = _startProb h `U.unsafeIndex` i

    a' h i j = _transitionProb h `MU.unsafeIndex` (i,j)

    b' h ob t i = logPDF (_emission h `V.unsafeIndex` i) o
      where
        o = ob `MU.takeRow` t

    setInitProb x h = h {_startProb=x}

    setTransProb x h = h {_transitionProb=MU.fromVector (nSt h,nSt h) x}

    updateEmission ob ps h opt = h {_emission=em'}
      where
        em' = case _covEstimator opt of
            FullCov -> V.generate r $ \i ->
                let ws = U.slice (i*c) c ps
                    (mean, cov) = weightedMeanCovMatrix ws ob
                in mvn mean cov
            DiagonalCov -> V.generate r $ \i ->
                let ws = U.slice (i*c) c ps
                    (mean, cov) = weightedMeanCovMatrix ws ob
                    n = U.length mean
                in mvn mean $ U.imap (\x v -> if x `div` n /= x `mod` n then 0 else v) cov
            _ -> undefined
        r = nSt h
        c = len h ob

    randHMM g n p = do
        let randProbVector = liftM normalize . uniformVector g
            normalize xs = U.map log $ U.map (/ U.sum xs) xs
            randMVN = do
                m <- uniformVector g p
                d <- uniformVector g p
                return $ mvn m $ MU.flatten $ MU.diag (d :: V.Vector Double)
        startProb <- randProbVector n
        transProb <- liftM MU.fromRows $ replicateM n $ randProbVector n
        em <- liftM V.fromList $ replicateM n randMVN
        return $ GaussianHMM startProb transProb em

instance Default (HMMOpt GaussianHMM) where
    def = HMMOpt
        { _seed = toSeed $ U.fromList [22]
        , _initialization = Random
        , _nIter = 50
        , _nStates = 2
        , _nObs = 1
        , _nMixture = 1
        , _mStepOpt = MStepOpt FullCov
        }

example :: (GaussianHMM, MU.Matrix Double)
example = (hmm, obs)
  where
    hmm = GaussianHMM (U.fromList $ map log [0.5,0.5])
                      (MU.fromLists $ (map.map) log [[0.1,0.9],[0.5,0.5]])
                      (V.fromList [m1,m2])
    m1 = mvn (U.fromList [1]) (U.fromList [1])
    m2 = mvn (U.fromList [-1]) (U.fromList [1])
    obs = MU.fromLists $ map return [ -1.6835, 0.0635, -2.1688, 0.3043, -0.3188
                                   , -0.7835, 1.0398, -1.3558, 1.0882, 0.4050 ]

test :: IO ()
test = do
    let h = fitHMM (snd example) def {_initialization=Fixed (fst example)}:: GaussianHMM
    print h
