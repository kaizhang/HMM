{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE TypeFamilies #-}
module AI.HMM.Gaussian
    ( GaussianHMM(..)
    , module AI.HMM.Class
    ) where

import Algorithms.GLasso (glasso)
import Control.Lens
import Control.Monad (liftM, replicateM)
import Data.Default.Class
import qualified Data.Vector as V
import qualified Data.Vector.Generic as G
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
        , _lambda :: !Double  -- ^ glasso parameter
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
        em' = V.generate r $ \i -> estimateMVN ob $ U.slice (i*c) c ps
        estimateMVN o ws = case _covEstimator opt of
            FullCov -> let (mean, cov) = weightedMeanCovMatrix ws o
                       in mvn mean cov
            DiagonalCov -> let (mean, cov) = weightedMeanDiagCov ws o
                           in mvnDiag mean cov
            LassoCov -> let (mean, cov) = weightedMeanCovMatrix ws o
                            (cov', _) = glasso (G.length mean) cov $ _lambda opt
                        in mvn (U.convert mean) $ U.convert cov'
        r = nSt h
        c = len h ob
    {-# INLINE updateEmission #-}

    randHMM g opt = do
        let randProbVector = liftM normalize . uniformVector g
            normalize xs = U.map log $ U.map (/ U.sum xs) xs
            randMVN = do
                m <- uniformVector g p
                d <- uniformVector g p
                return $ mvn m $ MU.flatten $ MU.diag (d :: V.Vector Double)
        startProb <- randProbVector s
        transProb <- liftM MU.fromRows $ replicateM s $ randProbVector s
        em <- liftM V.fromList $ replicateM s randMVN
        return $ GaussianHMM startProb transProb em
      where
        s = _nStates opt
        p = _nObs opt

instance Default (HMMOpt GaussianHMM) where
    def = HMMOpt
        { _seed = toSeed $ U.fromList [22]
        , _initialization = Random
        , _nIter = 50
        , _nStates = 2
        , _nObs = 1
        , _nMix = 1
        , _mStepOpt = MStepOpt FullCov 0.01
        }
