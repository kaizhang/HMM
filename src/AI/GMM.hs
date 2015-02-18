{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE TemplateHaskell #-}
module AI.GMM where

import Control.Arrow (first)
import Control.Lens (makeLenses)
import Control.Monad
import Control.Monad.ST
import Data.Default.Class
import qualified Data.Matrix.Unboxed as MU
import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as U
import System.Random.MWC

import AI.Function
import AI.HMM.Type

import Debug.Trace

data Initialization = Fixed MixMVN
                    | Random

data GMMOpt = GMMOpt
    { _covEstimator :: CovEstimator
    , _nIter :: Int
    , _nMix :: Int
    , _initialization :: Initialization
    , _seed :: Seed
    }

makeLenses ''GMMOpt

instance Default GMMOpt where
    def = GMMOpt
        { _covEstimator = DiagonalCov
        , _nIter = 20
        , _nMix = 2
        , _initialization = Random
        , _seed = toSeed $ U.fromList [22]
        }

fitGMM :: MU.Matrix Double -> GMMOpt -> MixMVN
fitGMM xs opt = loop (_nIter opt) initM
  where
    loop !iter !m | iter > 0 = traceShow (sum . map (logPDFMix m) $ MU.toRows xs) loop (iter-1) . mStep . eStep $ m
                  | otherwise = m 
    -- k x n
    eStep :: MixMVN -> U.Vector Double
    eStep (MixMVN mvns) = normalize $ U.generate (s*n) $ \i ->
        let k = i `div` n
            x = xs `MU.takeRow` (i `mod` n)
            (w, mvn) = mvns V.! k
        in w + logPDF mvn x
      where
        normalize v =
            let sums = U.generate n $ \j -> logSumExp . U.map (\i -> v U.! (i*n+j)) . U.enumFromN 0 $ s
            in U.imap (\i x -> exp $ x - (sums U.! (i `mod` n))) v

    mStep :: U.Vector Double -> MixMVN
    mStep ws = case _covEstimator opt of
        FullCov -> f . V.generate s $ \k ->
            let w = U.slice (k*n) n ws
                (mean, cov) = weightedMeanCovMatrix w xs
            in (U.sum w, mvn mean cov)

        DiagonalCov -> f . V.generate s $ \k ->
                let w = U.slice (k*n) n ws
                    (mean, cov) = weightedMeanDiagCov w xs
                in (U.sum w, mvnDiag mean cov)
        _ -> undefined

    initM = case _initialization opt of
        Fixed m -> m
        Random -> runST $ do
            g <- restore $ _seed opt
            v <- uniformVector g (s*n)
            let sums = U.generate n $ \j -> U.sum . U.map (\i -> v U.! (i*n+j)) . U.enumFromN 0 $ s
                v' = U.imap (\i x -> x / (sums U.! (i `mod` n))) v
            return $ mStep v'

    f :: V.Vector (Double, MVN) -> MixMVN
    f xs = MixMVN $ V.map (first (log . (/t))) xs
      where t = V.foldl' (\acc x -> fst x + acc) 0 xs
    s = _nMix opt
    n = MU.rows xs
    p = MU.cols xs

test :: IO ()
test = do
    let -- obs = MU.fromLists $ map return [ -1.6835, 0.0635, -2.1688, 0.3043, -0.3188
        --                            , -0.7835, 1.0398, -1.3558, 1.0882, 0.4050, 1,2,-0.5 ]
        obs = MU.fromLists $ [ [ -31,-2]
                             , [3.2,0.5]
                             , [100,10]
                             , [-13.2,-100]
                             , [-40,-4]
                             , [0,100]
                             ]
        h = fitGMM obs def
    print h
