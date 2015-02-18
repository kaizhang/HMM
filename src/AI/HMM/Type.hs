{-# LANGUAGE BangPatterns #-}
module AI.HMM.Type where

import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Generic as G
import Numeric.LinearAlgebra.HMatrix

import AI.Function

-- constant: log (2*pi)
m_log_2_pi :: Double
m_log_2_pi = 1.8378770664093453

-- | multivariate normal distribution
data MVN = MVN
    { _mean :: !(U.Vector Double)
    , _cov :: !(U.Vector Double)  -- ^ row majored covariance matrix
    , _invcov :: !(U.Vector Double)
    , _logdet :: !Double  -- ^ log determinant of covariance matrix
    , _dim :: !Int
    , _covType :: !CovType
    } deriving (Show)

data CovType = Full
             | Diagonal
    deriving (Show)

data CovEstimator = FullCov
                  | DiagonalCov
                  | LassoCov

mvn :: U.Vector Double -> U.Vector Double -> MVN
mvn m cov | d*d /= U.length cov = error "incompatible dimemsion of mean and covariance"
          | otherwise = MVN m cov (G.convert $ flatten invcov) logdet d Full
  where
    (invcov, (logdet, _)) = invlndet $ reshape d $ G.convert cov
    d = U.length m
{-# INLINE mvn #-}

mvnDiag :: U.Vector Double -> U.Vector Double -> MVN
mvnDiag m cov | d /= U.length cov = error "incompatible dimemsion of mean and covariance"
              | otherwise = MVN m cov invcov logdet d Diagonal
  where
    invcov = U.map (1/) cov
    logdet = U.sum . U.map log $ cov
    d = U.length m
{-# INLINE mvnDiag #-}

{-
-- | log probability of MVN
logPDF :: MVN -> Vector Double -> Double
logPDF (MVN m _ invcov logdet) x = -0.5 * ( d * log (2*pi) + logdet
                                        + (x' <> invcov <> tr x') ! 0 ! 0
                                        )
  where
    x' = asRow $ x - m
    d = fromIntegral . G.length $ m
{-# INLINE logPDF #-}
-}

logPDF :: MVN -> U.Vector Double -> Double
logPDF (MVN m _ invcov logdet d t) x = -0.5 * (fromIntegral d * m_log_2_pi + logdet + quadTerm)
  where
    quadTerm = case t of
        Full -> loop 0 0
        Diagonal -> U.sum $ U.zipWith (*) invcov $ U.map (**2) x'
      where
        loop !acc !i
            | i < d*d = let r = i `div` d
                            c = i `mod` d
                        in acc + U.unsafeIndex invcov i * U.unsafeIndex x' r * U.unsafeIndex x' c
            | otherwise = acc
        x' = U.zipWith (-) x m
{-# INLINE logPDF #-}


-- | mixture of multivariate normal variables, weight is in log form
newtype MixMVN = MixMVN (V.Vector (Double, MVN)) deriving (Show)

getWeight :: MixMVN -> Int -> Double
getWeight (MixMVN xs) i = fst $ xs V.! i
{-# INLINE getWeight #-}

logPDFMix :: MixMVN -> U.Vector Double -> Double
logPDFMix (MixMVN v) xs = logSumExp . V.map (\(logw, m) -> logw + logPDF m xs) $ v
{-# INLINE logPDFMix #-}
