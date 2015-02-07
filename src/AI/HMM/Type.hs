{-# LANGUAGE BangPatterns #-}
module AI.HMM.Type where

import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Generic as G
import Numeric.LinearAlgebra.HMatrix

-- | multivariate normal distribution
data MVN = MVN
    { _mean :: !(U.Vector Double)
    , _cov :: !(U.Vector Double)  -- ^ row majored covariance matrix
    , _invcov :: !(U.Vector Double)
    , _logdet :: !Double  -- ^ log determinant of covariance matrix
    , _dim :: !Int
    } deriving (Show)

mvn :: U.Vector Double -> U.Vector Double -> MVN
mvn m cov | d*d /= U.length cov = error "incompatible dimemsion of mean and covariance"
          | otherwise = MVN m cov (G.convert $ flatten invcov) logdet d
  where
    (invcov, (logdet, _)) = invlndet $ reshape d $ G.convert cov
    d = U.length m
{-# INLINE mvn #-}

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
logPDF (MVN m _ invcov logdet d) x = -0.5 * (fromIntegral d * log (2*pi) + logdet + quadTerm)
  where
    quadTerm = loop 0 0
      where
        x' = U.zipWith (-) x m
        loop !acc !i
            | i < d*d = let r = i `div` d
                            c = i `mod` d
                        in acc + U.unsafeIndex invcov i * U.unsafeIndex x' r * U.unsafeIndex x' c
            | otherwise = acc
{-# INLINE logPDF #-}
