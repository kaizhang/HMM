{-# LANGUAGE FlexibleContexts #-}
module AI.Function where

import Control.Monad.Primitive
import Numeric.LinearAlgebra.HMatrix (Matrix, Vector, (<>), invlndet, (!), tr, asRow)
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Generic.Mutable as GM

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

logSumExp :: G.Vector v Double => v Double -> Double
logSumExp xs = m + log (G.foldl' (\acc x -> acc + exp (x-m)) 0 xs)
  where
    m = G.maximum xs
{-# INLINE logSumExp #-}

logSumExpM :: (PrimMonad m, GM.MVector v Double) => v (PrimState m) Double -> m Double
logSumExpM xs = do
    m <- U.foldM' (\acc i -> return . max acc =<< GM.unsafeRead xs i) (-1/0) . U.enumFromN 0 $ n
    let f acc i = do
            x <- GM.unsafeRead xs i
            return $ acc + exp (x-m)
    s <- U.foldM' f 0 . U.enumFromN 0 $ n
    return $! m + log s
  where
    n = GM.length xs
{-# INLINE logSumExpM #-}
