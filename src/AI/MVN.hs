module AI.MVN where

import Data.Binary
import qualified Data.Vector.Generic as G
import Numeric.LinearAlgebra.HMatrix (Matrix, Vector, (<>), invlndet, (!), tr, asRow, toLists, fromLists)

-- | multivariate normal distribution
data MVN = MVN
    { _mean :: !(Vector Double)
    , _cov :: !(Matrix Double)
    , _invcov :: !(Matrix Double)
    , _logdet :: !Double  -- ^ log determinant of covariance matrix
    } deriving (Show)

instance Binary MVN where
    put (MVN a b c d) = do
        put $ G.toList a
        put $ toLists b
        put $ toLists c
        put d
    get = do
        a <- get
        b <- get
        c <- get
        d <- get
        return $ MVN (G.fromList a) (fromLists b) (fromLists c) d

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
