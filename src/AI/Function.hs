{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE BangPatterns #-}
module AI.Function where

import Control.Monad (forM_)
import Control.Monad.Primitive
import Numeric.LinearAlgebra.HMatrix (Matrix, Vector, (<>), invlndet, (!), tr, asRow)
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Generic.Mutable as GM
import qualified Data.Matrix.Generic as M
import qualified Data.Matrix.Generic.Mutable as MM
import Statistics.Sample (meanWeighted)

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

-- | normalize a matrix of which elements are in log form
normalizeRow :: (PrimMonad m, GM.MVector v Double) => M.MMatrix v (PrimState m) Double -> m ()
normalizeRow mat =
    forM_ [0..r-1] $ \i -> do
        let vec = mat `MM.takeRow` i
        x <- logSumExpM vec
        forM_ [0..c-1] $ \j ->
            MM.unsafeRead mat (i,j) >>= MM.unsafeWrite mat (i,j) . subtract x
  where
    (r,c) = MM.dim mat

covWeighted :: U.Vector Double -> (Vector Double, Double) -> (Vector Double, Double) -> Double
covWeighted ws (xs, mx) (ys, my) = g . G.foldl f (0,0) $ U.enumFromN 0 $ G.length ws
  where
    f (!a,!b) i = let w = ws G.! i
                      x = xs G.! i
                      y = ys G.! i
                  in (a + w * (x - mx) * (y - my), b+w)
    g (a,b) | b == 0 = 0
            | otherwise = a / b
{-# INLINE covWeighted #-}

--weightedMeanCovMatrix :: G.Vector v Double => U.Vector Double -> M.Matrix v Double -> (v Double, M.Matrix v Double)
weightedMeanCovMatrix ws xs | G.length ws /= M.rows xs = error $ (show $ M.rows xs) ++ "/" ++ (show $ G.length ws)
                            | otherwise = (means, covMat)
  where
    means = G.generate n $ \i -> meanWeighted $ G.zip (G.convert $ xs `M.takeColumn` i) ws
    covMat = MM.create $ do
        mat <- MM.new (n,n)
        forM_ [0..n-1] $ \i ->
            forM_ [i..n-1] $ \j -> do
                let cov = covWeighted ws (xs `M.takeColumn` i, means G.! i) (xs `M.takeColumn` j, means G.! j) 
                MM.unsafeWrite mat (i,j) cov
                MM.unsafeWrite mat (j,i) cov
        return mat

    n = M.cols xs
