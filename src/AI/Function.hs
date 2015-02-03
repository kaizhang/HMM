{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE BangPatterns #-}
module AI.Function where

import Control.Monad (forM_)
import Control.Monad.Primitive
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Generic.Mutable as GM
import qualified Data.Matrix.Generic as M
import qualified Data.Matrix.Generic.Mutable as MM
import Numeric.LinearAlgebra.HMatrix (Vector, (<>), invlndet, (!), tr, asRow)
import Statistics.Sample (meanWeighted)

import Debug.Trace

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
    g (a,b) | a == 0 || b == 0 = trace "zero covariance" $ 1e-200 -- pseudocount
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
