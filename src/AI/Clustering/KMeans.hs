{-# LANGUAGE FlexibleContexts  #-}

module AI.Clustering.KMeans
    ( kmeans
    , kmeansWith
    , forgyMethod
    ) where

import Control.Monad (forM_)
import Control.Monad.Primitive
import qualified Data.Matrix.Generic as M
import qualified Data.Matrix.Generic.Mutable as MM
import Data.Ord (comparing)
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Unboxed.Mutable as UM
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Generic.Mutable as GM
import Data.List (minimumBy)
import System.Random.MWC (uniformR, Gen)

-- | Lloyd's algorithm, also known as K-means algorithm
kmeans :: (G.Vector v1 Double, G.Vector v2 Int, Eq (v2 Int), PrimMonad m)
       => Gen (PrimState m)
       -> Int                 -- ^ number of clusters
       -> M.Matrix v1 Double  -- ^ each row represents a point
       -> m (v2 Int)          -- ^ membership vector
kmeans g k mat = do
    initial <- forgyMethod g k mat
    return $ kmeansWith mat initial

-- | Lloyd's algorithm, also known as K-means algorithm
kmeansWith :: (G.Vector v1 Double, G.Vector v2 Int, Eq (v2 Int))
           => M.Matrix v1 Double  -- ^ each row represents a point
           -> M.Matrix v1 Double  -- ^ initial set of k means
           -> v2 Int              -- ^ membership vector
kmeansWith mat initial | d /= M.cols initial || k > n = error "check input"
                   | otherwise = iter initial G.empty
  where
    iter means member = let member' = membership means
                        in if member' == member
                              then member
                              else iter (updateMeans member') member'
    membership means = G.generate n $ \i -> let x = M.takeRow mat i
                                            in fst $ minimumBy (comparing snd) $ zip [0..k-1] $ map (euclidean x) $ M.toRows means
    updateMeans member = MM.create $ do
        m <- MM.replicate (k,d) 0.0
        count <- UM.replicate k 0.0
        forM_ [0..n-1] $ \i -> do
            let x = member G.! i
            GM.unsafeRead count x >>= GM.unsafeWrite count x . (+1)
            forM_ [0..d-1] $ \j ->
                MM.unsafeRead m (x,j) >>= MM.unsafeWrite m (x,j) . (+ mat M.! (i,j))
        forM_ [0..k-1] $ \i -> do
            c <- GM.unsafeRead count i
            forM_ [0..d-1] $ \j ->
                MM.unsafeRead m (i,j) >>= MM.unsafeWrite m (i,j) . (/c)
        return m

    n = M.rows mat
    k = M.rows initial
    d = M.cols mat
{-# INLINE kmeansWith #-}

-- * Initialization methods

-- | The Forgy method randomly chooses k observations from the data set and uses
-- these as the initial means
forgyMethod :: (PrimMonad m, G.Vector v a) => Gen (PrimState m) -> Int -> M.Matrix v a -> m (M.Matrix v a)
forgyMethod g k mat | k > n = error "k is larger than sample size"
                    | otherwise = do
                        vec <- sample g n . U.enumFromN 0 $ n
                        return . M.fromRows . map (M.takeRow mat) . G.toList $ vec
  where
    n = M.rows mat

sample :: (PrimMonad m, G.Vector v a) => Gen (PrimState m) -> Int -> v a -> m (v a)
sample g n xs = do
    v <- G.thaw xs
    forM_ [0..n-1] $ \i -> do
        j <- uniformR (i, lst) g
        GM.unsafeSwap v i j
    G.unsafeFreeze . GM.take n $ v
  where
    lst = G.length xs - 1
{-# INLINE sample #-}
    
euclidean :: G.Vector v Double => v Double -> v Double -> Double
euclidean xs = sqrt . G.sum . G.zipWith (\x y -> (x - y)**2) xs
{-# INLINE euclidean #-}
