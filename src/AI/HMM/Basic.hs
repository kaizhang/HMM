{-# LANGUAGE MultiParamTypeClasses #-}
module AI.HMM.Basic where

import AI.HMM.Class
import Control.Monad (forM_)
import Control.Monad.ST (ST)
import qualified Data.Matrix.Unboxed as M
import qualified Data.Matrix.Unboxed.Mutable as MM
import Data.STRef
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Unboxed as U

data BasicHMM = BasicHMM
    { _startProb :: !(U.Vector Double)
    , _transitonProb :: !(M.Matrix Double)
    , _emission :: !(M.Matrix Double)
    } deriving (Show)

instance HMM BasicHMM Int where
    nSt (BasicHMM x _ _) = U.length x
    {-# INLINE nSt #-}

    π (BasicHMM x _ _) s = x `U.unsafeIndex` s
    {-# INLINE π #-}

    a (BasicHMM _ x _) i j = x `M.unsafeIndex` (i,j)
    {-# INLINE a #-}

    b (BasicHMM _ _ x) i o = x `M.unsafeIndex` (i,o)
    {-# INLINE b #-}

    baumWelch ob h@(BasicHMM _ trans em) = BasicHMM ini' trans' em'
      where
        ini' = G.map (/(scales `G.unsafeIndex` 0)) $ G.zipWith (*) (fw `M.takeColumn` 0) (bw `M.takeColumn` 0)
        trans' = MM.create $ do
            mat <- MM.replicate (n,n) 0
            forM_ [1 .. G.length ob - 1] $ \t -> do
                let o = ob `G.unsafeIndex` t
                forM_ [0..n-1] $ \i -> do
                    let α_it' = fw `M.unsafeIndex` (i,t-1)
                    forM_ [0..n-1] $ \j -> do
                        let a_ij = trans `M.unsafeIndex` (i,j)
                            b_jo = em `M.unsafeIndex` (j,o)
                            β_jt = bw `M.unsafeIndex` (j,t)
                        MM.unsafeRead mat (i,j) >>= MM.unsafeWrite mat (i,j) .
                            (+) (α_it' * a_ij * b_jo * β_jt)
            normalizeByRow n n mat
            return mat

        em' = MM.create $ do
            mat <- MM.replicate (n, M.cols em) 0
            forM_ [0 .. G.length ob - 1] $ \t -> do
                let o = ob `G.unsafeIndex` t
                    sc = scales `G.unsafeIndex` t
                forM_ [0 .. n-1] $ \i -> do
                    let α_it = fw `M.unsafeIndex` (i,t)
                        β_it = bw `M.unsafeIndex` (i,t)
                    MM.unsafeRead mat (i,o) >>= MM.unsafeWrite mat (i,o) .
                        (+) (α_it * β_it / sc)
            normalizeByRow n (M.cols em) mat
            return mat

        (fw, scales) = forward h ob
        bw = backward h ob scales
        n = nSt h

        normalizeByRow :: Int -> Int -> MM.MMatrix s Double -> ST s ()
        normalizeByRow r c mat = forM_ [0..r-1] $ \i -> do
            temp <- newSTRef 0
            forM_ [0..c-1] $ \j -> MM.unsafeRead mat (i,j) >>= modifySTRef' temp . (+) 
            s <- readSTRef temp
            forM_ [0..c-1] $ \j -> MM.unsafeRead mat (i,j) >>= MM.unsafeWrite mat (i,j) . (/s)

