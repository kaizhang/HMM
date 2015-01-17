{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FunctionalDependencies #-}
{-# LANGUAGE BangPatterns #-}
module AI.HMM.Class
    ( HMM(..)
    , loglikFromScales
    ) where

import Control.Monad.ST (runST)
import Control.Monad (forM_, foldM)
import Data.List (maximumBy)
import Data.Ord (comparing)
import Data.STRef
import qualified Data.Matrix.Unboxed as M
import qualified Data.Matrix.Unboxed.Mutable as MM
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Generic.Mutable as GM
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Unboxed.Mutable as UM

import AI.Function (logSumExpM)

-- discrete state and discrete time hidden markov model
class HMM hmm observ | hmm -> observ where
    -- | number of states
    nSt :: hmm -> Int

    -- | initial state distribution
    π :: hmm -> Int -> Double

    -- | transition probability
    a :: hmm -> Int -> Int -> Double
    
    -- | emission probability and its log form
    b, b' :: hmm -> Int -> observ -> Double
    b' hmm s = log . b hmm s

    baumWelch :: G.Vector v observ => v observ -> hmm -> hmm

    viterbi :: G.Vector v observ => hmm -> v observ -> ([Int], Double)
    viterbi h ob = foldl track ([],undefined) . init . G.foldl' f [] $ ob
      where
        f [] o = [U.generate n $ \i -> (log (π h i) + b' h i o, -1)]
        f acc@(l_t:_) o =
            let score i j = fst (G.unsafeIndex l_t i) + log (a h i j)
                vec = U.generate n $ \j -> let (i, x) = maximumBy (comparing snd) . map (\i' -> (i', score i' j)) $ [0..n-1]
                                           in (x + b' h j o, i)
            in vec : acc
        track ([], _) v = let ((p, prev), i) = G.maximumBy (comparing (fst.fst)) $ G.zip v $ G.enumFromN 0 n
                          in ([prev, i], p)
        track (path@(i:_), p) v = (snd (G.unsafeIndex v i) : path, p)
        n = nSt h

    forward :: G.Vector v observ => hmm -> v observ -> (M.Matrix Double, U.Vector Double)
    forward h ob = runST $ do
        mat <- MM.new (r,c)
        scales <- GM.replicate c 0
        
        -- update first column
        forM_ [0..r-1] $ \i -> do 
            let x = π h i * b h i (G.head ob)
            MM.unsafeWrite mat (i,0) x
            GM.unsafeRead scales 0 >>= GM.unsafeWrite scales 0 . (+x)
        GM.unsafeRead scales 0 >>= GM.unsafeWrite scales 0 . (1/)
        -- normalize
        sc0 <- GM.unsafeRead scales 0
        forM_ [0..r-1] $ \i -> MM.unsafeRead mat (i,0) >>= MM.unsafeWrite mat (i,0) . (*sc0)

        -- update the rest of columns
        forM_ [1..c-1] $ \t -> do
            forM_ [0..r-1] $ \j -> do
                temp <- newSTRef 0
                forM_ [0..r-1] $ \i -> do
                    α_it' <- MM.unsafeRead mat (i,t-1)
                    modifySTRef' temp (+ α_it' * a h i j)
                modifySTRef' temp (* b h j (ob `G.unsafeIndex` t))
                x <- readSTRef temp
                MM.unsafeWrite mat (j,t) x
                GM.unsafeRead scales t >>= GM.unsafeWrite scales t . (+x)
            GM.unsafeRead scales t >>= GM.unsafeWrite scales t . (1/)
            -- normalize
            sc <- GM.unsafeRead scales t
            forM_ [0..r-1] $ \j -> MM.unsafeRead mat (j,t) >>= MM.unsafeWrite mat (j,t) . (*sc)

        mat' <- MM.unsafeFreeze mat
        scales' <- G.unsafeFreeze scales
        return (mat', scales')
      where
        r = nSt h
        c = G.length ob
    {-# INLINE forward #-}

    -- | forward algorithm in log scale
    -- because in high dimension settings the pdf can be super small, the probability
    -- must be in log scale to prevent underflow. Moreover, even the probability is
    -- in log scale, it still can underflow in a long run (a long sequence of obervations).
    -- Therefore, we need to properly scale the probability at each iteration.
    forward' :: G.Vector v observ => hmm -> v observ -> (M.Matrix Double, U.Vector Double)
    forward' h ob = runST $ do
        mat <- MM.new (r,c)
        scales <- GM.new c

        -- update first column
        temp0 <- UM.new r
        forM_ [0..r-1] $ \i -> do
            let x = log (π h i) + b' h i (G.head ob)
            GM.unsafeWrite temp0 i x
        s0 <- fmap negate . logSumExpM $ temp0
        GM.unsafeWrite scales 0 s0
        -- normalize
        forM_ [0..r-1] $ \i -> GM.unsafeRead temp0 i >>= MM.unsafeWrite mat (i,0) . (+s0)

        -- update the rest of columns
        forM_ [1..c-1] $ \t -> do
            temp <- UM.new r
            forM_ [0..r-1] $ \j -> do
                sum_α_a <- foldM ( \acc i -> do
                    α_it' <- MM.unsafeRead mat (i,t-1)
                    return $! acc + exp (α_it' + log (a h i j)) ) 0 [0..r-1]
                GM.unsafeWrite temp j $ log sum_α_a + b' h j (ob `G.unsafeIndex` t)

            s <- fmap negate . logSumExpM $ temp
            GM.unsafeWrite scales t s
            -- normalize
            forM_ [0..r-1] $ \i -> GM.unsafeRead temp i >>= MM.unsafeWrite mat (i,t) . (+s)
        
        mat' <- MM.unsafeFreeze mat
        scales' <- G.unsafeFreeze scales
        return (mat', scales')
      where
        r = nSt h
        c = G.length ob
    {-# INLINE forward' #-}

    backward :: G.Vector v observ => hmm -> v observ -> U.Vector Double -> M.Matrix Double
    backward h ob scales = MM.create $ do
        mat <- MM.new (r,c)
        -- fill in last column
        forM_ [0..r-1] $ \i -> MM.unsafeWrite mat (i,c-1) $ G.last scales
        
        forM_ [c-2,c-3..0] $ \t -> do
            let sc = scales `G.unsafeIndex` t
            forM_ [0..r-1] $ \i -> do
                let f !acc j = do
                        let b_jo = b h j $ ob `G.unsafeIndex` (t+1)
                            a_ij = a h i j
                        β_jt' <- MM.unsafeRead mat (j,t+1)
                        return $ acc + b_jo * a_ij * β_jt'
                x <- foldM f 0 [0..r-1]
                MM.unsafeWrite mat (i,t) $ sc * x
        return mat
      where
        r = nSt h
        c = G.length ob
    {-# INLINE backward #-}

    -- | backward in log scale
    backward' :: G.Vector v observ => hmm -> v observ -> U.Vector Double -> M.Matrix Double
    backward' h ob scales = MM.create $ do
        mat <- MM.new (r,c)
        -- fill in last column
        forM_ [0..r-1] $ \i -> MM.unsafeWrite mat (i,c-1) $ G.last scales
        
        forM_ [c-2,c-3..0] $ \t -> do
            let sc = scales `G.unsafeIndex` t
            forM_ [0..r-1] $ \i -> do
                temp <- UM.new r
                forM_ [0..r-1] $ \j -> do
                    let b_jo = b' h j $ ob `G.unsafeIndex` (t+1)
                        a_ij = log $ a h i j
                    β_jt' <- MM.unsafeRead mat (j,t+1)
                    UM.unsafeWrite temp j $! b_jo + a_ij + β_jt'
                x <- logSumExpM temp
                MM.unsafeWrite mat (i,t) $! x + sc
        return mat
      where
        r = nSt h
        c = G.length ob
    {-# INLINE backward' #-}

    -- | log likelihood of obervations
    loglik :: G.Vector v observ => hmm -> v observ -> Double
    loglik h = loglikFromScales . snd . forward h

    {-# MINIMAL nSt, π, a, b, baumWelch #-}

-- | compute log likelihood from scales of forward probabilities vector
loglikFromScales :: U.Vector Double -> Double
loglikFromScales = G.sum . G.map (negate . log)
{-# INLINE loglikFromScales #-}
