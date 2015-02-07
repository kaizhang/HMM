{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FunctionalDependencies #-}

module AI.HMM.Class
    ( HMMLike(..)
    , HMMOpt(..)
    , Initialization(..)
    ) where

import Control.Monad (forM_, foldM)
import Control.Monad.Primitive (PrimMonad, PrimState)
import Control.Monad.ST (runST)
import Data.STRef
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Unboxed.Mutable as UM
import System.Random.MWC
import Data.Default.Class

import AI.Function

import Debug.Trace

class HMMLike h ob | h -> ob where
    -- | number of states
    nSt :: h -> Int

    len :: h -> ob -> Int

    -- | probability is in log form
    isLogProb :: h -> Bool

    -- | get and set initial probability vector
    setInitProb :: U.Vector Double -> h -> h

    -- | get and set transition probability matrix
    setTransProb :: U.Vector Double -> h -> h

    setEmProb :: ob -> U.Vector Double -> h -> h

    -- | initial probability at given state
    π :: h -> Int -> Double
    π h = exp . π' h
    π' :: h -> Int -> Double
    π' h = log . π h

    -- | transition probability from state i to j, a i j = P(s_t = j | s_t-1 = i)
    a :: h -> Int -> Int -> Double
    a h i = exp . a' h i
    a' :: h -> Int -> Int -> Double
    a' h i = log . a h i

    -- | emission probability, b o t i = P(O_t | s_t = i)
    b :: h -> ob -> Int -> Int -> Double
    b h o i = exp . b' h o i
    b' :: h -> ob -> Int -> Int -> Double
    b' h o i = log . b h o i

    forward :: h -> ob -> (U.Vector Double, U.Vector Double)
    forward h ob = runST $ do
        mat <- UM.new $ r * c
        scales <- UM.replicate c 0
        
        -- update first column
        forM_ [0..r-1] $ \i -> do 
            let x = π h i * b h ob 0 i
            UM.unsafeWrite mat (i*c) x
            UM.unsafeRead scales 0 >>= UM.unsafeWrite scales 0 . (+x)
        UM.unsafeRead scales 0 >>= UM.unsafeWrite scales 0 . (1/)
        -- normalize
        sc0 <- UM.unsafeRead scales 0
        forM_ [0..r-1] $ \i -> UM.unsafeRead mat (i*c) >>= UM.unsafeWrite mat (i*c) . (*sc0)

        -- update the rest of columns
        forM_ [1..c-1] $ \t -> do
            forM_ [0..r-1] $ \j -> do
                temp <- newSTRef 0
                forM_ [0..r-1] $ \i -> do
                    α_it' <- UM.unsafeRead mat $ i*c+(t-1)
                    modifySTRef' temp (+ α_it' * a h i j)
                modifySTRef' temp (* b h ob t j)
                x <- readSTRef temp
                UM.unsafeWrite mat (j*c+t) x
                UM.unsafeRead scales t >>= UM.unsafeWrite scales t . (+x)
            UM.unsafeRead scales t >>= UM.unsafeWrite scales t . (1/)
            -- normalize
            sc <- UM.unsafeRead scales t
            forM_ [0..r-1] $ \j -> UM.unsafeRead mat (j*c+t) >>= UM.unsafeWrite mat (j*c+t) . (*sc)

        mat' <- U.unsafeFreeze mat
        scales' <- U.unsafeFreeze scales
        return (mat', scales')
      where
        r = nSt h
        c = len h ob
    {-# INLINE forward #-}

    -- | log form
    forward' :: h -> ob -> (U.Vector Double, U.Vector Double)
    forward' h ob = runST $ do
        mat <- UM.new (r*c)
        scales <- UM.new c

        -- update first column
        temp0 <- UM.new r
        forM_ [0..r-1] $ \i -> do
            let x = π' h i + b' h ob 0 i
            UM.unsafeWrite temp0 i x
        s0 <- fmap negate . logSumExpM $ temp0
        UM.unsafeWrite scales 0 s0
        -- normalize
        forM_ [0..r-1] $ \i -> UM.unsafeRead temp0 i >>= UM.unsafeWrite mat (i*c) . (+s0)

        -- update the rest of columns
        forM_ [1..c-1] $ \t -> do
            temp <- UM.new r
            forM_ [0..r-1] $ \j -> do
                temp' <- UM.new r
                forM_ [0..r-1] $ \i -> do
                    α_it' <- UM.unsafeRead mat $ i*c+t-1
                    UM.unsafeWrite temp' i $ α_it' + a' h i j
                    
                sum_α_a <- logSumExpM temp'

                UM.unsafeWrite temp j $ sum_α_a + b' h ob t j

            s <- fmap negate . logSumExpM $ temp
            UM.unsafeWrite scales t s
            -- normalize
            forM_ [0..r-1] $ \i -> UM.unsafeRead temp i >>= UM.unsafeWrite mat (i*c+t) . (+s)
        
        mat' <- U.unsafeFreeze mat
        scales' <- U.unsafeFreeze scales
        return (mat', scales')
      where
        r = nSt h
        c = len h ob
    {-# INLINE forward' #-}

    backward :: h -> ob -> U.Vector Double -> U.Vector Double
    backward h ob scales = U.create $ do
        mat <- UM.new (r*c)
        -- fill in last column
        forM_ [0..r-1] $ \i -> UM.unsafeWrite mat (i*c+c-1) $ U.last scales
        
        forM_ [c-2,c-3..0] $ \t -> do
            let sc = scales `U.unsafeIndex` t
            forM_ [0..r-1] $ \i -> do
                let f !acc j = do
                        let b_jo = b h ob (t+1) j
                            a_ij = a h i j
                        β_jt' <- UM.unsafeRead mat (j*c+t+1)
                        return $ acc + b_jo * a_ij * β_jt'
                x <- foldM f 0 [0..r-1]
                UM.unsafeWrite mat (i*c+t) $ sc * x
        return mat
      where
        r = nSt h
        c = len h ob
    {-# INLINE backward #-}

    backward' :: h -> ob -> U.Vector Double -> U.Vector Double
    backward' h ob scales = U.create $ do
        mat <- UM.new (r*c)
        -- fill in last column
        forM_ [0..r-1] $ \i -> UM.unsafeWrite mat (i*c+c-1) $ U.last scales
        
        forM_ [c-2,c-3..0] $ \t -> do
            let sc = scales `U.unsafeIndex` t
            forM_ [0..r-1] $ \i -> do
                temp <- UM.new r
                forM_ [0..r-1] $ \j -> do
                    let b_jo = b' h ob (t+1) j
                        a_ij = a' h i j
                    β_jt' <- UM.unsafeRead mat (j*c+t+1)
                    UM.unsafeWrite temp j $! b_jo + a_ij + β_jt'
                x <- logSumExpM temp
                UM.unsafeWrite mat (i*c+t) $! x + sc
        return mat
      where
        r = nSt h
        c = len h ob
    {-# INLINE backward' #-}

    baumWelch :: ob -> h -> (h, U.Vector Double)
    baumWelch ob h = (setInitProb ini' . setTransProb trans' . setEmProb ob γ $ h, U.map log scales)
      where
        ini' = U.generate r $ \i -> γ `U.unsafeIndex` (i*c)
        trans' = U.create $ do
            mat <- UM.replicate (r*r) 0
            forM_ [1 .. c-1] $ \t ->
                forM_ [0..r-1] $ \i -> do
                    let α_it' = fw `U.unsafeIndex` (i*c+t-1)
                    forM_ [0..r-1] $ \j -> do
                        let a_ij = a h i j
                            b_jo = b h ob t j
                            β_jt = bw `U.unsafeIndex` (j*c+t)
                        UM.unsafeRead mat (i*r+j) >>= UM.unsafeWrite mat (i*r+j) .
                            (+) (α_it' * a_ij * b_jo * β_jt)
            normalize mat
            return mat

        (fw, scales) = forward h ob
        bw = backward h ob scales
        γ = U.generate (r*c) $ \i -> fw `U.unsafeIndex` i * bw `U.unsafeIndex` i / scales `U.unsafeIndex` (i `mod` c)
        r = nSt h
        c = len h ob

        normalize mat = forM_ [0..r-1] $ \i -> do
            temp <- newSTRef 0
            forM_ [0..r-1] $ \j -> UM.unsafeRead mat (i*r+j) >>= modifySTRef' temp . (+) 
            s <- readSTRef temp
            forM_ [0..r-1] $ \j -> UM.unsafeRead mat (i*r+j) >>= UM.unsafeWrite mat (i*r+j) . (/s)
    {-# INLINE baumWelch #-}

    baumWelch' :: ob -> h -> (h, U.Vector Double)
    baumWelch' ob h = (setInitProb ini' . setTransProb trans' . setEmProb ob (U.map exp γ) $ h, scales)
      where
        ini' = U.generate r $ \i -> γ `U.unsafeIndex` (i*c)

        trans' = U.create $ do
            mat <- UM.new (r*r)
            forM_ [0..r-1] $ \i ->
                forM_ [0..r-1] $ \j -> do
                    let a_ij = a' h i j
                    temp <- UM.new $ c - 1
                    forM_ [1 .. c-1] $ \t -> do
                        let b_jo = b' h ob t j
                            α_it' = fw `U.unsafeIndex` (i*c+t-1)
                            β_jt = bw `U.unsafeIndex` (j*c+t)
                        UM.unsafeWrite temp (t-1) $ b_jo + α_it' + β_jt
                    x <- logSumExpM temp
                    UM.unsafeWrite mat (i*r+j) $ a_ij + x
            normalize mat
            return mat

        (fw, scales) = forward' h ob
        bw = backward' h ob scales
        γ = U.generate (r*c) $ \i -> fw `U.unsafeIndex` i + bw `U.unsafeIndex` i - scales `U.unsafeIndex` (i `mod` c)
        r = nSt h
        c = len h ob
        normalize mat = forM_ [0..r-1] $ \i -> do
            s <- logSumExpM $ UM.slice (i*r) r mat
            forM_ [0..r-1] $ \j -> UM.unsafeRead mat (i*r+j) >>= UM.unsafeWrite mat (i*r+j) . subtract s
    {-# INLINE baumWelch' #-}

    randHMM :: PrimMonad m => Gen (PrimState m) -> Int -> Int -> m h

    fitHMM :: ob -> Int -> Int -> HMMOpt h -> h
    fitHMM ob m n opt = runST $ do
        g <- restore $ _seed opt
        initialHMM <- case _initialization opt of
            Fixed inithmm -> return inithmm
            Random -> randHMM g m n
        let updateFn | isLogProb initialHMM = baumWelch'
                     | otherwise = baumWelch
        loop updateFn initialHMM 0
      where
        loop f !h !i | i >= _nIter opt = return h
                     | otherwise = do
                         let (h', scales) = f ob h
                         traceShow (negate . U.sum $ scales) $ return ()
                         loop f h' $ i+1

data HMMOpt h = HMMOpt
    { _seed :: !Seed
    , _initialization :: !(Initialization h)
    , _nIter :: !Int
    }

data Initialization h = Fixed h
                      | Random

instance Default (HMMOpt h) where
    def = HMMOpt
        { _seed = toSeed $ U.fromList [22]
        , _initialization = Random
        , _nIter = 50
        }
