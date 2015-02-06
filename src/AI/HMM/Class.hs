{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FunctionalDependencies #-}

module AI.HMM.Class where

import Control.Monad (forM_, foldM)
import Control.Monad.ST (runST)
import Data.STRef
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Unboxed.Mutable as UM

class Length a where
    len :: a -> Int

class Length ob => HMM h ob | h -> ob where


    -- | number of states
    nSt :: h -> Int

    -- | get and set initial probability vector
    getInitProb :: h -> U.Vector Double
    setInitProb :: U.Vector Double -> h -> h

    -- | get and set transition probability matrix
    getTransProb :: h -> U.Vector Double
    setTransProb :: U.Vector Double -> h -> h

    setEmProb :: ob -> U.Vector Double -> h -> h

    -- | initial probability at given state
    π :: h -> Int -> Double
    π' :: h -> Int -> Double

    -- | transition probability from state i to j, a i j = P(s_t = j | s_t-1 = i)
    a :: h -> Int -> Int -> Double
    a' :: h -> Int -> Int -> Double

    -- | emission probability, b o t i = P(O_t | s_t = i)
    b :: h -> ob -> Int -> Int -> Double
    b' :: h -> ob -> Int -> Int -> Double

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
        c = len ob
    {-# INLINE forward #-}

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
        c = len ob
    {-# INLINE backward #-}

    baumWelch :: ob -> h -> h
    baumWelch ob h = setInitProb ini' . setTransProb trans' . setEmProb ob γ $ h
      where
        ini' = U.slice 0 r γ
        trans' = U.create $ do
            mat <- UM.replicate (r*r) 0
            forM_ [1 .. c-1] $ \t ->
                forM_ [0..r-1] $ \i -> do
                    let α_it' = fw `U.unsafeIndex` (i*c+t-1)
                    forM_ [0..r-1] $ \j -> do
                        let a_ij = a h i j
                            b_jo = b h ob t j
                            β_jt = bw `U.unsafeIndex` (j*c+t)
                        UM.unsafeRead mat (i*c+j) >>= UM.unsafeWrite mat (i*c+j) .
                            (+) (α_it' * a_ij * b_jo * β_jt)
            normalize mat
            return mat

        (fw, scales) = forward h ob
        bw = backward h ob scales
        γ = U.generate (r*c) $ \i -> fw `U.unsafeIndex` i * bw `U.unsafeIndex` i / scales `U.unsafeIndex` (i `mod` r)
        r = nSt h
        c = len ob 

        normalize mat = forM_ [0..r-1] $ \i -> do
            temp <- newSTRef 0
            forM_ [0..c-1] $ \j -> UM.unsafeRead mat (i*c+j) >>= modifySTRef' temp . (+) 
            s <- readSTRef temp
            forM_ [0..c-1] $ \j -> UM.unsafeRead mat (i*c+j) >>= UM.unsafeWrite mat (i*c+j) . (/s)
