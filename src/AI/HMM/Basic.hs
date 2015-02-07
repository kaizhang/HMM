{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE MultiParamTypeClasses #-}
module AI.HMM.Basic
    ( BasicHMM
    ) where

import Control.Monad (forM_, liftM, replicateM)
import Data.STRef
import qualified Data.Vector.Unboxed as U
import qualified Data.Matrix.Unboxed as MU
import qualified Data.Matrix.Unboxed.Mutable as MUM
import System.Random.MWC

import AI.HMM.Class

data BasicHMM = BasicHMM
    { _startProb :: !(U.Vector Double)
    , _transitionProb :: !(MU.Matrix Double)
    , _emission :: !(MU.Matrix Double)
    , _boundOb :: !Int
    } deriving (Show)

instance HMMLike BasicHMM (U.Vector Int) where
    nSt = U.length . _startProb

    len _ = U.length

    isLogProb = const False

    π h i = _startProb h `U.unsafeIndex` i

    a h i j = _transitionProb h `MU.unsafeIndex` (i,j)

    b h ob t i | j < 0 || j >= _boundOb h = error "undefined observation"
               | otherwise = _emission h `MU.unsafeIndex` (i,j)
      where
        j = ob `U.unsafeIndex` t

    setInitProb x h = h {_startProb=x}

    setTransProb x h = h {_transitionProb=MU.fromVector (nSt h,nSt h) x}

    setEmProb ob ps h = h {_emission=em'}
      where
        em' = MUM.create $ do
            mat <- MUM.replicate (r, c) 0
            forM_ [0 .. l-1] $ \t -> do
                let o = ob `U.unsafeIndex` t
                forM_ [0 .. r-1] $ \i -> do
                    let γ_it = ps `U.unsafeIndex` (i*l+t)
                    MUM.unsafeRead mat (i,o) >>= MUM.unsafeWrite mat (i,o) . (+γ_it)
            forM_ [0..r-1] $ \i -> do
                temp <- newSTRef 0
                forM_ [0..c-1] $ \j -> MUM.unsafeRead mat (i,j) >>= modifySTRef' temp . (+) 
                s <- readSTRef temp
                forM_ [0..c-1] $ \j -> MUM.unsafeRead mat (i,j) >>= MUM.unsafeWrite mat (i,j) . (/s)
            return mat

        r = nSt h
        c = _boundOb h
        l = len h ob

    randHMM g states obs = do
        let randProbVector = liftM normalize . uniformVector g
            normalize xs = U.map (/ U.sum xs) xs
        startProb <- randProbVector states
        transProb <- liftM MU.fromRows $ replicateM states $ randProbVector states
        em <- liftM MU.fromRows $ replicateM states $ randProbVector obs
        return $ BasicHMM startProb transProb em obs

exampleData :: U.Vector Int
exampleData = U.fromList [1,2,3,0,1,1,1,0,2,3,2,2,2,3,0,0,0,1,1,1,2,2,2,0,1,1]

{-
test :: IO ()
test = do h <- fitHMM exampleData 2 4 200 :: IO BasicHMM
          print h
-}
