{-# LANGUAGE DeriveGeneric #-}
module AI.GaussianHMM.Types where

import Data.Binary

import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Storable as S

import qualified Data.Matrix.Unboxed as MU
import qualified Data.Matrix.Storable as MS
import qualified Data.Matrix.Generic as M

import GHC.Generics (Generic)
import AI.MVN

type Observation = MS.Matrix Double   -- ^ n x p matrix

-- log form
data GaussianHMM = GaussianHMM
    { _startProb :: !(U.Vector Double)
    , _transitonProb :: !(MU.Matrix Double)
    , _emission :: !(V.Vector MVN)
    } deriving (Show, Generic)

instance Binary GaussianHMM

π :: GaussianHMM -> Int -> Double
π (GaussianHMM start _ _)  = U.unsafeIndex start
{-# INLINE π #-}

a :: GaussianHMM -> (Int, Int) -> Double
a (GaussianHMM _ trans _) = M.unsafeIndex trans
{-# INLINE a #-}

b :: GaussianHMM -> Int -> S.Vector Double -> Double
b (GaussianHMM _ _ e) i = logPDF (V.unsafeIndex e i)
{-# INLINE b #-}

nSt :: GaussianHMM -> Int
nSt (GaussianHMM s _ _) = U.length s
{-# INLINE nSt #-}
