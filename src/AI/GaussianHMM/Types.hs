module AI.GaussianHMM.Types where

import Numeric.LinearAlgebra.HMatrix (Vector)
import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as U
import qualified Data.Matrix.Unboxed as MU
import qualified Data.Matrix.Generic as M
import qualified Data.Matrix.Storable as S

import AI.Function

type Observation = S.Matrix Double   -- ^ n x p matrix

-- log form
data GaussianHMM = GaussianHMM
    { _startProb :: !(U.Vector Double)
    , _transitonProb :: !(MU.Matrix Double)
    , _emission :: !(V.Vector MVN)
    } deriving (Show)

π :: GaussianHMM -> Int -> Double
π (GaussianHMM start _ _)  = U.unsafeIndex start
{-# INLINE π #-}

a :: GaussianHMM -> (Int, Int) -> Double
a (GaussianHMM _ trans _) = M.unsafeIndex trans
{-# INLINE a #-}

b :: GaussianHMM -> Int -> Vector Double -> Double
b (GaussianHMM _ _ e) i = logPDF (V.unsafeIndex e i)
{-# INLINE b #-}

nSt :: GaussianHMM -> Int
nSt (GaussianHMM s _ _) = U.length s
{-# INLINE nSt #-}
