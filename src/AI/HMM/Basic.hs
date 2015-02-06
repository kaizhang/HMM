module HMM.Basic where

import qualified Data.Vector.Unboxed as U

data BasicHMM = BasicHMM
    { _startProb :: !(U.Vector Double)
    , _transitionProb :: !(U.Vector Double)
    , _emission :: !(U.Vector Double)
    }
