{-# LANGUAGE BangPatterns #-}

module AI.GaussianHMM where

import Control.Monad.ST (runST)
import Numeric.LinearAlgebra.HMatrix (Matrix, Vector, (<>), invlndet, (!), tr, asRow, vector, matrix, reshape)
import qualified Data.Vector as V
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Unboxed as U
import qualified Data.Matrix.Generic as M
import System.Random.MWC

import AI.Function
import AI.GaussianHMM.Types
import AI.GaussianHMM.Algorithms

import Debug.Trace

test = do
   let hmm = GaussianHMM (U.fromList $ map log [0.5,0.5]) (M.fromLists $ ((map.map) log) [[0.1,0.9],[0.5,0.5]]) (V.fromList [m1,m2])
       m1 = mvn (vector [1]) (matrix 1 [1])
       m2 = mvn (vector [-1]) (matrix 1 [1])
       obs = M.fromLists $ map return [-1.6835, 0.0635, -2.1688, 0.3043, -0.3188, -0.7835, 1.0398, -1.3558, 1.0882, 0.4050]
       (f, sc) = forward hmm obs
       b = backward hmm obs sc
--   print $ p_S_t M.! (1,5)
   print $ train' obs 2 30
--   print $ f M.! (1,5) + b M.! (1,5) - sc G.! 5
  where

train' :: Observation -> Int -> Int -> GaussianHMM
train' ob k n = run (initHMM,sc) 0
 where
   run (!hmm,s) !i | i > n = hmm
                   | otherwise =
                         let (hmm', sc') = baumWelch ob hmm
                         in traceShow (G.sum $ G.zipWith (-) s sc') $ run (hmm',sc') (i+1)
   initHMM = runST $ do
       g <- create
       kMeansInitial g ob k
   (_, sc) = forward initHMM ob
