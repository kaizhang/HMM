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

hmmExample :: (GaussianHMM, Observation)
hmmExample = (hmm, obs)
  where
    hmm = GaussianHMM (U.fromList $ map log [0.5,0.5])
                      (M.fromLists $ ((map.map) log) [[0.1,0.9],[0.5,0.5]])
                      (V.fromList [m1,m2])
    m1 = mvn (vector [1]) (matrix 1 [1])
    m2 = mvn (vector [-1]) (matrix 1 [1])
    obs = M.fromLists $ map return [ -1.6835, 0.0635, -2.1688, 0.3043, -0.3188
                                   , -0.7835, 1.0398, -1.3558, 1.0882, 0.4050 ]


test = do
    let (hmm, obs) = hmmExample
    loop obs hmm 0
  where
    loop o h i | i > 200 = 1
               | otherwise = let h' = fst $ baumWelch o h
                                 (f, sc) = forward h o
                             in traceShow (f,G.sum sc) $ loop o h' (i+1)

train :: Observation -> Int -> Int -> GaussianHMM
train ob k n = run (initHMM,sc) 0
 where
   run (!hmm,s) !i | i > n = hmm
                   | otherwise =
                         let (hmm', sc') = baumWelch ob hmm
                         in traceShow (G.sum sc') $ run (hmm',sc') (i+1)
   initHMM = runST $ do
       g <- create
       kMeansInitial g ob k
   (_, sc) = forward initHMM ob

main = do
    let x = U.fromList [-155761.70085822034,-0.26605989874656655,-137946.9747534489,-1.454116005045439,-17814.877322709308,-2.2297494485961433,-0.15121793788899573,-3.4178055548950157,-17833.409256092447,-0.26605990643805666,-18.683151321028085,-1.454116012736929,-17814.795374554655,-2.9702413653215753,-6.926978323351274e-2,-4.158297471620448,-17818.755760800184,-0.28400031479570864,-4.029656028763902,-1.472056421094581,-17817.539730861983,-0.3279215567753966,-2.813626090563075,-1.515977663074269,-17814.751485235454,-3.9524988239600782,-2.538046403471017e-2,-5.140554930258951,-17826.96557454405,-0.26606473453025936,-12.239469772629661,-1.4541208408291317,-17814.745747346227,-4.2059209943346945,-1.9642574808625746e-2,-5.393977100633567,-155774.46838165057]
    v <- U.thaw x
    y <- logSumExpM v
    print y
