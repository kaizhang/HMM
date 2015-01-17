import AI.HMM.Basic
import AI.HMM.Class
import qualified Data.HashMap.Strict as HM
import qualified Data.Matrix.Unboxed as M
import qualified Data.Vector.Unboxed as U
import Data.List (group)

import Debug.Trace

readHMM = do
    f1 <- readFile "data/transitionMatrix.txt"
    f2 <- readFile "data/emissionMatrix.txt"
    f3 <- readFile "data/initialStateDistribution.txt"

    let a = M.fromLists . map (map read . words) . lines $ f1
        b = M.fromLists . map (map read . words) . lines $ f2
        π = U.fromList . map read . lines $ f3

    return $ BasicHMM π a b

main = do
    fl <- readFile "data/observations.txt"
    let ob = U.fromList $ map read . words $ fl
    hmm <- readHMM
    loop ob hmm
    print $ baumWelch ob hmm
--    putStrLn $ map ((\x -> HM.lookupDefault undefined x states) . head) $ group $ viterbi ob hmm
  where
    loop o h = let (_, sc) = forward h o
                   h' = baumWelch o h
               in traceShow (loglikFromScales sc) $ loop o h'
