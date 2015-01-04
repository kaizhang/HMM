import AI.HMM.Basic
import qualified Data.HashMap.Strict as HM
import qualified Data.Matrix.Unboxed as M
import qualified Data.Vector.Unboxed as U
import Data.List (group)

readHMM :: IO (HMM Char Int)
readHMM = do
    f1 <- readFile "data/transitionMatrix.txt"
    f2 <- readFile "data/emissionMatrix.txt"
    f3 <- readFile "data/initialStateDistribution.txt"

    let a = M.fromLists . map (map read . words) . lines $ f1
        b = M.fromLists . map (map read . words) . lines $ f2
        π = U.fromList . map read . lines $ f3
        states = HM.fromList $ zip ['a'..'z'] [0..25]
        observs = HM.fromList [(0,0),(1,1)]

    return $ HMM states observs π a b

main = do
    fl <- readFile "data/observations.txt"
    let ob = map read . words $ fl :: [Int]
        states = HM.fromList $ zip [0..25] ['a'..'z']
    hmm <- readHMM
    putStrLn $ map ((\x -> HM.lookupDefault undefined x states) . head) $ group $ viterbi ob hmm
