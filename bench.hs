import AI.Clustering.KMeans
import System.Random.MWC
import qualified Data.Matrix.Unboxed as M
import Data.DataSets.Cluster
import qualified Data.Vector.Unboxed as U

main = do
    g <- createSystemRandom
    let xs = M.fromLists $ map (\x -> [xclaraV1 x, xclaraV2 x]) xclara :: M.Matrix Double
    r <- kmeans g 3 xs :: IO (U.Vector Int)
