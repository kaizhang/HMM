import AI.HMM.Gaussian
import qualified Data.Vector.Unboxed as U
import qualified Data.Matrix.Unboxed as MU
import qualified Data.Vector as V
import AI.HMM.Type
import Data.Default.Class
import Control.Lens

example :: (GaussianHMM, MU.Matrix Double)
example = (hmm, obs)
  where
    hmm = GaussianHMM (U.fromList $ map log [0.5,0.5])
                      (MU.fromLists $ (map.map) log [[0.1,0.9],[0.5,0.5]])
                      (V.fromList [m1,m2])
    m1 = mvn (U.fromList [1]) (U.fromList [1])
    m2 = mvn (U.fromList [-1]) (U.fromList [1])
    obs = MU.fromLists $ map return [ -1.6835, 0.0635, -2.1688, 0.3043, -0.3188
                                   , -0.7835, 1.0398, -1.3558, 1.0882, 0.4050 ]

test :: IO ()
test = do
    let h = fitHMM (snd example) $ def & initialization .~ Fixed (fst example) 
                                       & mStepOpt.covEstimator .~ DiagonalCov
    print (h :: GaussianHMM)

main = test
