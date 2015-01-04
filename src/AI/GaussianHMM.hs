{-# LANGUAGE GADTs #-}

module AI.GaussianHMM where

import Numeric.LinearAlgebra.HMatrix (Matrix, Vector, (<>), invlndet, (!), tr, asRow, vector, matrix)
import qualified Data.Vector as V
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Unboxed as U
import qualified Data.Matrix.Unboxed as M
import Data.List (maximumBy)
import Data.Ord (comparing)

-- | multivariate normal distribution
data MVN = MVN
    { _mean :: !(Vector Double)
--    , _cov :: Matrix Double
    , _invcov :: !(Matrix Double)
    , _logdet :: !Double  -- ^ log determinant of covariance matrix
    } deriving (Show)

mvn :: Vector Double -> Matrix Double -> MVN
mvn m cov = MVN m invcov logdet
  where
    (invcov, (logdet, _)) = invlndet cov

-- | log probability of MVN
logPDF :: MVN -> Vector Double -> Double
logPDF (MVN m invcov logdet) x = -0.5 * ( d * log (2*pi) + logdet
                                        + (x' <> invcov <> tr x') ! 0 ! 0
                                        )
  where
    x' = asRow $ x - m
    d = fromIntegral . G.length $ m

data GaussianHMM = GaussianHMM
    { _startProb :: !(U.Vector Double)
    , _transitonProb :: !(M.Matrix Double)
    , _emission :: !(V.Vector MVN)
    } deriving (Show)

prob :: V.Vector MVN -> Int -> Vector Double -> Double
prob vec i = logPDF (vec G.! i)

forward :: U.Vector Int -> HMM' -> (M.Matrix Double, U.Vector Double)
forward ob (HMM' π a b) = (M.fromColumns *** U.fromList) . unzip . reverse . G.foldl' f [] $ ob
  where
    f [] o = return . normalize . U.imap (\i p -> p * M.unsafeIndex b (i,o)) $ π
    f acc o = let (x, _) = head acc
                  α_t = U.generate n $ \i ->
                      U.sum (U.zipWith (*) x $ M.takeColumn a i) * M.unsafeIndex b (i,o)
              in normalize α_t : acc
    normalize xs = let sc = 1 / U.sum xs
                   in (U.map (*sc) xs, sc)
    n = U.length π

viterbi :: V.Vector (Vector Double) -> GaussianHMM -> [Int]
viterbi ob (GaussianHMM π a b) = tail . foldl track [] . G.foldl f [] $ ob
  where
    f [] o = [V.generate n $ \i -> (log (π `G.unsafeIndex` i) + prob b i o, -1)]
    f acc@(l_t:_) o =
        let score i j = fst (G.unsafeIndex l_t i) + log (M.unsafeIndex a (i,j))
            vec = V.generate n $ \j -> let (i, x) = maximumBy (comparing snd) . map (\i' -> (i', score i' j)) $ [0..n-1]
                                       in (x + prob b j o, i)
        in vec : acc
    track [] v = let ((_, prev), i) = G.maximumBy (comparing (fst.fst)) $ G.zip v $ G.enumFromN 0 n
                 in [prev, i]
    track path@(i:_) v = snd (G.unsafeIndex v i) : path
    n = U.length π

-- test:
-- Gaussian Hidden Markov Model with 2 States
--
-- Transition matrix:
--    0.1 0.9
--    0.5 0.5
--
-- Emission parameters:
-- [N(1.0, 1.0), N(-1.0, 1.0)]
--
-- Initial probabilities: [0.5000, 0.5000]
--
-- ([1, 0, 1, 0, 1, 1, 0, 1, 0, 1], -16.67738270170788)
test = do
   let hmm = GaussianHMM (U.fromList [0.5,0.5]) (M.fromLists [[0.1,0.9],[0.5,0.5]]) (V.fromList [m1,m2])
       m1 = mvn (vector [1]) (matrix 1 [1])
       m2 = mvn (vector [-1]) (matrix 1 [1])
       obs = V.fromList $ map (vector.return) [-1.6835, 0.0635, -2.1688, 0.3043, -0.3188, -0.7835, 1.0398, -1.3558, 1.0882, 0.4050]
   print $ viterbi obs hmm
