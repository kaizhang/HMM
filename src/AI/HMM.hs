{-# LANGUAGE GADTs #-}
{-# LANGUAGE BangPatterns #-}

module AI.HMM where

import Control.Arrow ((***))
import Control.Monad (forM_, foldM)
import Control.Monad.ST (ST, runST)
import Control.Monad.State.Lazy
import Data.Hashable (Hashable)
import qualified Data.HashMap.Strict as HM
import qualified Data.Matrix.Unboxed as M
import qualified Data.Matrix.Unboxed.Mutable as MM
import Data.STRef
import qualified Data.Vector as V
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Unboxed as U
import Data.List (foldl', maximumBy, iterate)
import Data.Ord (comparing)
import System.Random.MWC

data HMM s o where
    HMM :: (Hashable s, Hashable o)
        => { _stateToId :: !(HM.HashMap s Int)
           , _idToState :: !(HM.HashMap Int s)
           , _eventToId :: !(HM.HashMap o Int)
           , _idToEvent :: !(HM.HashMap Int o)
           , _hmm :: !HMM'
           }
        -> HMM s o

data HMM' = HMM'
    { _startProb :: !(U.Vector Double)
    , _transitonProb :: !(M.Matrix Double)
    , _emission :: !(M.Matrix Double)
    } deriving (Show)

-- example 1: biased coin
biasedCoin = let π = U.fromList [0.85, 0.15]
                 a = M.fromLists [ [0.3, 0.7]
                                 , [0.1, 0.9]
                                 ] 
                 b = M.fromLists [ [0.4, 0.6]
                                 , [0.5, 0.5]
                                 ]
             in HMM' π a b

test = do
    let s = toSeed $ U.fromList [2,3,45,2,3,2,3,3,23,4,5,23,2]
        initial = HMM' (U.fromList [0.1,0.9])
                       (M.fromLists [[0.2,0.8],[0.5,0.5]])
                       (M.fromLists [[0.1,0.9],[0.1,0.9]])
        (path, ob) = U.unzip $ U.fromList $ generate s biasedCoin 20
    print $ iterate (baumWelchIter ob) initial !! 10

generate :: Seed -> HMM' -> Int -> [(Int, Int)]
generate seed (HMM' π a b) n = runST $ do
    g <- restore seed
    initialState <- pick g π
    let st = do
            (s,t) <- get
            if t >= n
               then return []
               else do
                   s' <- lift $ pick g $ a `M.takeRow` s
                   put (s', t+1)
                   o <- lift $ pick g $ b `M.takeRow` s
                   os <- st
                   return $ (s, o) :os
    evalStateT st (initialState, 0)
  where
    pick g' vec = do p <- uniform g'
                     go 0 0 p
      where
        go acc i p' | acc + vec G.! i >= p' = return i
                    | otherwise = go (acc + vec G.! i) (i+1) p'

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

-- | compute log likelihood from scales of forward probabilities vector
loglikFromScales :: V.Vector Double -> Double
loglikFromScales = V.sum . V.map (negate . log)

backward :: U.Vector Int -> HMM' -> U.Vector Double -> M.Matrix Double
backward ob (HMM' π a b) scales = MM.create $ do
    mat <- MM.new (r,c)
    -- fill in last column
    forM_ [0..r-1] $ \i -> MM.unsafeWrite mat (i,c-1) $ G.last scales
    
    forM_ [c-2,c-3..0] $ \t -> do
        let sc = scales `G.unsafeIndex` t
        forM_ [0..r-1] $ \i -> do
            let f !acc j = do
                    let b_jo = b `M.unsafeIndex` (j, ob `G.unsafeIndex` (t+1))
                        a_ij = a `M.unsafeIndex` (i,j)
                    β_jt' <- MM.unsafeRead mat (j,t+1)
                    return $ acc + b_jo * a_ij * β_jt'
            x <- foldM f 0 [0..r-1]
            MM.unsafeWrite mat (i,t) $ sc * x
    return mat
  where
    r = G.length π
    c = G.length ob
    
baumWelchIter :: U.Vector Int -> HMM' -> HMM'
baumWelchIter ob h@(HMM' π a b) = HMM' π' a' b'
  where
    π' = G.map (/(scales `G.unsafeIndex` 0)) $ G.zipWith (*) (fw `M.takeColumn` 0) (bw `M.takeColumn` 0)
    a' = MM.create $ do
        mat <- MM.replicate (n,n) 0
        forM_ [1 .. G.length ob - 1] $ \t -> do
            let o = ob `G.unsafeIndex` t
            forM_ [0..n-1] $ \i -> do
                let α_it' = fw `M.unsafeIndex` (i,t-1)
                forM_ [0..n-1] $ \j -> do
                    let a_ij = a `M.unsafeIndex` (i,j)
                        b_jo = b `M.unsafeIndex` (j,o)
                        β_jt = bw `M.unsafeIndex` (j,t)
                    MM.unsafeRead mat (i,j) >>= MM.unsafeWrite mat (i,j) .
                        (+) (α_it' * a_ij * b_jo * β_jt)
        normalizeByRow n n mat
        return mat

    b' = MM.create $ do
        mat <- MM.replicate (n, M.cols b) 0
        forM_ [0 .. G.length ob - 1] $ \t -> do
            let o = ob `G.unsafeIndex` t
                sc = scales `G.unsafeIndex` t
            forM_ [0 .. n-1] $ \i -> do
                let α_it = fw `M.unsafeIndex` (i,t)
                    β_it = bw `M.unsafeIndex` (i,t)
                MM.unsafeRead mat (i,o) >>= MM.unsafeWrite mat (i,o) .
                    (+) (α_it * β_it / sc)
        normalizeByRow n (M.cols b) mat
        return mat

    (fw, scales) = forward ob h
    bw = backward ob h scales
    n = G.length π

normalizeByRow :: Int -> Int -> MM.MMatrix s Double -> ST s ()
normalizeByRow r c mat = forM_ [0..r-1] $ \i -> do
    temp <- newSTRef 0
    forM_ [0..c-1] $ \j -> MM.unsafeRead mat (i,j) >>= modifySTRef' temp . (+) 
    s <- readSTRef temp
    forM_ [0..c-1] $ \j -> MM.unsafeRead mat (i,j) >>= MM.unsafeWrite mat (i,j) . (/s)

viterbi :: U.Vector Int -> HMM' -> [Int]
viterbi ob (HMM' π a b) = tail . foldl track [] . G.foldl f [] $ ob
  where
    f [] o = [V.generate n $ \i -> (log (π `G.unsafeIndex` i) + log (b `M.unsafeIndex` (i,o)), -1)]
    f acc@(l_t:_) o =
        let score i j = fst (G.unsafeIndex l_t i) + log (M.unsafeIndex a (i,j))
            vec = V.generate n $ \j -> let (i, x) = maximumBy (comparing snd) . map (\i' -> (i', score i' j)) $ [0..n-1]
                                       in (x + log (M.unsafeIndex b (j, o)), i)
        in vec : acc
    track [] v = let ((_, prev), i) = G.maximumBy (comparing (fst.fst)) $ G.zip v $ G.enumFromN 0 n
                 in [prev, i]
    track path@(i:_) v = snd (G.unsafeIndex v i) : path
    n = U.length π


{-
forward :: [Int] -> HMM s o -> M.Matrix Double
forward ob (HMM _ _ π a b) = MM.create $ do
    let r = U.length π
        c = length ob
    mat <- MM.new (r,c)
    scales <- UM.new c
    temp <- newSTRef 0  -- used to store sum
    
    -- update first column
    forM_ [0..r-1] $ \i -> do 
        let x = U.unsafeIndex π i * M.unsafeIndex b (i, head ob)
        MM.unsafeWrite mat (i,0) x
        modifySTRef' temp (+x)
    sc₀ <- (1/) <$> readSTRef temp
    UM.unsafeWrite scales 0 sc₀
    -- normalize
    forM_ [0..r-1] $ \i -> do
        x <- MM.unsafeRead mat (i,0)
        MM.unsafeWrite mat (i,0) $ x * sc₀
    
    return mat
    -}
