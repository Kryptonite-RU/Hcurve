{-# LANGUAGE TypeApplications #-}

-- (c) 2020 Kryptonite

module Curves
    ( greyEncode
    , greyDecode
    , hcurveEncode
    , hcurveDecode
    , hilbertEncode
    , hilbertDecode
    , zorderEncode
    , zorderDecode
    ) where

import Control.Lens.Operators
import Data.Bits
import Data.List

-- Composing a number from non-zero bit numbers
numFromBits :: (Num a, Bits a) => [Int] -> a
numFromBits = foldl setBit zeroBits

-- Grey decode : position -> index
greyDecode :: (Integral a, Bits a) => Int -> a -> a
greyDecode 0 _ = 0
greyDecode _ 0 = 0
greyDecode n i = xor i $ greyDecode (pred n) (shiftR i 1)

-- Grey encode : index -> position
greyEncode :: (Integral a, Bits a) => Int -> a -> a
greyEncode d i = mod (xor i $ div i 2) (2^d)

-- get list of n bits (as booleans) of a number 'a'
binaryRepr :: (Num a, Bits a) => Int -> a -> [Bool]
binaryRepr n a = [0..n-1] <&> testBit a

fromBinaryRepr :: (Num a, Bits a) => [Bool] -> a
fromBinaryRepr bs = zip [0..] bs & filter snd <&> fst & numFromBits

-- update head of list
updateHead :: (a -> a) -> [] a -> [] a
updateHead f (x:xs) = f x : xs

updateListElt :: Int -> (a -> a) -> [] a -> [] a
updateListElt 0 f l      = updateHead f l
updateListElt n f (x:xs) = x : updateListElt (pred n) f xs

-- entry point of a codepth 1 subcube in HCurve
entryPoint :: Int -> Int -> [Bool] -> [Integer]
entryPoint dim depth mainBits =
  let upd = updateHead (if odd . length . filter id $ mainBits then id else xor 1)
  in  mainBits <&> (\b -> (if b then 0 else 1) * (2 ^ pred depth - 1)) & upd

-- Slice list to list of k-length lists
by :: Int -> [a] -> [[a]]
by k = unfoldr (\l -> if null l then Nothing else Just $ splitAt k l)

-- index -> position in cube
hcurveEncode :: Int -> Int -> Integer -> [Integer]
hcurveEncode 0 _ _ = []
hcurveEncode d 0 _ = replicate d 0
hcurveEncode dim depth n =
  let n'             = mod (div n subCubeSize) (2^dim)
      mainBits       = binaryRepr dim $ greyEncode dim n'
      directionSign  = if odd dim && depth == 2 then -1 else 1
      subCubeSize    = 2 ^(dim * pred depth)
      ep             = entryPoint dim depth mainBits
      pointShift     = hcurveDecode dim (pred depth) ep
      subCubePoint   = hcurveEncode dim (pred depth) (mod (directionSign * n + pointShift) subCubeSize)
  in zipWith (+) (mainBits <&> \ b -> if b then 2^pred depth else 0) subCubePoint

-- position in cube -> index
hcurveDecode :: Int -> Int -> [Integer] -> Integer
hcurveDecode 0 _ _ = 0
hcurveDecode _ 0 _ = 0
hcurveDecode dim depth xs =
  let mainBits        = xs <&> flip testBit (pred depth)
      directionSign   = if odd dim && depth == 2 then -1 else 1
      subCubeSize     = 2 ^(dim * pred depth)
      ep              = entryPoint dim depth mainBits
      pointShift      = hcurveDecode dim (pred depth) ep
      pointNum        = hcurveDecode dim (pred depth) xs - pointShift
      n'              = mainBits
         & flip zip [0..] & filter fst <&> snd
         & numFromBits @Integer
         & greyDecode dim
  in  (mod (pointNum * directionSign) subCubeSize + subCubeSize * n')
      & flip mod (2 ^ (dim * depth))


hilbertEncode :: Int -> Int -> Integer -> [Integer]
hilbertEncode dim depth n =
  let x              = n & binaryRepr (dim*depth) & by dim & transpose <&> fromBinaryRepr @Integer & (reverse $!)
      t              = shiftR (last x) 1
      x'             = zipWith xor x (t:x)
      exchHead i j l =
        let t        = (xor (head l) (l !! j)) .&. (2^i-1)
        in             updateHead (xor t) . updateListElt j (xor t) $ l
      upd' i xs j    = xs & if testBit (xs !! j) i then updateHead (xor (2^i-1)) else exchHead i j
      upd xs i       = foldl (upd' i) xs $ reverse [0..dim-1]
      x''            = foldl upd x' [1..depth]
      x'''           = last x'' : init x''
  in  x'''

hilbertDecode :: Int -> Int -> [Integer] -> Integer
hilbertDecode dim depth x''' =
  let x''            = tail x''' ++ [head x''']
      x'             = foldr upd x'' [1..depth]
      upd i xs       = foldr (upd' i) xs $ reverse [0..dim-1]
      upd' i j xs    = xs & if testBit (xs !! j) i then updateHead (xor (2^i-1)) else exchHead i j
      exchHead i j l =
        let t        = (xor (head l) (l !! j)) .&. (2^i-1)
        in             updateHead (xor t) . updateListElt j (xor t) $ l
      ge _ []        = []
      ge acc (x:xs)  = let y = xor acc x in y : ge y xs
      x              = ge 0 x'
      t              = greyDecode depth $ (flip div 2) $ last x
      xs             = xor t <$> x
  in  xs <&> binaryRepr depth & reverse
      & transpose & concat
      & fromBinaryRepr @ Integer


zorderEncode :: Int -> Int -> Integer -> [Integer]
zorderEncode dim depth n = n & binaryRepr (dim*depth) & by dim & transpose <&> fromBinaryRepr

zorderDecode :: Int -> Int -> [Integer] -> Integer
zorderDecode dim depth xs = xs <&> binaryRepr depth & transpose & concat & fromBinaryRepr
