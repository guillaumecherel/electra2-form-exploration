{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TupleSections #-}
{-# LANGUAGE NoImplicitPrelude #-}

module Main where

import Protolude

import qualified Data.ByteString.Lazy as BL
import           Graphics.Gloss (play)
import           Electra2Shadow

main :: IO ()
main = do
  csv <- BL.readFile csvFilePath
  let (names, ctrlMap) = controlMapFromCSV csv
  putStrLn $ (show names :: Text)
  play
   displayMode
   backgroundColor
   fps
   (initialWorld $ controlsFromInputValues 1 (-2) 0.5 0 0 0)
   -- (initialWorld $ controlsFromMap ctrlMap 0 0 0 0 0)
   view
   updateInputs
   updateTime

csvFilePath :: FilePath
csvFilePath = "data/resultsDirectSampling3.csv"
