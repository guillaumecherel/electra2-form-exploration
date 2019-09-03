{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TupleSections #-}
{-# LANGUAGE NoImplicitPrelude #-}

module Main where

import Protolude hiding (option)

import qualified Data.ByteString.Lazy as BL
import           Graphics.Gloss (play)
import           Electra2Shadow
import           Electra2Shadow.Options

main :: IO ()
main = do
  op <- options
  controls <- case (csvPath op) of
    Nothing -> return $controlsFromInputValues 1 (-2) 0.5 0 0 0
    Just p -> do
      csv <- BL.readFile p
      let (names, ctrlMap) = controlMapFromCSV csv
      putStrLn $ (show names :: Text)
      return $ controlsFromMap ctrlMap 0 0 0 0 0
  play
    (displayMode op)
    backgroundColor
    (fps op)
    (initialWorld op controls)
    -- (initialWorld $ )
    view
    (updateInputs op)
    (updateTime op)

csvFilePath :: FilePath
csvFilePath = "data/resultsDirectSampling3.csv"
