{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TupleSections #-}
{-# LANGUAGE NoImplicitPrelude #-}

module Main where

import Protolude hiding (option)

import qualified Data.ByteString.Lazy as BL
import qualified Data.Vector as Vector
import           Graphics.Gloss (play)
import           Electra2Shadow
import           Electra2Shadow.Options
import           Electra2Shadow.Config
import qualified Electra2Shadow.Control as Control

main :: IO ()
main = do
  op <-options
  config <- getConfig op
  controls <- case (Electra2Shadow.Config.csvPath config) of
    Nothing -> return $ Control.fromInputValues
                         (v1name config, 1)
                         (v2name config, (-2))
                         (v3name config, 0.5)
                         (phi1name config, 0)
                         (phi2name config, 0)
                         (phi3name config, 0)
    Just p -> do
      csv <- BL.readFile p
      let names = Control.getControlName identity
              <$> Vector.toList (initialControls config)
      let ctrlMap = Control.mapFromCSV
            names
            (v1name config)
            (v2name config)
            (v3name config)
            (phi1name config)
            (phi2name config)
            (phi3name config)
            csv
      putStrLn $ (show names :: Text)
      return $ Control.fromMap ctrlMap $ initialControls config
  play
    (displayMode config)
    backgroundColor
    (fromIntegral $ Electra2Shadow.Config.fps config)
    (initialWorld config controls)
    -- (initialWorld $ )
    view
    (updateInputs config)
    (updateTime config)

csvFilePath :: FilePath
csvFilePath = "data/resultsDirectSampling3.csv"
