{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TupleSections #-}
{-# LANGUAGE NoImplicitPrelude #-}

module Electra2Shadow.Options
  ( Options
  , options
  , csvPath
  , configFile
  , hideControls
  , initialWindowSize
  , initialWindowPosition
  , fps
  , timeResolution
  , traceDuration
  , controlKeys
  ) where

import Protolude hiding (option)

import           Numeric.Natural (Natural)
import qualified Data.Text as Text
import           Options.Applicative

data Options = Options
  { configFile :: FilePath
  , csvPath :: Maybe FilePath
  , hideControls :: Maybe Bool
  , initialWindowSize :: Maybe (Natural,Natural)
  , initialWindowPosition :: Maybe (Natural,Natural)
  , fps :: Maybe Natural
  , timeResolution :: Maybe Double
  , traceDuration :: Maybe Double
  , controlKeys :: Maybe [Char]
  }

options :: IO Options
options = execParser optionsInfo

optionsInfo :: ParserInfo Options
optionsInfo = info (optionsParser <**> helper)
  fullDesc

optionsParser :: Parser Options
optionsParser = Options
  <$> (strOption
      ( long "config"
     <> metavar "PATH"
     <> value "config.example"
     <> help "Config file path."))
  <*> (optional $ strOption
      ( long "csv"
     <> metavar "PATH"
     <> help "Csv file path."))
  <*> (optional $ option auto
      ( long "hide-controls"))
  <*> (optional $ option auto
      ( long "window-size"
     <> metavar "(INT, INT)"
     <> help "Initial window dimension (width, height)."
     <> showDefault))
  <*> (optional $ option auto
      ( long "window-pos"
     <> metavar "(INT, INT)"
     <> help "Initial window position (x, y)."
     <> showDefault))
  <*> (optional $ option auto
      ( long "fps"
     <> metavar "INT"
     <> help "Frames per second."
     <> showDefault))
  <*> (optional $ option auto
      ( long "time-resolution"
     <> metavar "DOUBLE"
     <> help "Frames per second."
     <> showDefault))
  <*> (optional $ option auto
      ( long "trace-duration"
     <> metavar "DOUBLE"
     <> help "Frames per second."
     <> showDefault))
  <*> (optional $ option str
      ( long "control-keys"
     <> metavar "STRING"
     <> help "Keys to use for control."
     -- <> value ['u', 'i', 'e' ,'\195', 'p', 'o']
     <> showDefault))

