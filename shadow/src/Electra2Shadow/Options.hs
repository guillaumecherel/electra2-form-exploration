{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TupleSections #-}
{-# LANGUAGE NoImplicitPrelude #-}

module Electra2Shadow.Options
  ( Options
  , options
  , display
  , csvPath
  , configFile
  , hideControls
  , fps
  , timeResolution
  , traceDuration
  , controlKeys
  ) where

import Protolude hiding (option)

import           Numeric.Natural (Natural)
import           Options.Applicative

data Options = Options
  { configFile :: FilePath
  , display :: Maybe (Natural, Natural, Natural, Natural)
  , csvPath :: Maybe FilePath
  , hideControls :: Maybe Bool
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
  <*> (optional $ option auto
      ( long "display"
     <> metavar "(INT, INT, INT, INT)"
     <> help "Initial window dimension (width, height, x, y). If not specified, display is fullscreen."
     <> showDefault))
  <*> (optional $ strOption
      ( long "csv"
     <> metavar "PATH"
     <> help "Csv file path."))
  <*> (optional $ option auto
      ( long "hide-controls"))
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

