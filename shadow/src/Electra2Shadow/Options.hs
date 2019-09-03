{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TupleSections #-}
{-# LANGUAGE NoImplicitPrelude #-}

module Electra2Shadow.Options
  ( Options
  , options
  , csvPath
  , hideControls
  , initialWindowSize
  , initialWindowPosition
  , fps
  , timeResolution
  , traceDuration
  , controlKeys
  ) where

import Protolude hiding (option)

import qualified Data.Text as Text
import Options.Applicative

data Options = Options
  { csvPath :: Maybe FilePath
  , hideControls :: Bool
  , initialWindowSize :: (Int,Int)
  , initialWindowPosition :: (Int,Int)
  , fps :: Int
  , timeResolution :: Double
  , traceDuration :: Double
  , controlKeys :: [Char]
  }

options :: IO Options
options = execParser optionsInfo

optionsInfo :: ParserInfo Options
optionsInfo = info (optionsParser <**> helper)
  fullDesc

optionsParser :: Parser Options
optionsParser = Options
  <$> (optional $ strOption
      ( long "csv"
     <> metavar "PATH"
     <> help "Csv file path."))
  <*> switch
      ( long "hide-controls")
  <*> option auto
      ( long "window-size"
     <> metavar "(INT, INT)"
     <> value (400, 600)
     <> help "Initial window dimension (width, height)."
     <> showDefault)
  <*> option auto
      ( long "window-pos"
     <> metavar "(INT, INT)"
     <> value (0, 0)
     <> help "Initial window position (x, y)."
     <> showDefault)
  <*> option auto
      ( long "fps"
     <> metavar "INT"
     <> help "Frames per second."
     <> value 60
     <> showDefault)
  <*> option auto
      ( long "time-resolution"
     <> metavar "DOUBLE"
     <> help "Frames per second."
     <> value 0.005
     <> showDefault)
  <*> option auto
      ( long "trace-duration"
     <> metavar "DOUBLE"
     <> help "Frames per second."
     <> value (1 / 25)
     <> showDefault)
  <*> option str
      ( long "control-keys"
     <> metavar "STRING"
     <> help "Keys to use for control."
     <> value ['u', 'i', 'e' ,'\195', 'p', 'o']
     <> showDefault)

