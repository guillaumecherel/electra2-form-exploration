{-# LANGUAGE DeriveGeneric     #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TupleSections #-}
{-# LANGUAGE NoImplicitPrelude #-}

module Electra2Shadow.Config
  where

import Protolude hiding (Type)

import           Numeric.Natural (Natural)
import           Data.Vector (Vector)
import           Dhall
import           Electra2Shadow.Control (Control)
import qualified Electra2Shadow.Options as Options
import           Electra2Shadow.Options (Options)

data Config = Config
  { csvPath :: Maybe FilePath
  , hideControls :: Bool
  , initialWindowSize :: Size
  , initialWindowPosition :: Position
  , fps :: Natural
  , timeResolution :: Double
  , traceDuration :: Double
  , controlKeys :: [Char]
  , initialControls :: Vector Control
  , v1name :: Text
  , v2name :: Text
  , v3name :: Text
  , phi1name :: Text
  , phi2name :: Text
  , phi3name :: Text
  , lightBname :: Text
  , lightCname :: Text
  , lightDname :: Text
  , lightEname :: Text
  , lightFname :: Text
  , lightGname :: Text
  }
  deriving (Generic, Show, Eq)

data Size = Size {width :: Natural, height :: Natural}
  deriving (Generic, Show, Eq)
data Position = Position {x :: Natural, y :: Natural}
  deriving (Generic, Show, Eq)

instance Interpret Size
instance Interpret Position
instance Interpret Config

getConfig :: Options -> IO Config
getConfig opt = do
  conf <- (Dhall.inputFile Dhall.auto (Options.configFile opt))
  return $ Config
    { csvPath = (Options.csvPath opt) <|> (csvPath conf)
    , hideControls = fromMaybe
        (hideControls conf)
        (Options.hideControls opt)
    , initialWindowSize = fromMaybe
        (initialWindowSize conf)
        ((\(w, h) -> Size w h) <$> Options.initialWindowSize opt)
    , initialWindowPosition = fromMaybe
        (initialWindowPosition conf)
        ((\(x', y') -> Position x' y') <$> Options.initialWindowPosition opt)
    , fps = fromMaybe
        (fps conf)
        (Options.fps opt)
    , timeResolution = fromMaybe
        (timeResolution conf)
        (Options.timeResolution opt)
    , traceDuration = fromMaybe
        (traceDuration conf)
        (Options.traceDuration opt)
    , controlKeys = fromMaybe
        (controlKeys conf)
        (Options.controlKeys opt)
    , initialControls = initialControls conf
    , v1name = v1name conf
    , v2name = v2name conf
    , v3name = v3name conf
    , phi1name = phi1name conf
    , phi2name = phi2name conf
    , phi3name = phi3name conf
    , lightBname = lightBname conf
    , lightCname = lightCname conf
    , lightDname = lightDname conf
    , lightEname = lightEname conf
    , lightFname = lightFname conf
    , lightGname = lightGname conf
    }

