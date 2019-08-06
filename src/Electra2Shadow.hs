{-# LANGUAGE NoImplicitPrelude #-}

module Electra2Shadow
    ( displayMode
    , backgroundColor
    , fps
    , initialWorld
    , view
    , updateInputs
    , updateTime
    ) where

import Protolude

import GHC.Float (double2Float)
import qualified Graphics.Gloss as Gloss
import qualified Graphics.Gloss.Interface.IO.Interact as Gloss

displayMode = (Gloss.InWindow "Nice Window" (200, 200) (10, 10))

backgroundColor = Gloss.black

fps :: Int
fps = 60

data World = World
  { worldWidth :: Int
  , worldHeight :: Int
  , worldV1 :: Double
  , worldV2 :: Double
  , worldV3 :: Double
  , worldA1 :: Double
  , worldA2 :: Double
  , worldA3 :: Double
  } deriving (Eq, Show)

initialWorld :: World
initialWorld = World
  { worldV1 = 0.0
  , worldV2 = 0.0
  , worldV3 = 0.0
  , worldA1 = 0.0
  , worldA2 = 0.0
  , worldA3 = 0.0
  }

view :: World -> Gloss.Picture
view w = Gloss.pictures
  [ position 0 (-0.5) w $ viewCommands w
  , position 0 (0.5) w $ viewModel w]

viewCommands :: World -> Gloss.Picture
viewCommands w = Gloss.color Gloss.white $ Gloss.circle 25

viewModel :: World -> Gloss.Picture
viewModel w = Gloss.color Gloss.white $ Gloss.circle 5

position :: Double -> Double -> World -> Gloss.Picture -> Gloss.Picture
position x y world = Gloss.translate moveX moveY
  where moveX = double2Float $ (0.5 - x) * fromIntegral (worldWidth world)
        moveY = double2Float $ (0.5 - y) * fromIntegral (worldHeight world)

updateInputs :: Gloss.Event -> World -> World
updateInputs event world = case event of
  Gloss.EventResize (height,width) ->
    world{worldWidth = width
         , worldHeight= height}
  _ -> world

updateTime :: Float -> World -> World
updateTime dt w = w

