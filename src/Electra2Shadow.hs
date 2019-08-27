{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TupleSections #-}
{-# LANGUAGE NoImplicitPrelude #-}

module Electra2Shadow
    ( displayMode
    , backgroundColor
    , fps
    , initialWorld
    , updateInputs
    , updateTime
    , view
    ) where

import Protolude

import qualified Data.List as List
import qualified Data.Set as Set
import           Data.Set (Set)
import qualified Data.Vector as Vector
import           Data.Vector (Vector)
import           GHC.Float (double2Float, float2Double)
import qualified Graphics.Gloss as Gloss
import qualified Graphics.Gloss.Interface.IO.Interact as Gloss

import qualified Electra2Shadow.GUI as GUI
import qualified Electra2Shadow.Model as Model

displayMode = (Gloss.InWindow "Nice Window" (worldWindow initialWorld) (0, 0))

backgroundColor = Gloss.black

fps :: Int
fps = 60

-- Control name lowerBound upperBound value
data Control = Control Text Double Double Double
  deriving (Show, Eq)

controlSpecs :: Control -> (Text, Double, Double, Double)
controlSpecs (Control name lower upper val) = (name, lower, upper, val)

toInput :: Vector Control -> Model.Input
toInput controls = Model.inputDefault v1 v2 v3 phi1 phi2 phi3
  where [v1, v2, v3, phi1, phi2, phi3] =
          fmap (\(Control _ _ _ v) -> v)$ Vector.toList controls

keyControl :: Char -> Maybe Int
keyControl k = List.elemIndex k
  ['u', 'i', 'e'
  ,'y', 'x', '.']

data World = World
  { worldWindow :: (Int, Int)
  , worldControls :: Vector Control
  , worldTime :: Double
  , worldTrajectories :: [[(Double, Double)]]
  , worldLayout :: GUI.Layout
  , worldMouseGrabControl :: Maybe Int
  , worldKeyboardGrabControl :: Set Int
  , worldMousePos :: (Float, Float)
  } deriving (Eq, Show)

initialWorld :: World
initialWorld =
  let controls = Vector.fromList
        [ Control "v1" (-5) 5 1
        , Control "v2" (-5) 5 (-2)
        , Control "v3" (-5) 5 0.5
        , Control "phi1" 0 (2 * pi) 0
        , Control "phi2" 0 (2 * pi) 0
        , Control "phi3" 0 (2 * pi) 0
        ]
      window = (400, 600)
  in World
  { worldWindow = window
  , worldControls = controls
  , worldTime = 0
  , worldTrajectories = []
  , worldLayout = GUI.layoutWithControl
      window
      (fmap controlSpecs $ Vector.toList controls)
      []
  , worldMouseGrabControl = Nothing
  , worldKeyboardGrabControl = mempty
  , worldMousePos = (0, 0)
  }

view :: World -> Gloss.Picture
view world =
  Gloss.pictures
    [ GUI.viewLayout (worldLayout world)
    , Gloss.color Gloss.white
      $ Gloss.rectangleWire
          (fromIntegral $ fst (worldWindow world) - 1)
          (fromIntegral $ snd (worldWindow world) - 1)
    ]

data Event =
  ResizeWindow (Int, Int)
  | DragControl Int Double
  | MouseGrabControl Int
  | KeyboardGrabControl Int
  | MouseReleaseControl Int
  | KeyboardReleaseControl Int
  | SetMousePos (Float, Float)
  deriving (Show, Eq)

events :: Gloss.Event -> World -> [Event]
events event world =
  trace (show event :: Text) $
  case event of
    Gloss.EventResize (width, height) -> [ResizeWindow (width, height)]
    Gloss.EventKey (Gloss.MouseButton Gloss.LeftButton) Gloss.Up _ (x, y) ->
      case worldMouseGrabControl world of
        Nothing -> []
        Just i -> [MouseReleaseControl i]
    Gloss.EventKey (Gloss.MouseButton Gloss.LeftButton) Gloss.Down _ (x, y) ->
      case GUI.layoutQuery (worldLayout world) (x, y) of
        GUI.SliderAnswer i v -> [MouseGrabControl i]
        GUI.NoAnswer -> []
    Gloss.EventKey (Gloss.Char k) Gloss.Down _ _ ->
      toList $ KeyboardGrabControl <$> keyControl k
    Gloss.EventKey (Gloss.Char k) Gloss.Up _ _ ->
      toList $ KeyboardReleaseControl <$> keyControl k
    Gloss.EventMotion (x, y) ->
      let grabbedControls =
            worldKeyboardGrabControl world
            <> Set.fromList (toList $ worldMouseGrabControl world)
          slidersHeight = fmap float2Double $ Vector.fromList
            $ GUI.slidersHeight (worldLayout world)
          dragControls = flip fmap
            (Set.toList grabbedControls)
            (\i ->
              let dx = x - fst (worldMousePos world)
                  dy = y - snd (worldMousePos world)
                  (Control _ l u v) = worldControls world Vector.! i
                  h = slidersHeight Vector.! i
                  dragAmount = float2Double dy / h
              in DragControl i dragAmount)
          setMousePos = SetMousePos (x, y)
      in setMousePos : dragControls
    _ -> []

updateEvent :: Event -> World -> World
updateEvent event world =
  trace (show $ event :: Text) $
  case event of
    ResizeWindow (width, height) ->
      world{ worldWindow = (width, height) }
    MouseGrabControl i -> world
      {worldMouseGrabControl = Just i}
    MouseReleaseControl i -> world
      { worldMouseGrabControl = Nothing }
    KeyboardGrabControl i -> world
      { worldKeyboardGrabControl =
          Set.insert i $ worldKeyboardGrabControl world
      }
    KeyboardReleaseControl i -> world
      { worldKeyboardGrabControl =
          Set.delete i $ worldKeyboardGrabControl world
      }
    SetMousePos pos -> world { worldMousePos = pos }
    DragControl i amount ->
      let (Control n l u v) = worldControls world Vector.! i
          v' = (max l $ min u $ v + (u - l) * amount)
      in world
        { worldControls = worldControls world Vector.//
           [(i, Control n l u v')]
        }

updateInputs :: Gloss.Event -> World -> World
updateInputs event world =
  traceShow event $
  (foldl' (.) identity
  $ fmap updateEvent (events event world))
  $ world


updateTime :: Float -> World -> World
updateTime dt world =
  world{ worldTime = t
       , worldTrajectories = trajectories
       , worldLayout = GUI.layoutWithControl
          (worldWindow world)
          (controlSpecs <$> Vector.toList ctrl)
          ((fmap . fmap) (bimap double2Float double2Float) trajectories)
       }
  where t = worldTime world + float2Double dt
        trajectories =
           (fmap . fmap) (bimap scale scale)
           $ Model.trajectoriesToList
           $ Model.trajectories input tStart tResolution t
        scale x = x / (2 * Model.span input)
        input = toInput ctrl
        ctrl = worldControls world
        tStart = t - 1/25
        tResolution = 0.001
