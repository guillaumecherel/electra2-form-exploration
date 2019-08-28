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
    , quadraticToLinear
    , quadraticFromLinear
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

displayMode :: Gloss.Display
displayMode = (Gloss.InWindow "Nice Window" (worldWindow initialWorld) (0, 0))

backgroundColor :: Gloss.Color
backgroundColor = Gloss.black

fps :: Int
fps = 60

-- Control name lowerBound upperBound value
data Control =
  LinearControl
    { controlLabel :: Text
    , controlLowerBound :: Double
    , controlUpperBound :: Double
    , controlValue :: Double }
  | QuadraticControl
    { controlLabel :: Text
    , controlRadius :: Double
    , controlValue :: Double }
  deriving (Show, Eq)

quadraticToLinear :: Double -> Double -> Double
quadraticToLinear radius value =
  let value' = if value >= 0
                  then sqrt (value * radius)
                  else - sqrt (-value * radius)
  in bounded (- radius) radius value'

quadraticFromLinear :: Double -> Double -> Double
quadraticFromLinear radius value =
  let value' = if value >= 0
                  then value ** 2 / radius
                  else - value ** 2 / radius
  in bounded (-radius) radius value'

bounded :: Double -> Double -> Double -> Double
bounded lower upper x = min upper $ max lower x

controlSpecs :: Control -> (Text, Double, Double, Double)
controlSpecs (LinearControl name lower upper val) = (name, lower, upper, val)
controlSpecs (QuadraticControl name radius val) =
  (name, -radius, radius, quadraticToLinear radius val)

toInput :: Vector Control -> Model.Input
toInput controls =
  let values = controlValue <$> controls
  in Model.inputDefault
       (values Vector.! 0)
       (values Vector.! 1)
       (values Vector.! 2)
       (values Vector.! 3)
       (values Vector.! 4)
       (values Vector.! 5)

keyControl :: Char -> Maybe Int
keyControl k = List.elemIndex k
  ['u', 'i', 'e'
  ,'y', 'x', '.']

data World = World
  { worldWindow :: (Int, Int)
  , worldControls :: Vector Control
  , worldTime :: Double
  , worldTrajectories :: [[(Double, Double)]]
  , worldAngles :: Model.Angles
  , worldLayout :: GUI.Layout
  , worldMouseGrabControl :: Maybe Int
  , worldKeyboardGrabControl :: Set Int
  , worldMousePos :: (Float, Float)
  } deriving (Eq, Show)

initialWorld :: World
initialWorld =
  let controls = Vector.fromList
        [ QuadraticControl "v1" 5 1
        , QuadraticControl "v2" 5 (-2)
        , QuadraticControl "v3" 5 0.5
        , LinearControl "phi1" 0 (2 * pi) 0
        , LinearControl "phi2" 0 (2 * pi) 0
        , LinearControl "phi3" 0 (2 * pi) 0
        ]
      window = (400, 600)
  in World
  { worldWindow = window
  , worldControls = controls
  , worldTime = 0
  , worldTrajectories = replicate 8 []
  , worldAngles = Model.Angles 0 0 0
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
  -- trace (show event :: Text) $
  case event of
    Gloss.EventResize (width, height) -> [ResizeWindow (width, height)]
    Gloss.EventKey (Gloss.MouseButton Gloss.LeftButton) Gloss.Up _ _ ->
      case worldMouseGrabControl world of
        Nothing -> []
        Just i -> [MouseReleaseControl i]
    Gloss.EventKey (Gloss.MouseButton Gloss.LeftButton) Gloss.Down _ (x, y) ->
      case GUI.layoutQuery (worldLayout world) (x, y) of
        GUI.SliderAnswer i _ -> [MouseGrabControl i]
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
              let dy = y - snd (worldMousePos world)
                  h = slidersHeight Vector.! i
                  dragAmount = float2Double dy / h
              in DragControl i dragAmount)
          setMousePos = SetMousePos (x, y)
      in setMousePos : dragControls
    _ -> []

updateEvent :: Event -> World -> World
updateEvent event world =
  -- trace (show $ event :: Text) $
  case event of
    ResizeWindow (width, height) ->
      world{ worldWindow = (width, height) }
    MouseGrabControl i -> world
      {worldMouseGrabControl = Just i}
    MouseReleaseControl _ -> world
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
      let newCtrl = case worldControls world Vector.! i of
            (LinearControl n l u v) ->
              (LinearControl n l u (bounded l u $ v + (u - l) * amount))
            (QuadraticControl n r v) ->
              let x = quadraticToLinear r v
                  v' = quadraticFromLinear r (x + amount * 2 * r)
              in QuadraticControl n r v'
      in world { worldControls = worldControls world Vector.// [(i, newCtrl)] }

updateInputs :: Gloss.Event -> World -> World
updateInputs event world =
  -- traceShow event $
  (foldl' (.) identity
  $ fmap updateEvent (events event world))
  $ world


updateTime :: Float -> World -> World
updateTime dt world =
  let newAngles = Model.angles
        (Model.inputSetAngles input (worldAngles world))
        (float2Double dt)
      newTrajectories =
        (fmap . fmap) (bimap scale scale)
        $ Model.trajectoriesToList
        $ Model.trajectories inputShiftAngles 0 tResolution (float2Double dt)
      inputShiftAngles = Model.inputAddAngles input newAngles
      tResolution = 0.005
      t = worldTime world + float2Double dt
      scale x = x / (2 * Model.span input)
      trajectories = zipWith
        (\new prev ->
          take (ceiling $ traceDuration / tResolution)
            (new <> prev))
        newTrajectories (worldTrajectories world)
      traceDuration = 1 -- 1.0 / 25.0
      input = toInput ctrl
      ctrl = worldControls world
  in 
     world
     { worldTime = t
     , worldTrajectories = trajectories
     , worldAngles = newAngles
     , worldLayout = GUI.layoutWithControl
        (worldWindow world)
        (controlSpecs <$> Vector.toList ctrl)
        ((fmap . fmap) (bimap double2Float double2Float) trajectories)
     }
