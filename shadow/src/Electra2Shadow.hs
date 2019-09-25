{-# LANGUAGE DeriveGeneric     #-}
{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TupleSections #-}
{-# LANGUAGE NoImplicitPrelude #-}

module Electra2Shadow
    ( displayMode
    , backgroundColor
    , initialWorld
    , updateInputs
    , updateTime
    , view
    ) where

import Protolude 

import           Data.Char (isDigit)
import qualified Data.Set as Set
import           Data.Set (Set)
import qualified Data.Text as Text
import           Data.Text.Read (double)
import qualified Data.Vector as Vector
import           GHC.Float (double2Float, float2Double)
import qualified Graphics.Gloss as Gloss
import qualified Graphics.Gloss.Interface.IO.Interact as Gloss

import qualified Electra2Shadow.Config as Config
import           Electra2Shadow.Config (Config)
import qualified Electra2Shadow.Control as Control
import           Electra2Shadow.Control (Controls)
import qualified Electra2Shadow.GUI as GUI
import qualified Electra2Shadow.Model as Model

displayMode :: Config -> Gloss.Display
displayMode options =
  (Gloss.InWindow "Nice Window"
    (case Config.initialWindowSize options of
       Config.Size {Config.width = w, Config.height = h} ->
         (fromIntegral w, fromIntegral h))
    (case Config.initialWindowPosition options of
       Config.Position {Config.x = x, Config.y = y} ->
         (fromIntegral x, fromIntegral y)))

backgroundColor :: Gloss.Color
backgroundColor = Gloss.black

data World = World
  { worldWindow :: (Int, Int)
  , worldControls :: Controls
  , worldTime :: Double
  , worldTrajectories :: [[(Double, Double)]]
  , worldAngles :: Model.Angles
  , worldLayout :: GUI.Layout
  , worldMouseGrabControl :: Maybe Int
  , worldKeyboardGrabControl :: Set Int
  , worldMousePos :: (Float, Float)
  , worldNumberBuffer :: [Char]
  } deriving (Show)

getWorldWindow :: ((Int, Int) -> a) -> World -> a
getWorldWindow f w = f $  worldWindow w
getWorldControls :: (Controls -> a) -> World -> a
getWorldControls f w = f $ worldControls w
getWorldTime :: (Double -> a) -> World -> a
getWorldTime f w = f $  worldTime w
getWorldTrajectories :: ([[(Double, Double)]] -> a) -> World -> a
getWorldTrajectories f w = f $  worldTrajectories w
getWorldAngles :: (Model.Angles -> a) -> World -> a
getWorldAngles f w = f $  worldAngles w
getWorldLayout :: (GUI.Layout -> a) -> World -> a
getWorldLayout f w = f $  worldLayout w
getWorldMouseGrabControl :: (Maybe Int -> a) -> World -> a
getWorldMouseGrabControl f w = f $  worldMouseGrabControl w
getWorldKeyboardGrabControl :: (Set Int -> a) -> World -> a
getWorldKeyboardGrabControl f w = f $  worldKeyboardGrabControl w
getWorldMousePos :: ((Float, Float) -> a) -> World -> a
getWorldMousePos f w = f $  worldMousePos w
getWorldNumberBuffer :: ([Char] -> a) -> World -> a
getWorldNumberBuffer f w = f $  worldNumberBuffer w

setWorldWindow :: ((Int, Int) -> (Int, Int)) -> World -> World
setWorldWindow f w = w { worldWindow = f $ worldWindow w }
setWorldControls :: (Controls -> Controls) -> World -> World
setWorldControls f w = w { worldControls = f $ worldControls w }
setWorldMouseGrabControl :: (Maybe Int -> Maybe Int) -> World -> World
setWorldMouseGrabControl f w =
  w { worldMouseGrabControl = f $ worldMouseGrabControl w }
setWorldKeyboardGrabControl :: (Set Int -> Set Int) -> World -> World
setWorldKeyboardGrabControl f w =
  w { worldKeyboardGrabControl = f $ worldKeyboardGrabControl w }
setWorldMousePos :: ((Float, Float) -> (Float, Float)) -> World -> World
setWorldMousePos f w = w { worldMousePos = f $ worldMousePos w }
setWorldNumberBuffer :: ([Char] -> [Char]) -> World -> World
setWorldNumberBuffer f w = w { worldNumberBuffer = f $ worldNumberBuffer w }

grabbedControls :: World -> Set Int
grabbedControls world =
  getWorldKeyboardGrabControl identity world
  <> Set.fromList (toList $ getWorldMouseGrabControl identity world)

initialWorld :: Config -> Controls -> World
initialWorld options controls =
  let windowSize = case displayMode options of
        (Gloss.InWindow _ xy _) -> xy
        Gloss.FullScreen -> case Config.initialWindowSize options of
          Config.Size {Config.width = w, Config.height = h} ->
            (fromIntegral w, fromIntegral h)
  in World
      { worldWindow = windowSize
      , worldControls = controls
      , worldTime = 0
      , worldTrajectories = replicate 8 []
      , worldAngles = Model.Angles 0 0 0
      , worldLayout = GUI.layout
          (Config.hideControls options)
          windowSize
          (zip (Control.specs <$> Control.toList controls)
               (const 0 <$> Control.toList controls))
          []
          []
      , worldMouseGrabControl = Nothing
      , worldKeyboardGrabControl = mempty
      , worldMousePos = (0, 0)
      , worldNumberBuffer = []
      }

view :: World -> Gloss.Picture
view world =
  Gloss.pictures
    [ GUI.viewLayout (getWorldLayout identity world)
    , Gloss.color Gloss.white
      $ Gloss.rectangleWire
          (fromIntegral $ fst (getWorldWindow identity world) - 1)
          (fromIntegral $ snd (getWorldWindow identity world) - 1)
    ]

data Event =
  ResizeWindow (Int, Int)
  | MouseGrabControl Int
  | KeyboardGrabControl Int
  | DragControl Int Double
  | MouseReleaseControl Int
  | KeyboardReleaseControl Int
  | SetMousePos (Float, Float)
  | ReadDigit Char
  | SetControl Int Double
  deriving (Show, Eq)

events :: Config -> Gloss.Event -> World -> [Event]
events options event world =
  -- trace (show event :: Text) $
  case event of
    Gloss.EventResize (width, height) -> [ResizeWindow (width, height)]
    Gloss.EventKey (Gloss.MouseButton Gloss.LeftButton) Gloss.Up _ _ ->
      case getWorldMouseGrabControl identity world of
        Nothing -> []
        Just i -> [MouseReleaseControl i]
    Gloss.EventKey (Gloss.MouseButton Gloss.LeftButton) Gloss.Down _ (x, y) ->
      case GUI.layoutQuery (getWorldLayout identity world) (x, y) of
        GUI.SliderAnswer i r ->
          let newValue = case (getWorldControls . Control.getControlAt i) identity world of
                Just (Control.LinearControl _ l u _) ->
                  Control.bounded l u
                  $ r * (u - l)
                Just (Control.QuadraticControl _ l c u _) ->
                  Control.bounded l u
                  $ Control.quadraticFromLinear l c u (r * (u - l))
                Nothing -> panic $ "events: No control at index " <> show i
          in [MouseGrabControl i, SetControl i newValue ]
        GUI.NoAnswer -> []
    Gloss.EventKey (Gloss.Char k) Gloss.Down _ _ ->
      if isDigit k || elem k ['.', '-', '+']
        then [ReadDigit k]
        else
          toList $ KeyboardGrabControl
          <$> Control.keyControl (Config.controlKeys options) k
    Gloss.EventKey (Gloss.Char k) Gloss.Up _ _ ->
      if isDigit k || elem k ['.', '-', '+']
        then [] 
        else case Control.keyControl (Config.controlKeys options) k of
          Nothing -> []
          Just i ->
            let value = double (Text.pack $ reverse $ getWorldNumberBuffer identity world)
                setctrls = case value of
                  Right (n, _) -> [SetControl i n]
                  Left err -> [] &
                    trace ("getWorldNumberBuffer cannot be read as a number: "
                      <> show (getWorldNumberBuffer identity world)
                      <> " " <> Text.pack err :: Text)
                releases = [KeyboardReleaseControl i]
            in setctrls <> releases
    Gloss.EventMotion (x, y) ->
      let slidersHeight = fmap float2Double $ Vector.fromList
            $ GUI.slidersHeight (getWorldLayout identity world)
          -- Controls grabbed with the keyboard get dragged more slowly than those grabbed with the mouse, allowing higher accuracy.
          kDragControls = flip fmap
            (Set.toList $ getWorldKeyboardGrabControl identity world)
            (\i ->
              let dy = y - snd (getWorldMousePos identity world)
                  h = fromIntegral $ snd $ worldWindow world
                  dragAmount = float2Double dy / h
              in DragControl i dragAmount)
          mDragControls = flip fmap
            (toList $ getWorldMouseGrabControl identity world)
            (\i ->
              let dy = y - snd (getWorldMousePos identity world)
                  h = slidersHeight Vector.! i
                  dragAmount = float2Double dy / h
              in DragControl i dragAmount)
          setMousePos = SetMousePos (x, y)
      in setMousePos : kDragControls <> mDragControls
    _ -> []

updateEvent :: Event -> World -> World
updateEvent event world = 
  -- traceShow event $
  case event of
    ResizeWindow size -> setWorldWindow (const size) world
    MouseGrabControl i ->  setWorldMouseGrabControl (const (Just i)) world
    MouseReleaseControl _ -> setWorldMouseGrabControl (const Nothing) world
    KeyboardGrabControl i ->  setWorldKeyboardGrabControl (Set.insert i) world
    KeyboardReleaseControl i ->
      let keyGrab = Set.delete i $ getWorldKeyboardGrabControl identity world
          numBuf = if null keyGrab
                      then mempty
                      else (getWorldNumberBuffer identity world)
      in setWorldKeyboardGrabControl (const keyGrab)
         $ setWorldNumberBuffer (const numBuf)
         $ world
    SetMousePos pos -> setWorldMousePos (const pos) world
    DragControl i amount ->
      let newCtrl = case (getWorldControls . Control.getControlAt i) identity world of
            Just (Control.LinearControl n l u v) ->
              (Control.LinearControl n l u (Control.bounded l u $ v + (u - l) * amount))
            Just (Control.QuadraticControl n l c u v) ->
              let x = Control.quadraticToLinear l c u v
                  v' = Control.quadraticFromLinear l c u (x + amount * (u - l))
              in Control.QuadraticControl n l c u v'
            Nothing -> panic $ "updateEvent: No control at index " <> show i
      in (setWorldControls . Control.setControlAt i) (const newCtrl) world
    ReadDigit d -> setWorldNumberBuffer (d :) world
    SetControl i n ->
      (setWorldControls . Control.setControlAt i . Control.setControlValue) (const n) world

updateInputs :: Config -> Gloss.Event -> World -> World
updateInputs options event world =
  -- traceShow event $
  (foldl' (.) identity
  $ fmap updateEvent (events options event world))
  $ world


updateTime :: Config -> Float -> World -> World
updateTime options dt world =
  let tResolution = Config.timeResolution options
      traceDuration = Config.traceDuration options
      newAngles = Model.angles
        (Model.inputSetAngles input (getWorldAngles identity world))
        (float2Double dt)
      newTrajectories =
        (fmap . fmap) (bimap scale scale)
        $ Model.trajectoriesToList
        $ Model.trajectories inputShiftAngles 0 tResolution (float2Double dt)
      inputShiftAngles = Model.inputAddAngles input newAngles
      t = getWorldTime identity world + float2Double dt
      scale x = x / (2 * Model.span input)
      trajectories = zipWith
        (\new prev ->
          take (ceiling $ traceDuration / tResolution)
            (new <> prev))
        newTrajectories (getWorldTrajectories identity world)
      (nearestCtrlVals, input) = Control.toInput ctrls
      ctrls = getWorldControls identity world
      lightIntensities = [ Model.parameterValue $ Model.parametersLightB input
                         , Model.parameterValue $ Model.parametersLightC input
                         , Model.parameterValue $ Model.parametersLightD input
                         , Model.parameterValue $ Model.parametersLightE input
                         , Model.parameterValue $ Model.parametersLightF input
                         , Model.parameterValue $ Model.parametersLightG input
                         ]
      speeds :: [(Text, Double)]
      speeds = ((,) <$> Model.parameterName <*> Model.parameterValue)
               <$> [ Model.parametersV1 input
                   , Model.parametersV2 input
                   , Model.parametersV3 input
                   ]
  in 
     world
     { worldTime = t
     , worldTrajectories = trajectories
     , worldAngles = newAngles
     , worldLayout = GUI.layout
         (Config.hideControls options)
         (getWorldWindow identity world)
         (zip (Control.specs <$> Control.toList ctrls)
              (Control.controlsValuesList nearestCtrlVals))
         speeds
         (zip lightIntensities
           $ (fmap . fmap) (bimap double2Float double2Float) trajectories)
     }
