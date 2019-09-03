{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TupleSections #-}
{-# LANGUAGE NoImplicitPrelude #-}

module Electra2Shadow
    ( displayMode
    , backgroundColor
    , controlsFromInputValues
    , controlsFromMap
    , ModelInputValues (..)
    , ControlsValues (..)
    , controlMapFromCSV
    , initialWorld
    , controlsToInput
    , updateInputs
    , updateTime
    , view
    , quadraticToLinear
    , quadraticFromLinear
    ) where

import Protolude 

import qualified Data.ByteString.Lazy as BL
import           Data.Char (isDigit)
import qualified Data.Csv as CSV
import qualified Data.List as List
import qualified Data.KdMap.Static as KDM
import qualified Data.Set as Set
import           Data.Set (Set)
import qualified Data.Text as Text
import           Data.Text.Read (double)
import qualified Data.Vector as Vector
import           GHC.Float (double2Float, float2Double)
import qualified Graphics.Gloss as Gloss
import qualified Graphics.Gloss.Interface.IO.Interact as Gloss

import qualified Electra2Shadow.Options as Options
import           Electra2Shadow.Options (Options)
import qualified Electra2Shadow.GUI as GUI
import qualified Electra2Shadow.Model as Model

displayMode :: Options -> Gloss.Display
displayMode options =
  (Gloss.InWindow "Nice Window"
    (Options.initialWindowSize options)
    (Options.initialWindowPosition options))

backgroundColor :: Gloss.Color
backgroundColor = Gloss.black

-- Control name lowerBound upperBound value
data Control =
    LinearControl Text Double Double Double 
    -- label lowerBound upperBound value
  | QuadraticControl Text Double Double
    -- label radius value
  deriving (Show, Eq)

getControlValue :: (Double -> a) -> Control -> a
getControlValue f (LinearControl _ _ _ v) = f v
getControlValue f (QuadraticControl _ _ v) = f v

setControlValue :: (Double -> Double) -> Control -> Control
setControlValue f (LinearControl n l u v) = (LinearControl n l u (f v))
setControlValue f (QuadraticControl n r v) = (QuadraticControl n r (f v))

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

data Controls =
    ControlInput Control Control Control Control Control Control
    -- v1 v2 v3 phi1 phi2 phi3
  | ControlForm (KDM.KdMap Double ControlsValues ModelInputValues)
      Control Control Control Control Control
    -- ControlsInputsMap NumberSingular Density Moran MeanCurvature ReturnTime
  deriving (Show)


getControlsV1 :: (Control -> a) -> Controls -> Maybe a
getControlsV1 f (ControlInput v1 _ _ _ _ _) = Just $ f v1
getControlsV1 _ _ = Nothing

getControlsV2 :: (Control -> a) ->  Controls -> Maybe a
getControlsV2 f (ControlInput _ v2 _ _ _ _) = Just $ f v2
getControlsV2 _ _ = Nothing

getControlsV3 :: (Control -> a) -> Controls -> Maybe a
getControlsV3 f (ControlInput _ _ v3 _ _ _) = Just $ f v3
getControlsV3 _ _ = Nothing

getControlsPhi1 :: (Control -> a) -> Controls -> Maybe a
getControlsPhi1 f (ControlInput _  _  _ phi1 _  _) = Just $ f phi1
getControlsPhi1 _ _ = Nothing

getControlsPhi2 :: (Control -> a) -> Controls -> Maybe a
getControlsPhi2 f (ControlInput _  _  _  _ phi2 _) = Just $ f phi2
getControlsPhi2 _ _ = Nothing

getControlsPhi3 :: (Control -> a) -> Controls -> Maybe a
getControlsPhi3 f (ControlInput _  _  _  _  _ phi3) = Just $ f phi3
getControlsPhi3 _ _ = Nothing

getControlsNumberSingular :: (Control -> a) -> Controls -> Maybe a
getControlsNumberSingular f (ControlForm _ ns _  _  _  _) = Just $ f ns
getControlsNumberSingular _ _ = Nothing

getControlsDensity :: (Control -> a) -> Controls -> Maybe a
getControlsDensity f (ControlForm _  _ d _  _  _) = Just $ f d
getControlsDensity _ _ = Nothing

getControlsMoran :: (Control -> a) -> Controls -> Maybe a
getControlsMoran f (ControlForm _  _  _ m _  _) = Just $ f m
getControlsMoran _ _ = Nothing

getControlsMeanCurvature :: (Control -> a) -> Controls -> Maybe a
getControlsMeanCurvature f (ControlForm _  _  _  _ mc _) = Just $ f mc
getControlsMeanCurvature _ _ = Nothing

getControlsReturnTime :: (Control -> a) -> Controls -> Maybe a
getControlsReturnTime f (ControlForm _  _  _  _  _ rt) = Just $ f rt
getControlsReturnTime _ _ = Nothing

setControlsV1 :: (Control -> Control) -> Controls -> Controls
setControlsV1 f (ControlInput v1 v2 v3 phi1 phi2 phi3) =
 ControlInput (f v1) v2 v3 phi1 phi2 phi3
setControlsV1 _ c = c

setControlsV2 :: (Control -> Control) ->  Controls -> Controls
setControlsV2 f (ControlInput v1 v2 v3 phi1 phi2 phi3) =
 ControlInput v1 (f v2) v3 phi1 phi2 phi3
setControlsV2 _ c = c

setControlsV3 :: (Control -> Control) ->  Controls -> Controls
setControlsV3 f (ControlInput v1 v2 v3 phi1 phi2 phi3) =
 ControlInput v1 v2 (f v3) phi1 phi2 phi3
setControlsV3 _ c = c

setControlsPhi1 :: (Control -> Control) ->  Controls -> Controls
setControlsPhi1 f (ControlInput v1 v2 v3 phi1 phi2 phi3) =
 ControlInput v1 v2 v3 (f phi1) phi2 phi3
setControlsPhi1 _ c = c

setControlsPhi2 :: (Control -> Control) ->  Controls -> Controls
setControlsPhi2 f (ControlInput v1 v2 v3 phi1 phi2 phi3) =
 ControlInput v1 v2 v3 phi1 (f phi2) phi3
setControlsPhi2 _ c = c

setControlsPhi3 :: (Control -> Control) ->  Controls -> Controls
setControlsPhi3 f (ControlInput v1 v2 v3 phi1 phi2 phi3) =
 ControlInput v1 v2 v3 phi1 phi2 (f phi3)
setControlsPhi3 _ c = c

setControlsNumberSingular :: (Control -> Control) ->  Controls -> Controls
setControlsNumberSingular f (ControlForm kdm ns d m mc rt) =
 ControlForm kdm (f ns) d m mc rt
setControlsNumberSingular _ c = c

setControlsDensity :: (Control -> Control) ->  Controls -> Controls
setControlsDensity f (ControlForm kdm ns d m mc rt) =
 ControlForm kdm ns (f d) m mc rt
setControlsDensity _ c = c

setControlsMoran :: (Control -> Control) ->  Controls -> Controls
setControlsMoran f (ControlForm kdm ns d m mc rt) =
 ControlForm kdm ns d (f m) mc rt
setControlsMoran _ c = c

setControlsMeanCurvature :: (Control -> Control) ->  Controls -> Controls
setControlsMeanCurvature f (ControlForm kdm ns d m mc rt) =
 ControlForm kdm ns d m (f mc) rt
setControlsMeanCurvature _ c = c

setControlsReturnTime :: (Control -> Control) -> Controls -> Controls
setControlsReturnTime f (ControlForm kdm ns d m mc rt) =
 ControlForm kdm ns d m mc (f rt)
setControlsReturnTime _ c = c

controlsFromInputValues
  :: Double -> Double -> Double -> Double -> Double -> Double
  -> Controls
controlsFromInputValues v1 v2 v3 phi1 phi2 phi3 = ControlInput
  (QuadraticControl "v1" 5 v1)
  (QuadraticControl "v2" 5 v2)
  (QuadraticControl "v3" 5 v3)
  (LinearControl "phi1" 0 (2 * pi) phi1)
  (LinearControl "phi2" 0 (2 * pi) phi2)
  (LinearControl "phi3" 0 (2 * pi) phi3)

controlsFromMap
  :: [(ControlsValues, ModelInputValues)]
  -> Double -> Double -> Double -> Double -> Double
  -> Controls
controlsFromMap cim ns d m mc rt =
  ControlForm
    (KDM.build (controlsValuesList) cim)
    (LinearControl "singular" 0 10 ns)
    (LinearControl "density" 0 1 d)
    (LinearControl "moran" (-1) 1 m)
    (QuadraticControl "curv" 50 mc)
    (LinearControl "return time" 0 10 rt)

controlMapFromCSV :: BL.ByteString -> ([Text], [(ControlsValues, ModelInputValues)])
controlMapFromCSV csv = case CSV.decodeByName csv of
  Left err -> panic $ "Could not decode csv data.\n" <> show err
  Right (header, result) ->
    let names = decodeUtf8 <$> Vector.toList header
    in (names, Vector.toList result)

controlsValues :: Controls -> ControlsValues
controlsValues (ControlInput v1 v2 v3 phi1 phi2 phi3) =
  InputValues $ ModelInputValues
    (getControlValue identity v1)
    (getControlValue identity v2)
    (getControlValue identity v3)
    (getControlValue identity phi1)
    (getControlValue identity phi2)
    (getControlValue identity phi3)
controlsValues (ControlForm _ ns d m mc rt) =
  FormValues
    (getControlValue identity ns)
    (getControlValue identity d)
    (getControlValue identity m)
    (getControlValue identity mc)
    (getControlValue identity rt)

controlsValuesList :: ControlsValues -> [Double]
controlsValuesList (InputValues (ModelInputValues v1 v2 v3 phi1 phi2 phi3)) =
  [v1, v2, v3, phi1, phi2, phi3]
controlsValuesList (FormValues ns d m mc rt) = [ns, d, m, mc, rt]

controlsList :: Controls -> [Control]
controlsList (ControlInput v1 v2 v3 phi1 phi2 phi3) =
  [ v1, v2, v3, phi1, phi2, phi3 ]
controlsList (ControlForm _ ns d m mc rt) =
  [ ns, d, m, mc, rt ]
 
getControlAt :: Int -> (Control -> a) -> Controls -> Maybe a
getControlAt 0 f ctrls@ControlInput{} = getControlsV1 f ctrls
getControlAt 1 f ctrls@ControlInput{} = getControlsV2 f ctrls
getControlAt 2 f ctrls@ControlInput{} = getControlsV3 f ctrls
getControlAt 3 f ctrls@ControlInput{} = getControlsPhi1 f ctrls
getControlAt 4 f ctrls@ControlInput{} = getControlsPhi2 f ctrls
getControlAt 5 f ctrls@ControlInput{} = getControlsPhi3 f ctrls
getControlAt 0 f ctrls@ControlForm{} = getControlsNumberSingular f ctrls
getControlAt 1 f ctrls@ControlForm{} = getControlsDensity f ctrls
getControlAt 2 f ctrls@ControlForm{} = getControlsMoran f ctrls
getControlAt 3 f ctrls@ControlForm{} = getControlsMeanCurvature f ctrls
getControlAt 4 f ctrls@ControlForm{} = getControlsReturnTime f ctrls
getControlAt _ _ _ = Nothing

setControlAt :: Int -> (Control -> Control) -> Controls -> Controls
setControlAt 0 f ctrls@ControlInput{} = setControlsV1 f ctrls
setControlAt 1 f ctrls@ControlInput{} = setControlsV2 f ctrls
setControlAt 2 f ctrls@ControlInput{} = setControlsV3 f ctrls
setControlAt 3 f ctrls@ControlInput{} = setControlsPhi1 f ctrls
setControlAt 4 f ctrls@ControlInput{} = setControlsPhi2 f ctrls
setControlAt 5 f ctrls@ControlInput{} = setControlsPhi3 f ctrls
setControlAt 0 f ctrls@ControlForm{} = setControlsNumberSingular f ctrls
setControlAt 1 f ctrls@ControlForm{} = setControlsDensity f ctrls
setControlAt 2 f ctrls@ControlForm{} = setControlsMoran f ctrls
setControlAt 3 f ctrls@ControlForm{} = setControlsMeanCurvature f ctrls
setControlAt 4 f ctrls@ControlForm{} = setControlsReturnTime f ctrls
setControlAt i _ ctrls = panic
    $ "controlAt: No control at index " <> show i <> " for " <> show ctrls

controlsToInput :: Controls -> Model.Input
controlsToInput (ControlInput v1 v2 v3 phi1 phi2 phi3 ) =
  Model.inputDefault
    (getControlValue identity v1)
    (getControlValue identity v2)
    (getControlValue identity v3)
    (getControlValue identity phi1)
    (getControlValue identity phi2)
    (getControlValue identity phi3)
controlsToInput ctrl@(ControlForm kdm _ _ _ _ _) =
  let (_, ModelInputValues v1 v2 v3 phi1 phi2 phi3) =
        KDM.nearest kdm (controlsValues ctrl)
  in Model.inputDefault v1 v2 v3 phi1 phi2 phi3

data ControlsValues =
    FormValues Double Double Double Double Double
    -- numberSingular density moran meanCurvature returnTime
  | InputValues ModelInputValues
  deriving (Show, Eq)

data ModelInputValues =
  ModelInputValues Double Double Double Double Double Double
  -- (v1, v2, v3, phi1, phi2, phi3)
  deriving (Show, Eq)

instance CSV.FromNamedRecord ControlsValues where
  parseNamedRecord m = FormValues
    <$> m CSV..: "nbPointsSinguliers"
    <*> m CSV..: "densite"
    <*> m CSV..: "moran"
    <*> m CSV..: "courbureMoyenne"
    <*> m CSV..: "nbPointsRetour"

instance CSV.FromNamedRecord ModelInputValues where
  parseNamedRecord m = ModelInputValues
    <$> m CSV..: "v1"
    <*> m CSV..: "v2"
    <*> m CSV..: "v3"
    -- <*> m CSV..: "angleIni_B"
    <*> pure 0
    <*> m CSV..: "angleIni_D"
    <*> m CSV..: "angleIni_F"

instance CSV.FromNamedRecord (ControlsValues, ModelInputValues) where
  parseNamedRecord m = (,)
    <$> CSV.parseNamedRecord m
    <*> CSV.parseNamedRecord m

keyControl :: [Char] -> Char -> Maybe Int
keyControl l k = List.elemIndex k l

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

initialWorld :: Options -> Controls -> World
initialWorld options controls =
  let windowSize = case displayMode options of
        (Gloss.InWindow _ xy _) -> xy
        Gloss.FullScreen -> (Options.initialWindowSize options)
  in World
      { worldWindow = windowSize
      , worldControls = controls
      , worldTime = 0
      , worldTrajectories = replicate 8 []
      , worldAngles = Model.Angles 0 0 0
      , worldLayout = GUI.layout
          (Options.hideControls options)
          windowSize
          (controlSpecs <$> controlsList controls)
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
  | DragControl Int Double
  | MouseGrabControl Int
  | KeyboardGrabControl Int
  | MouseReleaseControl Int
  | KeyboardReleaseControl Int
  | SetMousePos (Float, Float)
  | ReadDigit Char
  | SetControl Int Double
  deriving (Show, Eq)

events :: Options -> Gloss.Event -> World -> [Event]
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
        GUI.SliderAnswer i _ -> [MouseGrabControl i]
        GUI.NoAnswer -> []
    Gloss.EventKey (Gloss.Char k) Gloss.Down _ _ ->
      if isDigit k || elem k ['.', '-', '+']
        then [ReadDigit k]
        else 
          toList $ KeyboardGrabControl
          <$> keyControl (Options.controlKeys options) k
    Gloss.EventKey (Gloss.Char k) Gloss.Up _ _ ->
      if isDigit k || elem k ['.', '-', '+']
        then [] 
        else case keyControl (Options.controlKeys options) k of
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
          dragControls = flip fmap
            (Set.toList $ grabbedControls world)
            (\i ->
              let dy = y - snd (getWorldMousePos identity world)
                  h = slidersHeight Vector.! i
                  dragAmount = float2Double dy / h
              in DragControl i dragAmount)
          setMousePos = SetMousePos (x, y)
      in setMousePos : dragControls
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
      let newCtrl = case (getWorldControls . getControlAt i) identity world of
            Just (LinearControl n l u v) ->
              (LinearControl n l u (bounded l u $ v + (u - l) * amount))
            Just (QuadraticControl n r v) ->
              let x = quadraticToLinear r v
                  v' = quadraticFromLinear r (x + amount * 2 * r)
              in QuadraticControl n r v'
            Nothing -> panic $ "updateEvent: No control at index " <> show i
      in (setWorldControls . setControlAt i) (const newCtrl) world
    ReadDigit d -> setWorldNumberBuffer (d :) world
    SetControl i n ->
      (setWorldControls . setControlAt i . setControlValue) (const n) world

updateInputs :: Options -> Gloss.Event -> World -> World
updateInputs options event world =
  -- traceShow event $
  (foldl' (.) identity
  $ fmap updateEvent (events options event world))
  $ world


updateTime :: Options -> Float -> World -> World
updateTime options dt world =
  let tResolution = Options.timeResolution options
      traceDuration = Options.traceDuration options
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
      input = controlsToInput ctrls
      ctrls = getWorldControls identity world
  in 
     world
     { worldTime = t
     , worldTrajectories = trajectories
     , worldAngles = newAngles
     , worldLayout = GUI.layout
         (Options.hideControls options)
        (getWorldWindow identity world)
        (controlSpecs <$> controlsList ctrls)
        ((fmap . fmap) (bimap double2Float double2Float) trajectories)
     }
