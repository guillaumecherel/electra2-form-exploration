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
    , fps
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

import Protolude hiding (get)

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

import qualified Electra2Shadow.GUI as GUI
import qualified Electra2Shadow.Model as Model

displayMode :: Gloss.Display
displayMode = (Gloss.InWindow "Nice Window" (400, 600) (0, 0))

backgroundColor :: Gloss.Color
backgroundColor = Gloss.black

fps :: Int
fps = 60

-- Control name lowerBound upperBound value
data Control =
    LinearControl Text Double Double Double 
    -- label lowerBound upperBound value
  | QuadraticControl Text Double Double
    -- label radius value
  deriving (Show, Eq)

getControlValue :: Control -> Double
getControlValue (LinearControl _ _ _ v) = v
getControlValue (QuadraticControl _ _ v) = v

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


getControlsV1 ::  Controls -> Maybe Control
getControlsV1 (ControlInput v1 _ _ _ _ _) = Just v1
getControlsV1 _ = Nothing

getControlsV2 ::  Controls -> Maybe Control
getControlsV2 (ControlInput _ v2 _ _ _ _) = Just v2
getControlsV2 _ = Nothing

getControlsV3 :: Controls -> Maybe Control
getControlsV3 (ControlInput _ _ v3 _ _ _) = Just v3
getControlsV3 _ = Nothing

getControlsPhi1 :: Controls -> Maybe Control
getControlsPhi1 (ControlInput _  _  _ phi1 _  _) = Just phi1
getControlsPhi1 _ = Nothing

getControlsPhi2 :: Controls -> Maybe Control
getControlsPhi2 (ControlInput _  _  _  _ phi2 _) = Just phi2
getControlsPhi2 _ = Nothing

getControlsPhi3 :: Controls -> Maybe Control
getControlsPhi3 (ControlInput _  _  _  _  _ phi3) = Just phi3
getControlsPhi3 _ = Nothing

getControlsNumberSingular :: Controls -> Maybe Control
getControlsNumberSingular (ControlForm _ ns _  _  _  _) = Just ns
getControlsNumberSingular _ = Nothing

getControlsDensity :: Controls -> Maybe Control
getControlsDensity (ControlForm _  _ d _  _  _) = Just d
getControlsDensity _ = Nothing

getControlsMoran :: Controls -> Maybe Control
getControlsMoran (ControlForm _  _  _ m _  _) = Just m
getControlsMoran _ = Nothing

getControlsMeanCurvature :: Controls -> Maybe Control
getControlsMeanCurvature (ControlForm _  _  _  _ mc _) = Just mc
getControlsMeanCurvature _ = Nothing

getControlsReturnTime :: Controls -> Maybe Control
getControlsReturnTime (ControlForm _  _  _  _  _ rt) = Just rt
getControlsReturnTime _ = Nothing

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
    (LinearControl "singular" 0 1000 ns)
    (LinearControl "density" 0 1 d)
    (LinearControl "moran" (-1) 1 m)
    (QuadraticControl "curv" 5 mc)
    (LinearControl "return time" 0 5 rt)

controlMapFromCSV :: BL.ByteString -> ([Text], [(ControlsValues, ModelInputValues)])
controlMapFromCSV csv = case CSV.decodeByName csv of
  Left err -> panic $ "Could not decode csv data.\n" <> show err
  Right (header, result) ->
    let names = decodeUtf8 <$> Vector.toList header
    in (names, Vector.toList result)

controlsValues :: Controls -> ControlsValues
controlsValues (ControlInput v1 v2 v3 phi1 phi2 phi3) =
  InputValues $ ModelInputValues
    (getControlValue v1)
    (getControlValue v2)
    (getControlValue v3)
    (getControlValue phi1)
    (getControlValue phi2)
    (getControlValue phi3)
controlsValues (ControlForm _ ns d m mc rt) =
  FormValues
    (getControlValue ns)
    (getControlValue d)
    (getControlValue m)
    (getControlValue mc)
    (getControlValue rt)

controlsValuesList :: ControlsValues -> [Double]
controlsValuesList (InputValues (ModelInputValues v1 v2 v3 phi1 phi2 phi3)) =
  [v1, v2, v3, phi1, phi2, phi3]
controlsValuesList (FormValues ns d m mc rt) = [ns, d, m, mc, rt]

-- controlsValues :: (ControlsValues -> (a, ControlsValues)) -> Controls -> (a, Controls)
-- controlsValues f (ControlForm kdm ns d m mc rt) =
--   (fst $ f vals, controlsFromMap kdm $ snd $ f vals)
--   where vals = FormValues
--                  (get controlValue ns)
--                  (get controlValue d)
--                  (get controlValue m)
--                  (get controlValue mc)
--                  (get controlValue rt)
-- controlsValues f (ControlInput v1 v2 v3 phi1 phi2 phi3) =
--   (fst $ f vals, controlsFromInputValues $ snd $ f vals)
--   where vals = 
      

controlsList :: Controls -> [Control]
controlsList (ControlInput v1 v2 v3 phi1 phi2 phi3) =
  [ v1, v2, v3, phi1, phi2, phi3 ]
controlsList (ControlForm _ ns d m mc rt) =
  [ ns, d, m, mc, rt ]
 
getControlAt :: Int -> Controls -> Maybe Control
getControlAt 0 ctrls@ControlInput{} = getControlsV1 ctrls
getControlAt 1 ctrls@ControlInput{} = getControlsV2 ctrls
getControlAt 2 ctrls@ControlInput{} = getControlsV3 ctrls
getControlAt 3 ctrls@ControlInput{} = getControlsPhi1 ctrls
getControlAt 4 ctrls@ControlInput{} = getControlsPhi2 ctrls
getControlAt 5 ctrls@ControlInput{} = getControlsPhi3 ctrls
getControlAt 0 ctrls@ControlForm{} = getControlsNumberSingular ctrls
getControlAt 1 ctrls@ControlForm{} = getControlsDensity ctrls
getControlAt 2 ctrls@ControlForm{} = getControlsMoran ctrls
getControlAt 3 ctrls@ControlForm{} = getControlsMeanCurvature ctrls
getControlAt 4 ctrls@ControlForm{} = getControlsReturnTime ctrls
getControlAt _ _ = Nothing

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
    (getControlValue v1)
    (getControlValue v2)
    (getControlValue v3)
    (getControlValue phi1)
    (getControlValue phi2)
    (getControlValue phi3)
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
    <$> m CSV..: "numberPointSinguliersTotal"
    <*> m CSV..: "totalDensity"
    <*> m CSV..: "totalmoran"
    <*> m CSV..: "meanCourbures"
    <*> m CSV..: "nbTempsRetour"

instance CSV.FromNamedRecord ModelInputValues where
  parseNamedRecord m = ModelInputValues
    <$> m CSV..: "v1"
    <*> m CSV..: "v2"
    <*> m CSV..: "v3"
    <*> m CSV..: "angleIni_B"
    <*> m CSV..: "angleIni_D"
    <*> m CSV..: "angleIni_F"

instance CSV.FromNamedRecord (ControlsValues, ModelInputValues) where
  parseNamedRecord m = (,)
    <$> CSV.parseNamedRecord m
    <*> CSV.parseNamedRecord m


keyControl :: Char -> Maybe Int
keyControl k = List.elemIndex k
  ['u', 'i', 'e'
  ,'\195', 'p', 'o']

data World = World
  { getWorldWindow :: (Int, Int)
  , getWorldControls :: Controls
  , getWorldTime :: Double
  , getWorldTrajectories :: [[(Double, Double)]]
  , getWorldAngles :: Model.Angles
  , getWorldLayout :: GUI.Layout
  , getWorldMouseGrabControl :: Maybe Int
  , getWorldKeyboardGrabControl :: Set Int
  , getWorldMousePos :: (Float, Float)
  , getWorldNumberBuffer :: [Char]
  } deriving (Show)

setWorldWindow :: ((Int, Int) -> (Int, Int)) -> World -> World
setWorldWindow f w = w { getWorldWindow = f $ getWorldWindow w }

setWorldControls :: (Controls -> Controls) -> World -> World
setWorldControls f w = w { getWorldControls = f $ getWorldControls w }

-- worldTimeSet :: (Double -> Double) -> World -> World
-- worldTimeSet f w = w { worldTime = f $ worldTime w }
-- worldTimeSet' :: Double -> World -> World
-- worldTimeSet' v = worldTimeSet (const v)

-- worldTrajectoriesSet :: ([[(Double, Double)]] -> [[(Double, Double)]]) -> World -> World
-- worldTrajectoriesSet f w = w { worldTrajectories = f $ worldTrajectories w }
-- worldTrajectoriesSet' :: [[(Double, Double)]] -> World -> World
-- worldTrajectoriesSet' v = worldTrajectoriesSet (const v)

-- worldAnglesSet :: (Model.Angles -> Model.Angles) -> World -> World
-- worldAnglesSet f w = w { worldAngles = f $ worldAngles w }
-- worldAnglesSet' :: Model.Angles -> World -> World
-- worldAnglesSet' v = worldAnglesSet (const v)

-- worldLayoutSet :: (GUI.Layout -> GUI.Layout) -> World -> World
-- worldLayoutSet f w = w { worldLayout = f $ worldLayout w }
-- worldLayoutSet' :: GUI.Layout -> World -> World
-- worldLayoutSet' v = worldLayoutSet (const v)

setWorldMouseGrabControl :: (Maybe Int -> Maybe Int) -> World -> World
setWorldMouseGrabControl f w =
  w { getWorldMouseGrabControl = f $ getWorldMouseGrabControl w }

setWorldKeyboardGrabControl :: (Set Int -> Set Int) -> World -> World
setWorldKeyboardGrabControl f w =
  w { getWorldKeyboardGrabControl = f $ getWorldKeyboardGrabControl w }

setWorldMousePos :: ((Float, Float) -> (Float, Float)) -> World -> World
setWorldMousePos f w = w { getWorldMousePos = f $ getWorldMousePos w }

setWorldNumberBuffer :: ([Char] -> [Char]) -> World -> World
setWorldNumberBuffer f w = w { getWorldNumberBuffer = f $ getWorldNumberBuffer w }

grabbedControls :: World -> Set Int
grabbedControls world =
  getWorldKeyboardGrabControl world
  <> Set.fromList (toList $ getWorldMouseGrabControl world)

initialWorld :: Controls -> World
initialWorld controls =
  let windowSize = case displayMode of
        (Gloss.InWindow _ xy _) -> xy
        Gloss.FullScreen -> (800, 800)
  in World
      { getWorldWindow = windowSize
      , getWorldControls = controls
      , getWorldTime = 0
      , getWorldTrajectories = replicate 8 []
      , getWorldAngles = Model.Angles 0 0 0
      , getWorldLayout = GUI.layoutWithControl
          windowSize
          (controlSpecs <$> controlsList controls)
          []
      , getWorldMouseGrabControl = Nothing
      , getWorldKeyboardGrabControl = mempty
      , getWorldMousePos = (0, 0)
      , getWorldNumberBuffer = []
      }

view :: World -> Gloss.Picture
view world =
  Gloss.pictures
    [ GUI.viewLayout (getWorldLayout world)
    , Gloss.color Gloss.white
      $ Gloss.rectangleWire
          (fromIntegral $ fst (getWorldWindow world) - 1)
          (fromIntegral $ snd (getWorldWindow world) - 1)
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

events :: Gloss.Event -> World -> [Event]
events event world =
  -- trace (show event :: Text) $
  case event of
    Gloss.EventResize (width, height) -> [ResizeWindow (width, height)]
    Gloss.EventKey (Gloss.MouseButton Gloss.LeftButton) Gloss.Up _ _ ->
      case getWorldMouseGrabControl world of
        Nothing -> []
        Just i -> [MouseReleaseControl i]
    Gloss.EventKey (Gloss.MouseButton Gloss.LeftButton) Gloss.Down _ (x, y) ->
      case GUI.layoutQuery (getWorldLayout world) (x, y) of
        GUI.SliderAnswer i _ -> [MouseGrabControl i]
        GUI.NoAnswer -> []
    Gloss.EventKey (Gloss.Char k) Gloss.Down _ _ ->
      if isDigit k || elem k ['.', '-', '+']
        then [ReadDigit k]
        else 
          toList $ KeyboardGrabControl <$> keyControl k
    Gloss.EventKey (Gloss.Char k) Gloss.Up _ _ ->
      if isDigit k || elem k ['.', '-', '+']
        then [] 
        else case keyControl k of
          Nothing -> []
          Just i ->
            let value = double (Text.pack $ reverse $ getWorldNumberBuffer world)
                setctrls = case value of
                  Right (n, _) -> [SetControl i n]
                  Left err -> [] &
                    trace ("getWorldNumberBuffer cannot be read as a number: "
                      <> show (getWorldNumberBuffer world)
                      <> " " <> Text.pack err :: Text)
                releases = [KeyboardReleaseControl i]
            in setctrls <> releases
    Gloss.EventMotion (x, y) ->
      let slidersHeight = fmap float2Double $ Vector.fromList
            $ GUI.slidersHeight (getWorldLayout world)
          dragControls = flip fmap
            (Set.toList $ grabbedControls world)
            (\i ->
              let dy = y - snd (getWorldMousePos world)
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
      let keyGrab = Set.delete i $ getWorldKeyboardGrabControl world
          numBuf = if null keyGrab
                      then mempty
                      else (getWorldNumberBuffer world)
      in setWorldKeyboardGrabControl (const keyGrab)
         $ setWorldNumberBuffer (const numBuf)
         $ world
    SetMousePos pos -> setWorldMousePos (const pos) world
    DragControl i amount ->
      let newCtrl = case getControlAt i (getWorldControls world) of
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

updateInputs :: Gloss.Event -> World -> World
updateInputs event world = 
  -- traceShow event $
  (foldl' (.) identity
  $ fmap updateEvent (events event world))
  $ world


updateTime :: Float -> World -> World
updateTime dt world =
  let newAngles = Model.angles
        (Model.inputSetAngles input (getWorldAngles world))
        (float2Double dt)
      newTrajectories =
        (fmap . fmap) (bimap scale scale)
        $ Model.trajectoriesToList
        $ Model.trajectories inputShiftAngles 0 tResolution (float2Double dt)
      inputShiftAngles = Model.inputAddAngles input newAngles
      tResolution = 0.005
      t = getWorldTime world + float2Double dt
      scale x = x / (2 * Model.span input)
      trajectories = zipWith
        (\new prev ->
          take (ceiling $ traceDuration / tResolution)
            (new <> prev))
        newTrajectories (getWorldTrajectories world)
      traceDuration = 1.0 / 25.0
      input = controlsToInput ctrls
      ctrls = getWorldControls world
  in 
     world
     { getWorldTime = t
     , getWorldTrajectories = trajectories
     , getWorldAngles = newAngles
     , getWorldLayout = GUI.layoutWithControl
        (getWorldWindow world)
        (controlSpecs <$> controlsList ctrls)
        ((fmap . fmap) (bimap double2Float double2Float) trajectories)
     }
