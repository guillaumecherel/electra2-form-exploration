{-# LANGUAGE DeriveGeneric     #-}
{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TupleSections #-}
{-# LANGUAGE NoImplicitPrelude #-}

module Electra2Shadow.Control where

import Protolude 

import qualified Dhall
import qualified Data.ByteString.Lazy as BL
import qualified Data.Csv as CSV
import qualified Data.List as List
import qualified Data.KdMap.Static as KDM
import qualified Data.Map.Strict as Map
import           Data.Map.Strict (Map)
import           Data.String (String)
import qualified Data.Vector as Vector
import           Data.Vector (Vector)

import qualified Electra2Shadow.Model as Model

-- Control name lowerBound upperBound value
data Control =
    LinearControl
      { name :: Text
      , lowerBound :: Double
      , upperBound :: Double
      , value :: Double
      }
  | QuadraticControl
      { name :: Text
      , radius :: Double
      , value :: Double
      }
  deriving (Generic, Show, Eq)

instance Dhall.Interpret Control

getControlName :: (Text -> a) -> Control -> a
getControlName f (LinearControl n _ _ _) = f n
getControlName f (QuadraticControl n _ _) = f n

getControlValue :: (Double -> a) -> Control -> a
getControlValue f (LinearControl _ _ _ v) = f v
getControlValue f (QuadraticControl _ _ v) = f v

setControlValue :: (Double -> Double) -> Control -> Control
setControlValue f (LinearControl n l u v) = (LinearControl n l u (f v))
setControlValue f (QuadraticControl n r v) = (QuadraticControl n r (f v))

quadraticToLinear :: Double -> Double -> Double
quadraticToLinear r v =
  let v' = if v >= 0
                  then sqrt (v * r)
                  else - sqrt (-v * r)
  in bounded (- r) r v'

quadraticFromLinear :: Double -> Double -> Double
quadraticFromLinear r v =
  let v' = if v >= 0
                  then v ** 2 / r
                  else - v ** 2 / r
  in bounded (-r) r v'

bounded :: Double -> Double -> Double -> Double
bounded lower upper x = min upper $ max lower x

specs :: Control -> (Text, Double, Double, Double)
specs (LinearControl n l u v) = (n, l, u, v)
specs (QuadraticControl n r v) =
  (n, -r, r, quadraticToLinear r v)

data Controls =
    ControlInput Control Control Control Control Control Control
    -- v1 v2 v3 phi1 phi2 phi3
  | ControlForm (KDM.KdMap Double ControlsValues ModelInputValues)
      (Vector Control)
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

fromInputValues
  :: (Text, Double) -> (Text, Double) -> (Text, Double)
  -> (Text, Double) -> (Text, Double) -> (Text, Double)
  -> Controls
fromInputValues (v1n, v1) (v2n, v2) (v3n, v3)
  (phi1n, phi1) (phi2n, phi2) (phi3n, phi3) =
    ControlInput
      (QuadraticControl v1n 5 v1)
      (QuadraticControl v2n 5 v2)
      (QuadraticControl v3n 5 v3)
      (LinearControl phi1n 0 (2 * pi) phi1)
      (LinearControl phi2n 0 (2 * pi) phi2)
      (LinearControl phi3n 0 (2 * pi) phi3)

fromMap
  :: [(ControlsValues, ModelInputValues)]
  -> Vector Control
  -> Controls
fromMap cim initialControls =
  ControlForm (KDM.build (controlsValuesList) cim) initialControls

mapFromCSV
  :: [Text]
  -> Text -> Text -> Text -> Text -> Text -> Text
  -> BL.ByteString -> [(ControlsValues, ModelInputValues)]
mapFromCSV ctrlNames v1n v2n v3n phi1n phi2n phi3n csv = case CSV.decodeByName csv :: Either String (CSV.Header, Vector (Map Text Double)) of
  Left err -> panic $ "Could not decode csv data.\n" <> show err
  Right (_, result) ->
    let getFromName :: Text -> Map Text Double -> Either Text Double
        getFromName n r = case Map.lookup n r of
          Nothing -> Left n
          Just a -> Right a
        recordToControlValues :: Map Text Double -> Either Text ControlsValues
        recordToControlValues r = FormValues
          <$> Vector.fromList
          <$> traverse (\n -> (n,) <$> getFromName n r) ctrlNames
        recordToInputValues :: Map Text Double -> Either Text ModelInputValues
        recordToInputValues r = ModelInputValues
          <$> getFromName v1n r
          <*> getFromName v2n r
          <*> getFromName v3n r
          <*> getFromName phi1n r
          <*> getFromName phi2n r
          <*> getFromName phi3n r
        recordToMapEntry r =
          (,) <$> recordToControlValues r <*> recordToInputValues r
    in case traverse recordToMapEntry (Vector.toList result) of
        Left n -> panic $ "controlFromMap: the control name " <> n <> " was not found in the csv file."
        Right a -> a

controlsValues :: Controls -> ControlsValues
controlsValues (ControlInput v1 v2 v3 phi1 phi2 phi3) =
  InputValues $ ModelInputValues
    (getControlValue identity v1)
    (getControlValue identity v2)
    (getControlValue identity v3)
    (getControlValue identity phi1)
    (getControlValue identity phi2)
    (getControlValue identity phi3)
controlsValues (ControlForm _ vs) =
  FormValues ((,) <$> getControlName identity <*> getControlValue identity
              <$> vs)

controlsValuesList :: ControlsValues -> [Double]
controlsValuesList (InputValues (ModelInputValues v1 v2 v3 phi1 phi2 phi3)) =
  [v1, v2, v3, phi1, phi2, phi3]
controlsValuesList (FormValues vs) = Vector.toList $ snd <$> vs

toList :: Controls -> [Control]
toList (ControlInput v1 v2 v3 phi1 phi2 phi3) =
  [ v1, v2, v3, phi1, phi2, phi3 ]
toList (ControlForm _ vs) = Vector.toList vs
 
getControlAt :: Int -> (Control -> a) -> Controls -> Maybe a
getControlAt 0 f ctrls@ControlInput{} = getControlsV1 f ctrls
getControlAt 1 f ctrls@ControlInput{} = getControlsV2 f ctrls
getControlAt 2 f ctrls@ControlInput{} = getControlsV3 f ctrls
getControlAt 3 f ctrls@ControlInput{} = getControlsPhi1 f ctrls
getControlAt 4 f ctrls@ControlInput{} = getControlsPhi2 f ctrls
getControlAt 5 f ctrls@ControlInput{} = getControlsPhi3 f ctrls
getControlAt _ _ ControlInput{} = Nothing
getControlAt i f (ControlForm _ vs) =  f <$> vs Vector.!? i

setControlAt :: Int -> (Control -> Control) -> Controls -> Controls
setControlAt 0 f ctrls@ControlInput{} = setControlsV1 f ctrls
setControlAt 1 f ctrls@ControlInput{} = setControlsV2 f ctrls
setControlAt 2 f ctrls@ControlInput{} = setControlsV3 f ctrls
setControlAt 3 f ctrls@ControlInput{} = setControlsPhi1 f ctrls
setControlAt 4 f ctrls@ControlInput{} = setControlsPhi2 f ctrls
setControlAt 5 f ctrls@ControlInput{} = setControlsPhi3 f ctrls
setControlAt i _ ctrls@ControlInput{} = panic
    $ "controlAt: No control at index " <> show i <> " for " <> show ctrls
setControlAt i f (ControlForm cim vs) =
  let update ctrl = ControlForm cim $ vs Vector.// [(i, f ctrl)]
  in maybe (ControlForm cim vs) update $ vs Vector.!? i
  
  

toInput :: Controls -> Model.Input
toInput (ControlInput v1 v2 v3 phi1 phi2 phi3 ) =
  Model.inputDefault
    (getControlValue identity v1)
    (getControlValue identity v2)
    (getControlValue identity v3)
    (getControlValue identity phi1)
    (getControlValue identity phi2)
    (getControlValue identity phi3)
toInput ctrl@(ControlForm kdm _) =
  let (_, ModelInputValues v1 v2 v3 phi1 phi2 phi3) =
        KDM.nearest kdm (controlsValues ctrl)
  in Model.inputDefault v1 v2 v3 phi1 phi2 phi3

data ControlsValues =
    FormValues (Vector (Text, Double))
    -- numberSingular density moran meanCurvature returnTime
  | InputValues ModelInputValues
  deriving (Show, Eq)


data ModelInputValues =
  ModelInputValues Double Double Double Double Double Double
  -- (v1, v2, v3, phi1, phi2, phi3)
  deriving (Show, Eq)

keyControl :: [Char] -> Char -> Maybe Int
keyControl l k = List.elemIndex k l

