{-# LANGUAGE DeriveGeneric     #-}
{-# LANGUAGE BangPatterns      #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TupleSections     #-}
{-# LANGUAGE NoImplicitPrelude #-}
{-# LANGUAGE NamedFieldPuns    #-}
{-# LANGUAGE RecordWildCards   #-}

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
import qualified Data.Text as Text
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
      , lowerBound :: Double
      , center :: Double
      , upperBound :: Double
      , value :: Double
      }
  deriving (Generic, Show, Eq)

instance Dhall.Interpret Control

getControlName :: (Text -> a) -> Control -> a
getControlName f (LinearControl n _ _ _) = f n
getControlName f (QuadraticControl n _ _ _ _) = f n

getControlValue :: (Double -> a) -> Control -> a
getControlValue f (LinearControl _ _ _ v) = f v
getControlValue f (QuadraticControl _ _ _ _ v) = f v

setControlValue :: (Double -> Double) -> Control -> Control
setControlValue f (LinearControl n l u v) = (LinearControl n l u (f v))
setControlValue f (QuadraticControl n l c u v) = (QuadraticControl n l c u (f v))

quadraticToLinear :: Double -> Double -> Double -> Double -> Double
quadraticToLinear l c u v 
  | v < l || v > u = bounded l u v
  | u == l && u == c = c
  | l == c = c + sqrt ((v - c) * (u - c))
  | u == c = c - sqrt ((v - c) * (l - c))
  | otherwise = 
     if v >= c  
	  then c + sqrt ((v - c) * (u - c))
	  else c - sqrt ((v - c) * (l - c))

quadraticFromLinear :: Double -> Double -> Double -> Double -> Double
quadraticFromLinear l c u v 
  | v < l || v > u = bounded l u v
  | u == l && u == c = c
  | l == c = 1 / (u - c) * (v - c) ** 2 + c
  | u == c = 1 / (l - c) * (v - c) ** 2 + c 
  | otherwise = 
     if v >= c  
	  then 1 / (u - c) * (v - c) ** 2 + c
	  else 1 / (l - c) * (v - c) ** 2 + c

bounded :: Double -> Double -> Double -> Double
bounded lower upper x = min upper $ max lower x

specs :: Control -> (Text, Double, Double, Double, Double)
specs (LinearControl n l u v) = (n, l, u, v, v)
specs (QuadraticControl n l c u v) =
  (n, l, u, quadraticToLinear l c u v, v)

data Controls =
    ControlInput
      { controlV1 :: Control
      , controlV2 :: Control
      , controlV3 :: Control
      , controlPhi1 :: Control
      , controlPhi2 :: Control
      , controlPhi3 :: Control
      , controlLightB :: Control
      , controlLightC :: Control
      , controlLightD :: Control
      , controlLightE :: Control
      , controlLightF :: Control
      , controlLightG :: Control
      }
  | ControlForm (KDM.KdMap Double ControlsValues ModelInputValues)
      (Vector Control)
  deriving (Show)


getControlsV1 :: (Control -> a) -> Controls -> Maybe a
getControlsV1 f ControlInput{controlV1, ..} = Just $ f controlV1
getControlsV1 _ _ = Nothing

getControlsV2 :: (Control -> a) ->  Controls -> Maybe a
getControlsV2 f ControlInput{controlV2, ..} = Just $ f controlV2
getControlsV2 _ _ = Nothing

getControlsV3 :: (Control -> a) -> Controls -> Maybe a
getControlsV3 f ControlInput{controlV3, ..} = Just $ f controlV3
getControlsV3 _ _ = Nothing

getControlsPhi1 :: (Control -> a) -> Controls -> Maybe a
getControlsPhi1 f ControlInput{controlPhi1, ..} = Just $ f controlPhi1
getControlsPhi1 _ _ = Nothing

getControlsPhi2 :: (Control -> a) -> Controls -> Maybe a
getControlsPhi2 f ControlInput{controlPhi2, ..} = Just $ f controlPhi2
getControlsPhi2 _ _ = Nothing

getControlsPhi3 :: (Control -> a) -> Controls -> Maybe a
getControlsPhi3 f ControlInput{controlPhi3, ..} = Just $ f controlPhi3
getControlsPhi3 _ _ = Nothing

getControlsLightB :: (Control -> a) -> Controls -> Maybe a
getControlsLightB f ControlInput{controlLightB, ..} = Just $ f controlLightB
getControlsLightB _ _ = Nothing

getControlsLightC :: (Control -> a) -> Controls -> Maybe a
getControlsLightC f ControlInput{controlLightC, ..} = Just $ f controlLightC
getControlsLightC _ _ = Nothing

getControlsLightD :: (Control -> a) -> Controls -> Maybe a
getControlsLightD f ControlInput{controlLightD, ..} = Just $ f controlLightD
getControlsLightD _ _ = Nothing

getControlsLightE :: (Control -> a) -> Controls -> Maybe a
getControlsLightE f ControlInput{controlLightE, ..} = Just $ f controlLightE
getControlsLightE _ _ = Nothing

getControlsLightF :: (Control -> a) -> Controls -> Maybe a
getControlsLightF f ControlInput{controlLightF, ..} = Just $ f controlLightF
getControlsLightF _ _ = Nothing

getControlsLightG :: (Control -> a) -> Controls -> Maybe a
getControlsLightG f ControlInput{controlLightG, ..} = Just $ f controlLightG
getControlsLightG _ _ = Nothing

setControlsV1 :: (Control -> Control) -> Controls -> Controls
setControlsV1 f ControlInput{controlV1, ..} =
  ControlInput {controlV1 = (f controlV1), ..}
setControlsV1 _ c = c

setControlsV2 :: (Control -> Control) -> Controls -> Controls
setControlsV2 f ControlInput{controlV2, ..} =
  ControlInput{controlV2 = (f controlV2), ..}
setControlsV2 _ c = c

setControlsV3 :: (Control -> Control) -> Controls -> Controls
setControlsV3 f ControlInput{controlV3, ..} =
  ControlInput{controlV3 = (f controlV3), ..}
setControlsV3 _ c = c

setControlsPhi1 :: (Control -> Control) -> Controls -> Controls
setControlsPhi1 f ControlInput{controlPhi1, ..} =
  ControlInput{controlPhi1 = (f controlPhi1), ..}
setControlsPhi1 _ c = c

setControlsPhi2 :: (Control -> Control) -> Controls -> Controls
setControlsPhi2 f ControlInput{controlPhi2, ..} =
  ControlInput{controlPhi2 = (f controlPhi2), ..}
setControlsPhi2 _ c = c

setControlsPhi3 :: (Control -> Control) -> Controls -> Controls
setControlsPhi3 f ControlInput{controlPhi3, ..} =
  ControlInput{controlPhi3 = (f controlPhi3), ..}
setControlsPhi3 _ c = c

setControlsLightB :: (Control -> Control) -> Controls -> Controls
setControlsLightB f ControlInput{controlLightB, ..} =
  ControlInput{controlLightB = (f controlLightB), ..}
setControlsLightB _ c = c

setControlsLightC :: (Control -> Control) -> Controls -> Controls
setControlsLightC f ControlInput{controlLightC, ..} =
  ControlInput{controlLightC = (f controlLightC), ..}
setControlsLightC _ c = c

setControlsLightD :: (Control -> Control) -> Controls -> Controls
setControlsLightD f ControlInput{controlLightD, ..} =
  ControlInput{controlLightD = (f controlLightD), ..}
setControlsLightD _ c = c

setControlsLightE :: (Control -> Control) -> Controls -> Controls
setControlsLightE f ControlInput{controlLightE, ..} =
  ControlInput{controlLightE = (f controlLightE), ..}
setControlsLightE _ c = c

setControlsLightF :: (Control -> Control) -> Controls -> Controls
setControlsLightF f ControlInput{controlLightF, ..} =
  ControlInput{controlLightF = (f controlLightF), ..}
setControlsLightF _ c = c

setControlsLightG :: (Control -> Control) -> Controls -> Controls
setControlsLightG f ControlInput{controlLightG, ..} =
  ControlInput{controlLightG = (f controlLightG), ..}
setControlsLightG _ c = c



fromInputValues
  :: (Text, Double) -> (Text, Double) -> (Text, Double)
  -> (Text, Double) -> (Text, Double) -> (Text, Double)
  -> (Text, Double) -> (Text, Double) -> (Text, Double)
  -> (Text, Double) -> (Text, Double) -> (Text, Double)
  -> Controls
fromInputValues
  (v1n, v1) (v2n, v2) (v3n, v3)
  (phi1n, phi1) (phi2n, phi2) (phi3n, phi3)
  (lightBn, lightB) (lightCn, lightC) (lightDn, lightD)
  (lightEn, lightE) (lightFn, lightF) (lightGn, lightG) =
    ControlInput
      (QuadraticControl v1n (-5) 0 5 v1)
      (QuadraticControl v2n (-5) 0 5 v2)
      (QuadraticControl v3n (-5) 0 5 v3)
      (LinearControl phi1n 0 (2 * pi) phi1)
      (LinearControl phi2n 0 (2 * pi) phi2)
      (LinearControl phi3n 0 (2 * pi) phi3)
      (LinearControl lightBn 0 1 lightB)
      (LinearControl lightCn 0 1 lightC)
      (LinearControl lightDn 0 1 lightD)
      (LinearControl lightEn 0 1 lightE)
      (LinearControl lightFn 0 1 lightF)
      (LinearControl lightGn 0 1 lightG)

fromMap
  :: [(ControlsValues, ModelInputValues)]
  -> Vector Control
  -> Controls
fromMap cim initialControls =
  ControlForm (KDM.build (controlsValuesList) cim) initialControls

mapFromCSV
  :: [Text]
  -> Text -> Text -> Text -> Text -> Text -> Text
  -> Text -> Text -> Text -> Text -> Text -> Text
  -> BL.ByteString -> [(ControlsValues, ModelInputValues)]
mapFromCSV
  ctrlNames
  v1n v2n v3n phi1n phi2n phi3n
  lightBn lightCn lightDn lightEn lightFn lightGn
  csv = case CSV.decodeByName csv :: Either String (CSV.Header, Vector (Map Text Double)) of
  Left err -> panic $ "Could not decode csv data.\n" <> show err
  Right (_, result) ->
    let getFromName :: Text -> Map Text Double -> Either String Double
        getFromName n r = case Map.lookup n r of
          Nothing -> Left $ Text.unpack n
          Just a -> Right a
        recordToControlValues :: Map Text Double -> Either String ControlsValues
        recordToControlValues r = FormValues
          <$> Vector.fromList
          <$> traverse (\n -> (n,) <$> getFromName n r) ctrlNames
        recordToInputValues :: Map Text Double -> Either String ModelInputValues
        recordToInputValues r = ModelInputValues
          <$> getFromName v1n r
          <*> getFromName v2n r
          <*> getFromName v3n r
          <*> (getFromName phi1n r <|> Right 0)
          <*> (getFromName phi2n r <|> Right 0)
          <*> (getFromName phi3n r <|> Right 0)
          <*> fmap (\x -> if x < 0.5 then 0 else 1) (getFromName lightBn r)
          <*> fmap (\x -> if x < 0.5 then 0 else 1) (getFromName lightCn r)
          <*> fmap (\x -> if x < 0.5 then 0 else 1) (getFromName lightDn r)
          <*> fmap (\x -> if x < 0.5 then 0 else 1) (getFromName lightEn r)
          <*> fmap (\x -> if x < 0.5 then 0 else 1) (getFromName lightFn r)
          <*> fmap (\x -> if x < 0.5 then 0 else 1) (getFromName lightGn r)
        recordToMapEntry r =
          (,) <$> recordToControlValues r <*> recordToInputValues r
    in case traverse recordToMapEntry (Vector.toList result) of
        Left n -> panic $ "controlFromMap: the control name " <> Text.pack n <> " was not found in the csv file."
        Right a -> a

controlsValues :: Controls -> ControlsValues
controlsValues (ControlInput v1 v2 v3 phi1 phi2 phi3 lightB lightC lightD lightE lightF lightG) =
  InputValues $ ModelInputValues
    (getControlValue identity v1)
    (getControlValue identity v2)
    (getControlValue identity v3)
    (getControlValue identity phi1)
    (getControlValue identity phi2)
    (getControlValue identity phi3)
    (getControlValue identity lightB)
    (getControlValue identity lightC)
    (getControlValue identity lightD)
    (getControlValue identity lightE)
    (getControlValue identity lightF)
    (getControlValue identity lightG)
controlsValues (ControlForm _ vs) =
  FormValues ((,) <$> getControlName identity <*> getControlValue identity
              <$> vs)

controlsValuesList :: ControlsValues -> [Double]
controlsValuesList (InputValues (ModelInputValues v1 v2 v3 phi1 phi2 phi3 lB lC lD lE lF lG)) =
  [v1, v2, v3, phi1, phi2, phi3, lB, lC, lD, lE, lF, lG]
controlsValuesList (FormValues vs) = Vector.toList $ snd <$> vs

toList :: Controls -> [Control]
toList (ControlInput v1 v2 v3 phi1 phi2 phi3 lightB lightC lightD lightE lightF lightG) =
  [ v1, v2, v3, phi1, phi2, phi3
  , lightB, lightC, lightD, lightE, lightF, lightG ]
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
  
  

toInput :: Controls -> (ControlsValues, Model.Input)
toInput (ControlInput v1 v2 v3 phi1 phi2 phi3 lB lC lD lE lF lG) =
  ( controlsValues (ControlInput v1 v2 v3 phi1 phi2 phi3 lB lC lD lE lF lG)
  , Model.inputDefault
    (getControlValue identity v1)
    (getControlValue identity v2)
    (getControlValue identity v3)
    (getControlValue identity phi1)
    (getControlValue identity phi2)
    (getControlValue identity phi3)
    (getControlValue identity lB)
    (getControlValue identity lC)
    (getControlValue identity lD)
    (getControlValue identity lE)
    (getControlValue identity lF)
    (getControlValue identity lG))
toInput ctrl@(ControlForm kdm _) =
  let (ctrlvals, ModelInputValues v1 v2 v3 phi1 phi2 phi3 lB lC lD lE lF lG) =
        KDM.nearest kdm (controlsValues ctrl)
  in (ctrlvals, Model.inputDefault v1 v2 v3 phi1 phi2 phi3 lB lC lD lE lF lG)

data ControlsValues =
    FormValues (Vector (Text, Double))
    -- numberSingular density moran meanCurvature returnTime
  | InputValues ModelInputValues
  deriving (Show, Eq)


data ModelInputValues =
  ModelInputValues
    -- (v1, v2, v3, phi1, phi2, phi3)
    Double Double Double Double Double Double
    -- light
    Double Double Double Double Double Double
  deriving (Show, Eq)

keyControl :: [Char] -> Char -> Maybe Int
keyControl l k = List.elemIndex k l

