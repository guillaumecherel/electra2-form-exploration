{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TupleSections #-}
{-# LANGUAGE NoImplicitPrelude #-}

module Electra2Shadow.Model
  ( Parameter (..)
  , Input
  , inputDefault
  , span
  , trajectories
  , trajectoriesToList
  ) where

import Protolude

import qualified Control.Foldl as Fold

type Point = (Double, Double)

data Parameter = Parameter
  { parameterName :: Text
  , parameterLowerBound :: Double
  , parameterUpperBound :: Double
  , parameterValue :: Double
  } deriving (Eq, Show)

parameterSpeed name value = Parameter
  { parameterName = name
  , parameterLowerBound = -5
  , parameterUpperBound = 5
  , parameterValue = value
  }

parameterAngle name value = Parameter
  { parameterName = name
  , parameterLowerBound = 0
  , parameterUpperBound = 2 * pi
  , parameterValue = value
  }

parameterRadius name value = Parameter
  { parameterName = name
  , parameterLowerBound = 0
  , parameterUpperBound = 10
  , parameterValue = value
  }

data Input = Input
  { parametersV1 :: Parameter
  , parametersV2 :: Parameter
  , parametersV3 :: Parameter
  , parametersPhiB :: Parameter
  , parametersPhiD :: Parameter
  , parametersPhiF :: Parameter
  , parametersAH :: Parameter
  , parametersAI :: Parameter
  , parametersRB :: Parameter
  , parametersRC :: Parameter
  , parametersRD :: Parameter
  , parametersRE :: Parameter
  , parametersRF :: Parameter
  , parametersRG :: Parameter
  , parametersRH :: Parameter
  , parametersRI :: Parameter
  }
  deriving (Eq, Show)

inputDefault :: Double -> Double -> Double -> Double -> Double -> Double -> Input
inputDefault v1 v2 v3 phiB phiD phiF = Input
  
  { -- Vitesses
    parametersV1 = parameterSpeed "v1" v1
  , parametersV2 = parameterSpeed "v2" v2
  , parametersV3 = parameterSpeed "v3" v3
    -- Angles
  , parametersPhiB = parameterAngle "phiB" phiB
  , parametersPhiD = parameterAngle "phiD" phiD
  , parametersPhiF = parameterAngle "phiF" phiF
  , parametersAH = parameterAngle "aH" (atan $ 0.05 / 0.5)
  , parametersAI = parameterAngle "aI" (atan $ 0.05 / 0.5)
    -- Rayons
  , parametersRB = parameterRadius "rB" 0.935
  , parametersRC = parameterRadius "rC" 0.935
  , parametersRD = parameterRadius "rD" 0.43
  , parametersRE = parameterRadius "rE" 0.38
  , parametersRF = parameterRadius "rF" 0.417
  , parametersRG = parameterRadius "rG" 0.417
  , parametersRH = parameterRadius "rH" (sqrt (0.5 ** 2 + 0.05 ** 2))
  , parametersRI = parameterRadius "rI" (sqrt (0.5 ** 2 + 0.05 ** 2))
  } 

data Positions = Positions
  { positionB :: Point
  , positionC :: Point
  , positionD :: Point
  , positionE :: Point
  , positionF :: Point
  , positionG :: Point
  , positionH :: Point
  , positionI :: Point
  }

positions :: Input -> Double -> Positions
positions parameters t =
  let v1 = parameterValue $ parametersV1 parameters
      v2 = parameterValue $ parametersV2 parameters
      v3 = parameterValue $ parametersV3 parameters
      rB = parameterValue $ parametersRB parameters
      rC = parameterValue $ parametersRC parameters
      rD = parameterValue $ parametersRD parameters
      rE = parameterValue $ parametersRE parameters
      rF = parameterValue $ parametersRF parameters
      rG = parameterValue $ parametersRG parameters
      rH = parameterValue $ parametersRH parameters
      rI = parameterValue $ parametersRI parameters
      aH = parameterValue $ parametersAH parameters
      aI = parameterValue $ parametersAI parameters
      phiB = parameterValue $ parametersPhiB parameters
      phiD = parameterValue $ parametersPhiD parameters
      phiF = parameterValue $ parametersPhiF parameters
      bX = rB * cos (2 * pi * v1 * t + phiB)
      bY = rB * sin (2 * pi * v1 * t + phiB)
      cX = rC * cos (2 * pi * v1 * t + phiB + pi)
      cY = rC * sin (2 * pi * v1 * t + phiB + pi)
      hX = rH * cos (2 * pi * v1 * t + phiB - aH)
      hY = rH * sin (2 * pi * v1 * t + phiB - aH)
      iX = rI * cos (2 * pi * v1 * t + phiB + pi + aI)
      iY = rI * sin (2 * pi * v1 * t + phiB + pi + aI)
      dX = hX + rD * cos (2 * pi * (v2 + v1) * t + phiD + phiB - aH)
      dY = hY + rD * sin (2 * pi * (v2 + v1) * t + phiD + phiB - aH)
      eX = hX + rE * cos (2 * pi * (v2 + v1) * t + phiD + phiB - aH + pi)
      eY = hY + rE * sin (2 * pi * (v2 + v1) * t + phiD + phiB - aH + pi)
      fX = iX + rF * cos (2 * pi * (v3 + v1) * t + phiF + phiB + aI + pi)
      fY = iY + rF * sin (2 * pi * (v3 + v1) * t + phiF + phiB + aI + pi)
      gX = iX + rG * cos (2 * pi * (v3 + v1) * t + phiF + phiB + aI + pi + pi)
      gY = iY + rG * sin (2 * pi * (v3 + v1) * t + phiF + phiB + aI + pi + pi)
  in Positions
      { positionB = (bX, bY)
      , positionC = (cX, cY)
      , positionD = (dX, dY)
      , positionE = (eX, eY)
      , positionF = (fX, fY)
      , positionG = (gX, gY)
      , positionH = (hX, hY)
      , positionI = (iX, iY)
      }

data Trajectories = Trajectories
  { trajectoryB :: [Point]
  , trajectoryC :: [Point]
  , trajectoryD :: [Point]
  , trajectoryE :: [Point]
  , trajectoryF :: [Point]
  , trajectoryG :: [Point]
  , trajectoryH :: [Point]
  , trajectoryI :: [Point]
  }

instance Semigroup Trajectories where
  a <> b = Trajectories
    { trajectoryB = trajectoryB a <> trajectoryB b
    , trajectoryC = trajectoryC a <> trajectoryC b
    , trajectoryD = trajectoryD a <> trajectoryD b
    , trajectoryE = trajectoryE a <> trajectoryE b
    , trajectoryF = trajectoryF a <> trajectoryF b
    , trajectoryG = trajectoryG a <> trajectoryG b
    , trajectoryH = trajectoryH a <> trajectoryH b
    , trajectoryI = trajectoryI a <> trajectoryI b
    }

instance Monoid Trajectories where
  mempty = Trajectories
    { trajectoryB = mempty
    , trajectoryC = mempty
    , trajectoryD = mempty
    , trajectoryE = mempty
    , trajectoryF = mempty
    , trajectoryG = mempty
    , trajectoryH = mempty
    , trajectoryI = mempty
    }

positionsToTrajectories :: Positions -> Trajectories
positionsToTrajectories positions = Trajectories
    { trajectoryB = [positionB positions]
    , trajectoryC = [positionC positions]
    , trajectoryD = [positionD positions]
    , trajectoryE = [positionE positions]
    , trajectoryF = [positionF positions]
    , trajectoryG = [positionG positions]
    , trajectoryH = [positionH positions]
    , trajectoryI = [positionI positions]
    }

trajectories :: Input -> Double -> Double -> Double -> Trajectories
trajectories parameters t0 dt tmax =
  let timesteps = [t0, t0 + dt .. tmax]
      states = fmap (positions parameters) timesteps
  in Fold.fold (Fold.foldMap positionsToTrajectories identity) states

trajectoriesToList :: Trajectories -> [[Point]]
trajectoriesToList trajectories =
    [ trajectoryB trajectories
    , trajectoryC trajectories
    , trajectoryD trajectories
    , trajectoryE trajectories
    , trajectoryF trajectories
    , trajectoryG trajectories
    , trajectoryH trajectories
    , trajectoryI trajectories
    ]

-- Radius of a circle within which the trajectories fit.
span :: Input -> Double
span parameters = maximum [rB, rC, rH + rD, rH + rE, rI + rF, rI + rG]
  where
    rB = parameterValue $ parametersRB parameters
    rC = parameterValue $ parametersRC parameters
    rD = parameterValue $ parametersRD parameters
    rE = parameterValue $ parametersRE parameters
    rF = parameterValue $ parametersRF parameters
    rG = parameterValue $ parametersRG parameters
    rH = parameterValue $ parametersRH parameters
    rI = parameterValue $ parametersRI parameters
