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

displayMode = (Gloss.InWindow "Nice Window" (worldWindow initialWorld) (0, 0))

backgroundColor = Gloss.black

fps :: Int
fps = 60

data World = World
  { worldWindow :: (Int, Int)
  , worldV1 :: Double
  , worldV2 :: Double
  , worldV3 :: Double
  , worldA1 :: Double
  , worldA2 :: Double
  , worldA3 :: Double
  } deriving (Eq, Show)

initialWorld :: World
initialWorld = World
  { worldWindow = (200, 200)
  , worldV1 = 0.0
  , worldV2 = 0.0
  , worldV3 = 0.0
  , worldA1 = 0.0
  , worldA2 = 0.0
  , worldA3 = 0.0
  }

viewSlider :: Float -> Float -> Float -> Gloss.Picture
viewSlider lowerBound upperBound value =
  Gloss.pictures [line, mark]
  where position = (value - lowerBound) / (upperBound - lowerBound)
        line = Gloss.color Gloss.green $ Gloss.line [(0, -0.5), (0,0.5)]
        mark = Gloss.color Gloss.green
               $ Gloss.translate 0 (position - 0.5)
               $ Gloss.scale 0.1 0.1
               $ Gloss.circleSolid 1

assignBox :: Float -> Float -> Float -> Float
          -> Gloss.Picture -> Gloss.Picture
assignBox x1 y1 x2 y2 =
  Gloss.translate moveX moveY . Gloss.scale scaleX scaleY
  where moveX = ((x1 + x2 - 1) / 2.0) 
        moveY = ((y1 + y2 - 1) / 2.0) 
        scaleX = (x2 - x1) 
        scaleY = (y2 - y1)

type Size = (Float, Float)
type Box = ((Float, Float), (Float, Float))

data Alignment = Center | Bottom | Top
  deriving (Show, Eq)

data Justification = Start | End | Middle | SpaceBetween | SpaceAround | SpaceEvenly
  deriving (Show, Eq)

data Direction = Horizontal | Vertical
  deriving (Show, Eq)

viewSeq :: Direction -> Justification -> [(Size, Alignment, Gloss.Picture)]
        -> Gloss.Picture
viewSeq direction justification pictures = finalPicture
  where countItems = length pictures
        sizeAlong = case direction of Horizontal -> fst
                                      Vertical -> snd
        sizeOrtho = case direction of Horizontal -> snd
                                      Vertical -> fst
        sumWidths = sum $ fmap (\(s,_,_) -> sizeAlong s) pictures
        spaceLeftOver = max 0.0 (1 - sumWidths)
        spaceBetween :: Float
        spaceBetween = case justification of
          Start -> 0
          End -> 0
          Middle -> 0
          SpaceBetween -> spaceLeftOver / (fromIntegral countItems - 1)
          SpaceAround -> spaceLeftOver / fromIntegral countItems
          SpaceEvenly -> spaceLeftOver / (fromIntegral countItems + 1)
        spaceBefore :: Float
        spaceBefore = case justification of
          Start -> 0
          End -> spaceLeftOver
          Middle -> spaceLeftOver / 2.0
          SpaceBetween -> 0
          SpaceAround -> spaceBetween / 2.0
          SpaceEvenly -> spaceBetween
        orthogonalSpan (size, alignment, picture) =
          case alignment of
            Center -> ((1 - sizeOrtho size) / 2.0,
                       1 - (1 - sizeOrtho size) / 2.0)
            Bottom -> (0, sizeOrtho size)
            Top -> (1 - sizeOrtho size, 1)
        sizesAlong = (\(s,_,_) -> sizeAlong s) <$> pictures
        sizesOrtho = (\(s,_,_) -> sizeOrtho s) <$> pictures
        boxes :: [(Box, Gloss.Picture)]
        (_, boxes) = foldr
          (\(s, alignment, picture) (curPos, boxesAcc) ->
            let (x1,x2) = (curPos, curPos + sizeAlong s)
                (y1,y2) = orthogonalSpan (s, alignment, picture)
                newPos = curPos + sizeAlong s + spaceBetween
            in (newPos, (((x1, y1), (x2, y2)), picture) : boxesAcc))
          (spaceBefore, [])
          pictures
        finalPicture = case direction of
          Horizontal -> boxes
            & fmap (\(((x1, y1), (x2, y2)), picture) ->
              assignBox x1 y1 x2 y2 picture)
            & Gloss.pictures
          Vertical -> boxes
            & (fmap . second) (Gloss.rotate (90))
            & fmap (\(((x1, y1), (x2, y2)), picture) ->
              assignBox x1 y1 x2 y2 picture)
            & Gloss.pictures
            & Gloss.rotate (-90)

viewCommands :: World -> Gloss.Picture
viewCommands w = Gloss.color Gloss.green $ Gloss.rectangleSolid 1 1

viewModel :: World -> Gloss.Picture
viewModel w = Gloss.color Gloss.yellow $ Gloss.rectangleSolid 1 1


scaleToWindow :: World -> Gloss.Picture -> Gloss.Picture
scaleToWindow world = Gloss.scale windowWidth windowHeight
  where (windowWidth, windowHeight) = bimap fromIntegral fromIntegral
          $ worldWindow world

view :: World -> Gloss.Picture
view w = let v = viewSeq Vertical SpaceAround
                    [ ((0.8, 0.7), Center, viewModel w)
                    , ((1, 0.3), Center, viewCommands w)
                    ]
         in Gloss.pictures
             [ scaleToWindow w v
             , Gloss.color Gloss.white
               $ Gloss.rectangleWire
                  (fromIntegral $ fst (worldWindow w) - 1)
                  (fromIntegral $ snd (worldWindow w) - 1)
             ]

updateInputs :: Gloss.Event -> World -> World
updateInputs event world = 
  case event of
    Gloss.EventResize (width, height) ->
      world{ worldWindow = (width, height) }
    _ -> world

updateTime :: Float -> World -> World
updateTime dt w = w

