{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TupleSections #-}
{-# LANGUAGE NoImplicitPrelude #-}

module Electra2Shadow.GUI
  ( Layout(..)
  , LayoutAnswer(..)
  , layoutQuery
  , layoutWithControl
  , viewLayout
  -- DEBUG
  , boxColumn
  , boxRow
  , box
  , extent
  , Dim(..)
  , Justification(..)
  , Length(..)
  , Alignment(..)
  ) where

import Protolude

import qualified Data.Text as Text
import GHC.Float (double2Float, float2Double)
import qualified Graphics.Gloss as Gloss
import qualified Graphics.Gloss.Interface.IO.Interact as Gloss

data Layout = Layout Dim Box LayoutElement
  deriving (Show, Eq)

data LayoutElement =
  Row Justification [(Layout, Alignment)]
  | Column Justification [(Layout, Alignment)]
  | Canvas [[(Gloss.Point)]]
  | Slider { sliderName :: Text
           , sliderLowerBound :: Float
           , sliderUpperBound :: Float
           , sliderValue :: Float
           }
  deriving (Show, Eq)

data Dim =
  Fixed (Length, Length)
  | Stretch
  | StretchVertically {dimWidth :: Length}
  | StretchHorizontally {dimHeight :: Length}
  | StretchRatio {aspectRatio :: Rational} -- width / height
  deriving (Show, Eq)

data Length = AbsoluteLength Float | RelativeLength Float
  deriving (Show, Eq)

absLength :: Float -> Length -> Float
absLength parentLength length = case length of
  (AbsoluteLength l) -> l
  (RelativeLength l) -> l * parentLength

-- Box: (x_min y_min x_max y_max)
type Box = (Float, Float, Float, Float)

row
  :: Justification -> [(Dim -> Box -> Layout, Dim, Alignment)]
  -> Dim -> Box -> Layout
row justification seq dim parent =
  trace ("ROW PARENT " <> show parent :: Text)
  trace ("ROW DIM " <> show dim :: Text)
  trace ("ROW BOX' " <> show box' :: Text)
  trace ("ROW BOXES " <> show boxes :: Text)
  trace ("ROW DAS " <> show das :: Text)
  Layout dim box' row'
  where box' = box dim parent
        row' = Row justification seq'
        seq' :: [(Layout, Alignment)]
        seq' = (\(b, (l, d, a)) -> (l (Stretch) b, a)) <$> zip boxes seq
        boxes = boxRow justification box' das
        das = ((\(_, d, a) -> (d, a)) <$> seq)

column
  :: Justification -> [(Dim -> Box -> Layout, Dim, Alignment)]
  -> Dim -> Box -> Layout
column justification seq dim parent =
  Layout dim box' column'
  where box' = box dim parent
        column' = Column justification seq'
        seq' :: [(Layout, Alignment)]
        seq' = (\(b, (l, d, a)) -> (l Stretch b, a)) <$> zip boxes seq
        boxes = boxColumn justification box' das
        das = ((\(_, d, a) -> (d, a)) <$> seq)

canvas :: [[(Gloss.Point)]] -> Dim -> Box -> Layout
canvas trajectories dim parent =
  Layout dim (box dim parent) (Canvas trajectories)

slider
  :: Text -> Float -> Float -> Float
  -> Dim -> Box -> Layout
slider name lowerBound upperBound value dim parent =
  Layout dim (box dim parent) (Slider name lowerBound upperBound value)

layoutWithControl
  :: (Int, Int)
  -> [(Text, Double, Double, Double)]
  -> [[Gloss.Point]]
  -> Layout
layoutWithControl (width, height) slidersSpecs trajectories =
  trace ("WINDOW " <> show (width, height) :: Text) $
  layout
  where layout = column SpaceAround
          [ ( canvas trajectories
            , StretchRatio 1
            , Center
            )
          , ( row SpaceAround
              (fmap (, StretchRatio (1 % 5), Center) sliders)
            , Fixed (RelativeLength 1, RelativeLength (1.0/6.0))
            , Center
            )
          ]
          (StretchRatio (6 % 10))
          ( - fromIntegral width / 2, - fromIntegral height / 2
          , fromIntegral width / 2, fromIntegral height / 2)
        sliders :: [Dim -> Box -> Layout]
        sliders = fmap (\(a,b,c,d) -> slider a (double2Float b)
                         (double2Float c) (double2Float d))
                       slidersSpecs

-- Compute the coordinates of the box specified by the given dimensions, inside
-- the parent box.
box :: Dim -> Box -> Box
box dim (x1, y1, x2, y2) =
  let (width, height) = extent dim (x2 - x1, y2 - y1)
      xc = (x1 + x2) / 2
      yc = (y1 + y2) / 2
  in (xc - width / 2, yc - height / 2, xc + width / 2, yc + height / 2)

-- Compute the width and height of a box, given the dimension specification and the parent width and height.
extent :: Dim -> (Float, Float) -> (Float, Float)
extent dim (parentWidth, parentHeight) = case dim of
  Fixed (width, height) -> ( absLength parentWidth width
                           , absLength parentHeight height
                           )
  Stretch -> (parentWidth, parentHeight)
  StretchVertically width -> ( absLength parentWidth width
                             , parentHeight)
  StretchHorizontally height -> ( parentWidth
                                , absLength parentHeight height)
  StretchRatio ratio ->
    if parentWidth < parentHeight * fromRational ratio
         then (parentWidth, parentWidth / fromRational ratio)
         else (parentHeight * fromRational ratio, parentHeight)

-- r = w / h
-- w = r * h
-- h = w / r

boxColumn :: Justification -> Box -> [(Dim, Alignment)] -> [Box]
boxColumn justification parent dims =
  let (x1p, y1p, x2p, y2p) = parent
      -- Transform to work in the horizontal direction
      symmetry (x1, y1, x2, y2) = (y1, x1, y2, x2)
      parent' = symmetry parent
      dim' dim = case dim of
        Fixed (width, height) -> Fixed (height, width)
        Stretch -> Stretch
        StretchVertically dimWidth -> StretchHorizontally dimWidth
        StretchHorizontally dimHeight -> StretchVertically dimHeight
        StretchRatio ratio ->
          StretchRatio (denominator ratio % numerator ratio)
      dims' = reverse $ first dim' <$> dims
      boxes' = boxRow justification parent' dims'
      -- rotateInv (x1', y1', x2', y2') = (y1', -x1', y2', -x2')
      boxes = reverse $ symmetry <$> boxes'
  in boxes

boxRow :: Justification -> Box -> [(Dim, Alignment)] -> [Box]
boxRow justification parent dims =
  let (x1p, y1p, x2p, y2p) = parent
      countItems = length dims
      widthParent = x2p - x1p
      heightParent = y2p - y1p
      extentPending :: Dim -> Either ((Float, Float) -> (Float, Float)) (Float, Float)
      extentPending dim = case dim of
        Fixed (width, height) -> Right $ extent dim (widthParent, heightParent)
        Stretch -> Left (extent dim)
        StretchVertically width -> Left (extent dim)
        StretchHorizontally height -> Left (extent dim)
        StretchRatio ratio -> Left (extent dim)
      extentsPending
        :: [Either  ((Float, Float) -> (Float, Float)) (Float, Float)]
      extentsPending = extentPending . fst <$> dims
      sumKnownWidths = sum $ fmap fst $ rights $ extentsPending
      spaceLeftForStretch = max 0.0 (widthParent - sumKnownWidths)
      countStretch = length $ filter isLeft extentsPending
      extents :: [(Float, Float)]
      extents = either (\f -> f (spaceLeftForStretch / fromIntegral countStretch, heightParent)) identity
                <$> extentsPending
      sumWidths = sum $ fst <$> extents
      spaceLeftOver = max 0.0 (widthParent - sumWidths)
      spaceBetween :: Float
      spaceBetween = case justification of
        Start -> 0
        End -> 0
        Middle -> 0
        SpaceBetween -> spaceLeftOver / (fromIntegral countItems - 1)
        SpaceAround -> spaceLeftOver / fromIntegral countItems
        SpaceEvenly -> spaceLeftOver / (fromIntegral countItems + 1)
      spaceAtEnd :: Float
      spaceAtEnd = case justification of
        Start -> spaceLeftOver
        End -> 0
        Middle -> spaceLeftOver / 2.0
        SpaceBetween -> 0
        SpaceAround -> spaceBetween / 2.0
        SpaceEvenly -> spaceBetween
      verticalPos (height, alignment) =
        case alignment of
          Center -> ((heightParent - height) / 2.0 + y1p,
                     y2p - (heightParent - height) / 2.0)
          Bottom -> (y1p, y1p + height)
          Top -> (y2p - height, y2p)
      boxes :: [Box]
      (_, boxes) = foldr
        (\(extent, alignment) (curPos, boxesAcc)->
          let (x1,x2) = (curPos - fst extent, curPos)
              (y1,y2) = verticalPos (snd extent, alignment)
              newPos = x1 - spaceBetween
          in (newPos, ((x1, y1, x2, y2)) : boxesAcc))
        (x2p - spaceAtEnd, [])
        (zip extents $ snd <$> dims)
  in
     boxes

data Alignment = Center | Bottom | Top
  deriving (Show, Eq)

data Justification = Start | End | Middle | SpaceBetween | SpaceAround | SpaceEvenly
  deriving (Show, Eq)

layoutQuery :: Layout -> (Float, Float) -> LayoutAnswer
layoutQuery layout (x, y) = panic "UNDEFINED !"

data LayoutAnswer =
  SliderAnswer Int Double
  | NoAnswer
  deriving (Show, Eq)

viewLayout :: Layout -> Gloss.Picture
viewLayout (Layout dim box elem) = case elem of
  (Row justification layouts) ->
    Gloss.pictures $ fmap viewLayout $ fmap fst layouts
  (Column justification layouts) ->
    Gloss.pictures $ fmap viewLayout $ fmap fst layouts
  (Canvas trajectories) -> viewTrajectories trajectories box
  (Slider name lower upper value) ->
    viewSlider name lower upper value box

viewTrajectories :: [[(Gloss.Point)]] -> Box -> Gloss.Picture
viewTrajectories trajectories parent =
  assignBox parent
  $ Gloss.color Gloss.yellow
  $ Gloss.pictures
  $ fmap Gloss.line trajectories

viewSlider :: Text -> Float -> Float -> Float -> Box -> Gloss.Picture
viewSlider label lowerBound upperBound value parent =
  Gloss.pictures $ zipWith assignBox boxes pics
  where pics = [upperBoundP, sliderP, lowerBoundP, labelP]
        boxes = boxColumn SpaceBetween parent
          [ (Fixed (RelativeLength 1, RelativeLength 0.1), Center)
          , (Fixed (RelativeLength 1, RelativeLength 0.6), Center)
          , (Fixed (RelativeLength 1, RelativeLength 0.1), Center)
          , (Fixed (RelativeLength 1, RelativeLength 0.1), Center)
          ]
        position = (value - lowerBound) / (upperBound - lowerBound)
        line = Gloss.color Gloss.green $ Gloss.line [(0, -0.5), (0,0.5)]
        mark = Gloss.color Gloss.green
               $ Gloss.translate 0 (position - 0.5)
               $ Gloss.scale 0.1 0.1
               $ Gloss.circleSolid 1
        sliderP = Gloss.pictures [line, mark]
        labelP = Gloss.color Gloss.green $ Gloss.rectangleWire 1 1
          -- Gloss.color Gloss.green $ Gloss.text (Text.unpack label)
        lowerBoundP = Gloss.color Gloss.green $ Gloss.rectangleWire 1 1
          -- Gloss.color Gloss.green $ Gloss.text (show lowerBound)
        upperBoundP = Gloss.color Gloss.green $ Gloss.rectangleWire 1 1
          -- Gloss.color Gloss.green $ Gloss.text (show upperBound)

assignBox :: Box -> Gloss.Picture -> Gloss.Picture
assignBox (x1, y1, x2, y2) =
  Gloss.translate moveX moveY . Gloss.scale scaleX scaleY
  where moveX = ((x1 + x2 - 1) / 2.0) 
        moveY = ((y1 + y2 - 1) / 2.0) 
        scaleX = (x2 - x1) 
        scaleY = (y2 - y1)

data Direction = Horizontal | Vertical
  deriving (Show, Eq)

viewSeq :: [(Box, Gloss.Picture)] -> Gloss.Picture
viewSeq seq = Gloss.pictures $ fmap (uncurry assignBox) seq

