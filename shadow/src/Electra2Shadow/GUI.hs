{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TupleSections #-}
{-# LANGUAGE NoImplicitPrelude #-}

module Electra2Shadow.GUI
  ( Layout(..)
  , LayoutAnswer(..)
  , layoutQuery
  , layout
  , slidersHeight
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

import GHC.Float (double2Float, float2Double)
import qualified Graphics.Gloss as Gloss

data Layout = Layout Dim Box LayoutElement
  deriving (Show, Eq)

data LayoutElement =
  Row Justification [(Layout, Alignment)]
  | Column Justification [(Layout, Alignment)]
  | Canvas [(Double, [(Gloss.Point)])] -- [(light intensity, trajectory)]
  | Slider { sliderIndex :: Int
           , sliderName :: Text
           , sliderLowerBound :: Float
           , sliderUpperBound :: Float
           , sliderValue :: Float
           , sliderNameBox :: Box
           , sliderLowerBox :: Box
           , sliderUpperBox :: Box
           , sliderSliderBox :: Box
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
absLength parentLength length' = case length' of
  (AbsoluteLength l) -> l
  (RelativeLength l) -> l * parentLength

-- Box: (x_min y_min x_max y_max)
type Box = (Float, Float, Float, Float)

root :: (Int, Int) -> Dim -> (Dim -> Box -> Layout) -> Layout
root (width, height) dim layout' = layout' dim box'
  where box' = box dim (-width' / 2, - height' / 2, width' / 2, height' / 2)
        width' = fromIntegral width
        height' = fromIntegral height
row
  :: Justification -> [(Dim -> Box -> Layout, Dim, Alignment)]
  -> Dim -> Box -> Layout
row justification elems dim box' =
  Layout dim box' row'
  where row' = Row justification seq'
        seq' :: [(Layout, Alignment)]
        seq' = (\(b, (l, d, a)) -> (l d b, a)) <$> zip boxes elems
        boxes = boxRow justification box' das
        das = ((\(_, d, a) -> (d, a)) <$> elems)

column
  :: Justification -> [(Dim -> Box -> Layout, Dim, Alignment)]
  -> Dim -> Box -> Layout
column justification elems dim box' =
  Layout dim box' column'
  where column' = Column justification seq'
        seq' :: [(Layout, Alignment)]
        seq' = (\(b, (l, d, a)) -> (l d b, a)) <$> zip boxes elems
        boxes = boxColumn justification box' das
        das = ((\(_, d, a) -> (d, a)) <$> elems)

canvas :: [(Double, [(Gloss.Point)])] -> Dim -> Box -> Layout
canvas trajectories dim box' =
  Layout dim box' (Canvas trajectories)

slider
  :: Int -> Text -> Float -> Float -> Float
  -> Dim -> Box -> Layout
slider index name lowerBound upperBound value dim box' =
  let boxes = boxColumn SpaceBetween box'
        [ (Fixed (RelativeLength 1, RelativeLength 0.1), Center)
        , (Fixed (RelativeLength 1, RelativeLength 0.6), Center)
        , (Fixed (RelativeLength 1, RelativeLength 0.1), Center)
        , (Fixed (RelativeLength 1, RelativeLength 0.1), Center)
        ]
      (upperBox, sliderBox, lowerBox, labelBox) = case boxes of
        [u, s, l, n] -> (u, s, l, n)
        _ -> panic "GUI.slider pattern matching error."
  in Layout dim box' (Slider
    { sliderIndex = index
    , sliderName = name
    , sliderLowerBound = lowerBound
    , sliderUpperBound = upperBound
    , sliderValue = value
    , sliderNameBox = labelBox
    , sliderLowerBox = lowerBox
    , sliderUpperBox = upperBox
    , sliderSliderBox = sliderBox
    })

layout
  :: Bool
  -> (Int, Int)
  -> [(Text, Double, Double, Double)]
  -> [(Double, [Gloss.Point])]
  -> Layout
layout showControls (width, height) slidersSpecs trajectories =
  root (width, height) (StretchRatio (6 % 10))
  $ column SpaceAround
          [ ( canvas trajectories
            , StretchRatio 1
            , Center
            )
          , ( row SpaceAround $ fmap (, StretchRatio (1 % 5), Center) sliders
            , Fixed (RelativeLength 1, RelativeLength (1.0/6.0))
            , Center
            )
          ]
  where sliders :: [Dim -> Box -> Layout]
        sliders = fmap (\(i,(a,b,c,d)) -> slider i a (double2Float b)
                         (double2Float c) (double2Float d))
                       (zip [0..] slidersSpecs)

layoutQuery :: Layout -> (Float, Float) -> LayoutAnswer
layoutQuery (Layout _ _ e) (x, y) = case e of
  Row _ las ->
    let target = head $ catMaybes $ foldr
          (\(l, _) -> (((l,) <$> inLayout (x, y) l) :)) [] las
    in case target  of
      Just (l,_) -> layoutQuery l (x, y)
      Nothing -> NoAnswer
  Column _ las ->
    let target = head $ catMaybes $ foldr
          (\(l, _) -> (((l,) <$> inLayout (x, y) l) :)) [] las
    in case target  of
      Just (l,_) -> layoutQuery l (x, y)
      Nothing -> NoAnswer
  Canvas _ -> NoAnswer
  Slider i _ _ _ _ _ _ _ sb -> case inBox (x,y) sb of
    Nothing -> NoAnswer
    Just (_, ry) -> SliderAnswer i ry


inLayout :: (Float, Float) -> Layout -> Maybe (Double, Double)
inLayout pos (Layout _ b _) = inBox pos b

inBox :: (Float, Float) -> Box -> Maybe (Double, Double)
inBox (x, y) (x1, y1, x2, y2) =
  if (x >= x1 && x <= x2 && y >= y1 && y <= y2)
    then Just ( float2Double $ (x - x1) / (x2 + x1)
              , float2Double $ (y - y1) / (y2 - y1))
    else Nothing

data LayoutAnswer =
  SliderAnswer Int Double
  | NoAnswer
  deriving (Show, Eq)

slidersHeight :: Layout -> [Float]
slidersHeight (Layout _ _ e) = case e of
  Row _ las -> las >>= slidersHeight . fst
  Column _ las -> las >>= slidersHeight . fst
  Slider{sliderSliderBox = (_, y1, _, y2) } -> [y2 - y1]
  _ -> []

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
  let 
      -- Transform to work in the horizontal direction
      symmetry (x1, y1, x2, y2) = (y1, x1, y2, x2)
      parent' = symmetry parent
      dim' dim = case dim of
        Fixed (width, height) -> Fixed (height, width)
        Stretch -> Stretch
        StretchVertically dimWidth' -> StretchHorizontally dimWidth'
        StretchHorizontally dimHeight' -> StretchVertically dimHeight'
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
        Fixed (_, _) -> Right $ extent dim (widthParent, heightParent)
        Stretch -> Left (extent dim)
        StretchVertically _ -> Left (extent dim)
        StretchHorizontally _ -> Left (extent dim)
        StretchRatio _ -> Left (extent dim)
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
        (\(extent', alignment) (curPos, boxesAcc)->
          let (x1,x2) = (curPos - fst extent', curPos)
              (y1,y2) = verticalPos (snd extent', alignment)
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

viewLayout :: Layout -> Gloss.Picture
viewLayout (Layout _ box' elem') = case elem' of
  (Row _ layouts) ->
    Gloss.pictures $ fmap viewLayout $ fmap fst layouts
  (Column _ layouts) ->
    Gloss.pictures $ fmap viewLayout $ fmap fst layouts
  (Canvas trajectories) -> viewTrajectories trajectories box'
  (Slider _ name lower upper value nameBox lowerBox upperBox sliderBox) ->
    viewSlider name lower upper value nameBox lowerBox upperBox sliderBox

viewTrajectories :: [(Double, [(Gloss.Point)])] -> Box -> Gloss.Picture
viewTrajectories trajectories parent =
  assignBox parent
  $ Gloss.pictures
  $ (\t -> dot t <> trace t)
  <$> trajectories
  where trace :: (Double, [Gloss.Point]) -> Gloss.Picture
        trace trajectory =
              Gloss.color
                (Gloss.withAlpha (double2Float $ fst trajectory) Gloss.yellow)
            $ Gloss.line (snd trajectory)
        dot :: (Double, [Gloss.Point]) -> Gloss.Picture
        dot trajectory =
              fromMaybe mempty
            $ Gloss.color
                (Gloss.withAlpha (double2Float $ fst trajectory) Gloss.yellow)
          <$> (\f -> f (Gloss.circleSolid (1/200)))
          <$> (\(x, y) -> Gloss.translate x y)
          <$> head (snd trajectory)

viewSlider :: Text -> Float -> Float -> Float -> Box -> Box -> Box -> Box -> Gloss.Picture
viewSlider _ lowerBound upperBound value labelBox lowerBox upperBox sliderBox =
  Gloss.pictures $ zipWith assignBox boxes pics
  where pics = [upperBoundP, sliderP, lowerBoundP, labelP]
        boxes = [upperBox, sliderBox, lowerBox, labelBox]
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

-- viewSeq :: [(Box, Gloss.Picture)] -> Gloss.Picture
-- viewSeq elems = Gloss.pictures $ fmap (uncurry assignBox) elems

