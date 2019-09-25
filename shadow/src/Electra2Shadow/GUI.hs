{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TupleSections #-}
{-# LANGUAGE NoImplicitPrelude #-}

module Electra2Shadow.GUI
  ( Layout(..)
  , LayoutAnswer(..)
  , layoutQuery
  , layout
  , slidersSpan
  , SliderSpan(..)
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

data Layout = Layout Dim Box LayoutElement
  deriving (Show, Eq)

data LayoutElement =
  Row Justification [(Layout, Alignment)]
  | Column Justification [(Layout, Alignment)]
  | Canvas [(Double, [(Gloss.Point)])] -- [(light intensity, trajectory)]
  | Slider { sliderVertical :: Bool
           , sliderIndex :: Int
           , sliderName :: Text
           , sliderLowerBound :: Float
           , sliderUpperBound :: Float
           , sliderValue :: Float
           , sliderSliderBox :: Box
           }
  | TextBox { scaleX :: Double, scaleY :: Double, content :: Text }
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
  :: Bool -> Int -> Text -> Float -> Float -> Float
  -> Dim -> Box -> Layout
slider vertical index name lowerBound upperBound value dim box' =
  Layout dim box' (Slider
    { sliderVertical = vertical
    , sliderIndex = index
    , sliderName = name
    , sliderLowerBound = lowerBound
    , sliderUpperBound = upperBound
    , sliderValue = value
    , sliderSliderBox = box'
    })

textBox :: Double -> Double -> Text -> Dim -> Box -> Layout
textBox scaleX' scaleY' txt dim box' = Layout dim box' (TextBox scaleX' scaleY' txt)

layout
  :: Bool
  -> (Int, Int)
  -> [((Text, Double, Double, Double, Double), Double)]
  -> [(Text, Double)]
  -> [(Double, [Gloss.Point])]
  -> Layout
layout hideControls (width, height) slidersSpecs speeds trajectories =
  root (width, height) (Stretch)
  $ row Middle
    [ ( canvas trajectories
      , StretchRatio 1
      , Center
      )
    , ( column Start
        (    mconcat ( fmap
               (\(index, ((txt, lower, upper, sliderVal, val), nearest)) ->
                 ( textBox
                     (1 / 10) (1 / 10)
                     (txt <> ": "
                      -- <> show (fromIntegral
                      --   (round (100 * sliderVal) :: Int) / 100 :: Double)
                      -- <> " // "
                      <> show (fromIntegral
                         (round (100 * val) :: Int) / 100 :: Double)
                      <> " (" <> show (fromIntegral
                         (round (100 * nearest) :: Int) / 100 :: Double)
                      <> ")")
                   , StretchHorizontally (AbsoluteLength 30)
                   , Bottom
                   )
                 : if hideControls then []
                     else [( slider False index txt (double2Float lower)
                               (double2Float upper) (double2Float sliderVal)
                           , StretchHorizontally (AbsoluteLength 20)
                           , Bottom )]
               )
               $ zip [0..] slidersSpecs)
          <> [( textBox (1 / 10) (1 / 10) ""
              , StretchHorizontally (AbsoluteLength 30), Bottom )]
          <> fmap
               (\(txt, val) ->
                 ( textBox
                     (1 / 10) (1 / 10)
                     (txt <> ": "<> show (fromIntegral (round (100 * val) :: Int) / 100 :: Double))
                 , StretchHorizontally (AbsoluteLength 30)
                 , Bottom
                 ))
               speeds
          <> [ ( textBox (1 / 10) (1 / 10)
                   (mconcat $ bool "0" "1" . (>= 0.5) . fst <$> trajectories)
               , StretchHorizontally (AbsoluteLength 30)
               , Bottom
               )
             ]
          )
      , StretchVertically (RelativeLength 0.3)
      , Center
      )
    ]

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
  Slider v i _ _ _ _ sb -> case inBox (x,y) sb of
    Nothing -> NoAnswer
    Just (rx, ry) -> SliderAnswer i (if v then ry else rx)
  TextBox _ _ _ -> NoAnswer


inLayout :: (Float, Float) -> Layout -> Maybe (Double, Double)
inLayout pos (Layout _ b _) = inBox pos b

inBox :: (Float, Float) -> Box -> Maybe (Double, Double)
inBox (x, y) (x1, y1, x2, y2) =
  if (x >= x1 && x <= x2 && y >= y1 && y <= y2)
    then Just ( float2Double $ (x - x1) / (x2 - x1)
              , float2Double $ (y - y1) / (y2 - y1))
    else Nothing

data LayoutAnswer =
  SliderAnswer Int Double
  | NoAnswer
  deriving (Show, Eq)

data SliderSpan =
    HorizontalSpan Float
  | VerticalSpan Float

slidersSpan :: Layout -> [SliderSpan]
slidersSpan (Layout _ _ e) = case e of
  Row _ las -> las >>= slidersSpan . fst
  Column _ las -> las >>= slidersSpan . fst
  Slider{sliderVertical = vertical, sliderSliderBox = (x1, y1, x2, y2) } ->
    if vertical then [VerticalSpan $ y2 - y1] else [HorizontalSpan $ x2 - x1]
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
  (Slider vertical _ name lower upper value sliderBox) ->
    viewSlider vertical name lower upper value sliderBox
  (TextBox scaleX' scaleY' txt) -> viewText scaleX' scaleY' txt box'

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

viewSlider :: Bool -> Text -> Float -> Float -> Float -> Box -> Gloss.Picture
viewSlider vertical _ lowerBound upperBound value sliderBox =
  assignBox sliderBox (Gloss.rotate 90 sliderP)
  where position = (value - lowerBound) / (upperBound - lowerBound)
        line = Gloss.color Gloss.azure $ Gloss.line [(0, -0.5), (0,0.5)]
        mark = Gloss.color Gloss.azure
               $ Gloss.translate 0 (position - 0.5)
               $ Gloss.line [(-0.25, 0), (0.25, 0)]
        sliderP = Gloss.pictures [line, mark]

viewText :: Double -> Double -> Text -> Box -> Gloss.Picture
viewText scaleX' scaleY' txt (x1, y1, x2, y2) =
    assignBox (x1 - (x2 - x1), y1, x2, y2)
  $ Gloss.color Gloss.white
  $ Gloss.scale (double2Float scaleX' / (x2 - x1))
                (double2Float scaleY' / (y2 - y1))
  $ Gloss.text (Text.unpack txt)

assignBox :: Box -> Gloss.Picture -> Gloss.Picture
assignBox (x1, y1, x2, y2) =
  Gloss.translate moveX moveY . Gloss.scale scaleX' scaleY'
  where moveX = ((x1 + x2 - 1) / 2.0) 
        moveY = ((y1 + y2 - 1) / 2.0) 
        scaleX' = (x2 - x1)
        scaleY' = (y2 - y1)

data Direction = Horizontal | Vertical
  deriving (Show, Eq)

-- viewSeq :: [(Box, Gloss.Picture)] -> Gloss.Picture
-- viewSeq elems = Gloss.pictures $ fmap (uncurry assignBox) elems

