module Main where

import qualified Graphics.Gloss as Gloss
import qualified Electra2Shadow

main :: IO ()
main = Gloss.play
         Electra2Shadow.displayMode
         Electra2Shadow.backgroundColor
         Electra2Shadow.fps
         Electra2Shadow.initialWorld
         Electra2Shadow.view
         Electra2Shadow.updateInputs
         Electra2Shadow.updateTime
