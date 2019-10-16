#!/usr/bin/env fish

set target $argv[1]

echo "Copying to "$target

rsync -avzP (stack exec -- which electra2shadow-exe) $target
rsync -avzP data $target
rsync -avzP config.example $target
