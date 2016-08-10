#!/usr/bin/env bash

gs -sDEVICE=pngalpha -r300 -dBATCH -dNOPAUSE -q -sOutputFile=%stdout fvcom-toolbox.pdf | \
    convert -trim -transparent white -background transparent -trim - fvcom-toolbox.png
rm fvcom-toolbox.pdf
