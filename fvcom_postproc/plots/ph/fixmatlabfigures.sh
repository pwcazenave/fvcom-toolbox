#!/bin/bash

# Use ImageMagick to convert ugly MATLAB pdfs to beautiful pngs.

inFiles=("$@")

for ((i=0; i<${#inFiles[@]}; i++)); do
    echo -n "Fixing ${inFiles[i]}... "
    convert \
        -trim \
        -density 600x600 \
        "${inFiles[i]}" \
        "${inFiles[i]%.*}".png
    echo "done."
done
