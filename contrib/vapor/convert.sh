#!/bin/bash

# This bash script converts all GHOST binary output files in a 
# directory to VDF format. VDF is the file format used by VAPOR, 
# a software for 3D interactive visualization: 
# https://www.vapor.ucar.edu

# File name for VAPOR VDC file and number of levels used to compress.
VDC='output.vdc'
LEV='3'
# Path to GHOST binary files.
DIR="../../bin/outs"

# Spatial resolution, number of snapshots, and variables to convert.
# The names of the variables must be the same as the name of the 
# output files from GHOST (or of the postprocessed files).
DIMS="512x512x256"
NUMS="10"
VARS="vx:vy:vz"

# Do not edit below this line
vdccreate \
    -dimension "$DIMS" \
    -numts "$NUMS" \
    -vars3d "$VARS" \
    "$VDC"

for FILE in "$DIR"/*.out
do
    STR=${FILE##*/}
    VAR=${STR%%.*}

    AUX=${STR#*.}
    NUM=${AUX%%.*}

    raw2vdc \
        -lod -1 \
        -ts "$((10#$NUM))" \
        -varname "$VAR" \
        "$VDC" \
        "$FILE"
done
