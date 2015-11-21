#!/bin/bash

# This bash script converts all GHOST binary output files in a 
# directory to VDF format. VDF is the file format used by VAPOR, 
# a software for 3D interactive visualization: 
# https://www.vapor.ucar.edu

# File name for VAPOR VDF file and number of levels used to compress.
VDF='output.vdf'
LEV='3'

# Path to GHOST binary files
DIR='../../3D/bin/outs'

# Spatial resolution, number of snapshots, and variables to convert.
# The names of the variables must be the same as the name of the 
# output files from GHOST (or of the postprocessed files).
DIMS='128x128x128'
NUMS='10'
VARS='vx:vy:vz'

# Do not edit below this line
vdfcreate -dimension $DIMS -level $LEV -numts $NUMS -vars3d $VARS $VDF

SIZESTR=${#DIR}
echo $SIZESTR

for FILE in $DIR/*.out
do
   STR=${FILE#$DIR/}
   VAR=${STR%%.*}
   AUX=${STR#$VAR.}
   NUM=${AUX%%.*}
   raw2vdf -level -1 -ts $NUM -varname $VAR $VDF $FILE
done
