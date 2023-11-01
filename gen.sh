#! /bin/bash

# creates a movie of the pdf of the wave function in the clips folder and
# outputs the potential function into the .used/datafiles/potential.dat file

make
if [ $? -eq 0 ]; then
    echo programs compiled...
else
    exit 1
fi
make clean

.used/otherexecs/datagen < .used/parameters.txt
if [ $? -eq 0 ]; then
    echo data generated...
else
    exit 1
fi

cd .used/pythonscripts/
rm images/*
python3 waveimgen.py
if [ $? -eq 0 ]; then
    echo images generated...
else
    exit 1
fi

if [ $# -eq 1 ]; then
    clipname=$1
else
    clipname='basicprop'
fi

python3 clipstitch.py $clipname
if [ $? -eq 0 ]; then
    echo clip stitched
fi
