#! /bin/bash

# clip of T over t through k and data of T over k

make
if [ $? -eq 0 ]; then
    echo programs compiled...
else
    exit 1
fi
make clean

.used/otherexecs/transmission < .used/parameters.txt
if [ $? -eq 0 ]; then
    echo data generated...
else
    exit 1
fi

cd .used/pythonscripts/
rm images/*
python3 transimgen.py
if [ $? -eq 0 ]; then
    echo images generated...
else
    exit 1
fi

if [ $# -eq 1 ]; then
    python3 clipstitch.py $1
else
    python3 clipstitch.py F
fi
if [ $? -eq 0 ]; then
    echo clip stitched
fi
