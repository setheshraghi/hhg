#! /bin/bash

# creates and opens images of plots that allow verification for the correctness
# of the datagen.c file

make
if [ $? -eq 0 ]; then
    echo programs compiled...
else
    exit 1
fi
make clean

.used/otherexecs/verify
if [ $? -eq 0 ]; then
    echo data generated...
else
    exit 1
fi

./vfiguresgen.plt
if [ $? -eq 0 ]; then
    echo figures generated
else
    exit 1
fi

open .used/figures/H.png .used/figures/norm.png &&
open .used/figures/p.png .used/figures/x.png
