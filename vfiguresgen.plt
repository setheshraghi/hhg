#! /usr/local/bin/gnuplot -persist

# plots generated verification data and saves the plots into images

set terminal png size 400,300
set output ".used/figures/norm.png"
plot ".used/datafiles/verify.txt" u 1:2 t 'norm' w lp
set output ".used/figures/x.png"
plot ".used/datafiles/verify.txt" u 1:3 t 'x' w lp
set output ".used/figures/H.png"
plot ".used/datafiles/verify.txt" u 1:4 t 'H' w lp
set output ".used/figures/p.png"
plot ".used/datafiles/verify.txt" u 1:5 t 'p' w lp
