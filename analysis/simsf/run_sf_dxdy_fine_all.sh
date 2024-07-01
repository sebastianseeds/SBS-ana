#!/bin/sh

# SSeeds - 11.10.23 - shell script to run all fine tune scale field settings from replayed g4sbs data to build dx plots

nsf=10
sf=(62 63 64 66 67 68 69 71 72 73)
sf_fine=(2653 5395 8137 879 3621 6363 9106 1848 4590 7332)

for ((i=0; i<$nsf; i++))
do

    root -l -b -q 'sf_dxdy.C('${sf[i]}','1','${sf_fine[i]}')'

wait

done
