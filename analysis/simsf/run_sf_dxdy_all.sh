#!/bin/sh

# SSeeds - 11.9.23 - shell script to run all scale field settings from replayed g4sbs data to build dx plots

nsf=11
sf=(0 10 20 30 40 50 60 70 80 90 100)

for ((i=0; i<$nsf; i++))
do

    root -l -b -q 'sf_dxdy.C('${sf[i]}')'

wait

done
