#!/bin/sh

# SSeeds - 10.23.23 - shell script to run several hde jobs on farm

nloop=5
kine=(8 8 9 9 9)
sigma=(3 4 2 3 4)
mag=(70 70 70 70 70)

for ((i=0; i<$nloop; i++))
do

    root -l -b -q 'hde_dataloop.C('${kine[i]}','${mag[i]}','${sigma[i]}')'

wait

done
