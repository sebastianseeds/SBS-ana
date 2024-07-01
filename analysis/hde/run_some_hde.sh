#!/bin/sh

# SSeeds - 10.23.23 - shell script to run several hde jobs on farm

nloop=4
kine=(8 8 9 9)
sigma=(2 3 2 3)
mag=(70 70 70 70)

for ((i=0; i<$nloop; i++))
do

    root -l -b -q 'hde_dataloop.C('${kine[i]}','${mag[i]}','${sigma[i]}')'

wait

done
