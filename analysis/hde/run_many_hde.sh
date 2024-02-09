#!/bin/sh

# SSeeds - 10.23.23 - shell script to run several hde jobs on farm

nloop=3
kine=(4 8 9)
esigma=(1 1 1)
dsigma=(3 3 3)
mag=(30 70 70)

for ((i=0; i<$nloop; i++))
do

    root -l -b -q 'hde_dataloop.C('${kine[i]}','${mag[i]}','${dsigma[i]}','${esigma[i]}')'

wait

done
