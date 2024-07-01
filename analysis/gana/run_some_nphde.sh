#!/bin/sh

# SSeeds - shell script to run nTPE kinematics through npHDE

nkine=6
kine=(4 4 8 8 8 9)
mag=(30 50 50 70 100 70)

for ((i=0; i<$nkine; i++))
do

    root -l -b -q 'npHDE.C('${kine[i]}','${mag[i]}',2,true)'

wait

done
