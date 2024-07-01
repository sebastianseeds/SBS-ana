#!/bin/sh

# SSeeds - shell script to run nTPE kinematics through plotwithcuts best cluster

nkine=4
kine=(8 8 8 9)
mag=(50 70 100 70)

for ((i=0; i<$nkine; i++))
do

    root -l -b -q 'rangeStability.C('${kine[i]}','${mag[i]}',2,8,0.1,2,true,true,false,true,true)'

wait

done
