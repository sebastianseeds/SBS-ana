#!/bin/sh

# SSeeds - shell script to run nTPE kinematics through plotwithcuts best cluster

nkine=6
kine=(4 4 8 8 8 9)
mag=(30 50 50 70 100 70)

for ((i=0; i<$nkine; i++))
do

    root -l -b -q 'plotWithCuts_mc.C('${kine[i]}','${mag[i]}',2,true,false,true,false,true,"alt")'

wait

done
