#!/bin/sh

# SSeeds - 2.3.24 - shell script to run nTPE kinematics through barebones parsing

pass=$1

nkine=4
kine=(4 4 8 9)
mag=(30 50 70 70)

for ((i=0; i<$nkine; i++))
do

    root -l -b -q 'plotWithCuts_p2.C('${kine[i]}','${mag[i]}','$pass',false)'

wait

done
