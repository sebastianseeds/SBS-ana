#!/bin/sh

# SSeeds - 2.2.24 - shell script to run nTPE kinematics through barebones parsing

nkine=3
kine=(4 8 9)

for ((i=0; i<$nkine; i++))
do

    root -l -b -q 'parse_barebones.C('${kine[i]}',2,3,true,false,false,true,true,false,false)'

wait

done
