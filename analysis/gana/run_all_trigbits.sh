#!/bin/sh

# SSeeds - 4.15.24 - shell script to run all kinematics through all_trigbits.C

nkine=6
kine=(4 7 11 14 8 9)

for ((i=0; i<$nkine; i++))
do

    root -l -b -q 'all_trigbits.C('${kine[i]}')'

    wait

done
