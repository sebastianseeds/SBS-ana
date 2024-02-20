#!/bin/sh

# SSeeds - 2.2.24 - shell script to run nTPE kinematics through barebones parsing

pass=$1
cluster_method=$2

nkine=3
kine=(4 8 9)

for ((i=0; i<$nkine; i++))
do

    root -l -b -q 'parse_barebones.C('${kine[i]}','$pass','$cluster_method',false)'

wait

done
