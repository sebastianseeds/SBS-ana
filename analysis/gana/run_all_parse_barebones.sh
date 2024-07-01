#!/bin/sh

# SSeeds - 2.2.24 - shell script to run all kinematics through barebones parsing

pass=$1
cluster_method=$2
verbose=$3

nkine=6
kine=(4 7 11 14 8 9)

for ((i=0; i<$nkine; i++))
do

    root -l -b -q 'parse_barebones.C('${kine[i]}','$pass','$cluster_method','$verbose')'

wait

done
