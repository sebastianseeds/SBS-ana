#!/bin/sh

# SSeeds - 4.23.23 - shell script to run all kinematics through parsing

## Usage
#./run_all_parse.sh <int: e' momentum method>

epm=$1

nkine=6
kine=(4 7 11 14 8 9)

for ((i=0; i<$nkine; i++))
do

    root -l -b -q 'parse_sh.C('${kine[i]}','$epm')'

wait

done
