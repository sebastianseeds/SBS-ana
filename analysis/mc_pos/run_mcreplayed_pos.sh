#!/bin/sh

# SSeeds - 4.19.23 - shell script to run mc hcal detection efficiency script quickly

## Usage
#./run_mcreplayed_pos.sh <int: kinematic> <int: e prime momentum calculation method> <bool: wide plot angle>

kine=$1
epm=$2
wide=$3

root -l 'mcreplayed_pos.C('$kine','$epm','$wide')'

