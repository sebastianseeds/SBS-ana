#!/bin/sh

# SSeeds - 4.19.23 - shell script to run mc hcal detection efficiency script quickly

## Usage
#./run_mcreplayed_pos.sh <kine>

kine=$1

root -l -b -q 'mcreplayed_pos.C('$kine')'

