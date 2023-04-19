#!/bin/sh

# SSeeds - 4.19.23 - shell script to run mc hcal detection efficiency script quickly

## Usage
#./run_mcreplayed_hde.sh <iter> <tfac>

iter=$1
tfac=$2

root -l 'mcreplayed_hde.C('$iter','$tfac')'

