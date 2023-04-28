#!/bin/sh

# SSeeds - 3.21.23 - simple shell script to run momentum_calib_simu.C

## Usage
#./run_pcalibsimu.sh <kine>

kine=$1

root -l 'momentum_calib_simu.C('$kine')'

