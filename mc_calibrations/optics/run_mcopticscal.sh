#!/bin/sh

# SSeeds - 4.23.23 - shell script to run mc optics calibration script quickly

## Usage
#./run_mcopticscal.sh <int: kinematic>

kine=$1
rd=$2

infilename='setup_OpticalCalib_MC_SBS'$kine'.txt'
outfilename='mcoptout_sbs'$kine'_rd'$rd'.root'

root -l -b -q 'Optics_GMN.C('$infilename','$outfilename')'

