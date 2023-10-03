#!/bin/sh

# shell script to run data_elastic.C for use with submit-data-elastic-jobs.sh

## Usage
#./run-data-elastic.sh <kinematic> <magnetic-field-setting>

kine=$1
mag=$2

cd /w/halla-scshelf2102/sbs/seeds/ana/analysis/gmn

root -l -b -q 'data_elastic.C('$kine','$mag')'
