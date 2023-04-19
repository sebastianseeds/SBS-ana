#!/bin/sh

# SSeeds - 4.1.23 - shell script to run mc hcal detection efficiency script quickly

## Usage
#./run_mchde.sh <kine> <field> <iter>

kine=$1
field=$2

root -l -b -q 'simhde.C('$kine','$field')'

