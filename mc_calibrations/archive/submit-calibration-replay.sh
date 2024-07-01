#!/bin/bash

inputfile=$1 # output of sbsdig. Don't add file extension
outdirpath='/lustre19/expphy/volatile/halla/sbs/seeds/simulation'

script='/work/halla/sbs/seeds/shell/batch_farm/run-digi-replay.sh'

$script $inputfile $outdirpath
