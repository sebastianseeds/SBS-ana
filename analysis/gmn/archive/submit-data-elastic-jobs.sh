#!/bin/bash

# swif2 job submit for data_elastic.C for gmn analysis

#Usage
#./submit-data-elastic-jobs.sh <int: kinematic> <int: sbs magnetic field in percent>

kine=$1
mag=$2

work_flow='elastic_select'

swif2 create $work_flow

jobname=$kine'-'$mag

script='/w/halla-scshelf2102/sbs/seeds/ana/analysis/gmn/run-data-elastic.sh'

#cd $SWIF_WORK_DIRECTORY
cd /w/halla-scshelf2102/sbs/seeds/ana/analysis/gmn

echo -e "\n Submitting " $script $kine $mag "\n"

swif2 add-job -workflow $work_flow -partition production -name elasticSelection-$jobname -cores 1 -disk 20GB -ram 2500MB $script $kine $mag

# run the workflow and then print status
swif2 run $work_flow
echo -e "\n Getting workflow status.. \n"
swif2 status $work_flow
