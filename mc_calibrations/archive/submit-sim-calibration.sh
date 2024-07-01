#!/bin/bash

#sseeds 4.23.23 - simple version of the submit-simu-digi-replay script written by pdatta. Modified to submit only replays for iterative calibrations on momentum and optics, mc data.

kine=$1 #kinematic {4, 7, 11, 14, 8, 9}
tar=$2 #target {lh2, ld2}
field=$3 #sbs magnetic field percent 
njobs=$4 #number of jobs
preinit='gmn_sbs'$kine'_'$tar'_'$field'p'
workflowname='seeds_calMCreplay_sbs'$kine'_'$tar'_'$field'p'
swif2 create $workflowname
# specify a directory on volatile to store g4sbs, sbsdig, & replayed files.
# Working on a single directory is convenient safe for the above mentioned
# three processes to run smoothly.
outdirpath='/lustre19/expphy/volatile/halla/sbs/seeds/sim042123'
infilepath='/lustre19/expphy/volatile/halla/sbs/seeds/sim042123'

# Valid targets
vt1="lh2"  # liquid hydrogen
vt2="ld2"  # liquid deuterium
vt3="parp"  # proton gun
vt4="parn"  # neutron gun

# Validating target parameter is correct
if [[ $tar != $vt1 ]] && [[ $tar != $vt2 ]] && [[ $tar != $vt3 ]] && [[ $tar != $vt4 ]] ; then
    echo -e "\n--!--\n Error: Target parameter should be in the set {lh2,ld2,parp,parn}."
    exit;
fi

# Validating the number of arguments provided
if [[ "$#" -ne 5 ]]; then
    echo -e "\n--!--\n Error: Incorrect number of arguments."
    echo -e "Script expecting 5 arguments: <kinematic> <target> <field> <njobs> \n"
    exit;
fi

for ((i=1; i<=$njobs; i++))
do
    #replay digitized data
    digireplayinfile=$infilepath'/'$preinit'_job_'$i
    digireplayjobname=$preinit'_digi_replay_job_'$i

    digireplayscript='/work/halla/sbs/seeds/shell/batch_farm/run-digi-replay.sh'
    
    swif2 add-job -workflow $workflowname -partition production -name $digireplayjobname -cores 1 -disk 5GB -ram 1500MB $digireplayscript $digireplayinfile $outdirpath $kine
done

# run the workflow and then print status
swif2 run $workflowname
echo -e "\n Getting workflow status.. [may take a few minutes!] \n"
swif2 status $workflowname
