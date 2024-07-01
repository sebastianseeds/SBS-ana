#!/bin/bash

# sseeds 2.2.24

# Quick check on Q2 reonstruction from p/n angles

echo Generating expected elastic parameters for GMn kinematic SBS $1..

if [ $1 -eq 1 ]
then 
    root -l -q -b 'GMnQ2recon.C(1.92,51.,1.85,1.,33.5,13.5)'

elif [ $1 -eq 4 ]
then
    root -l -q -b 'GMnQ2recon.C(3.7393,36.,1.8,3.,31.9,11.)'

elif [ $1 -eq 7 ]
then
    root -l -q -b 'GMnQ2recon.C(7.9308,40.,1.85,10.,16.1,14.)'

elif [ $1 -eq 8 ]
then
    root -l -q -b 'GMnQ2recon.C(5.9828,26.5,2.,4.45,29.4,11.)'

elif [ $1 -eq 9 ]
then
    root -l -q -b 'GMnQ2recon.C(4.0268,49.,1.55,4.45,22.,11.)'

elif [ $1 -eq 11 ]
then
    root -l -q -b 'GMnQ2recon.C(9.889,42.,1.55,13.6,13.3,14.5)'

elif [ $1 -eq 14 ]
then
    root -l -q -b 'GMnQ2recon.C(5.9828,46.5,1.85,7.5,17.3,14.)'

else
    echo Invalid GMn kinematic entered.

fi
