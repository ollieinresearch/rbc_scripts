#!/bin/bash


cp $SCRATCH/rbc_scripts/run_sim.sh .
cp $SCRATCH/rbc_scripts/post_process.sh .
cp $SCRATCH/rbc_scripts/plotting.sh .

if [ ! -f "$PWD/params.conf" ]; then
    cp $SCRATCH/rbc_scripts/default.conf .
fi