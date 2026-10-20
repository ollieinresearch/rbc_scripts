#!/bin/bash

script_path=$SCRATCH/rbc_scripts/shell
cp $script_path/post.sh .
cp $script_path/submit.sh .
if [[ ${CC_CLUSTER:-} == trillium ]]; then script=run_sim_trillium.sh
else script=run_sim.sh; fi
cp $script_path/$script .

if [ ! -f "$PWD/params.conf" ]; then
    cp $script_path/params.conf .
fi