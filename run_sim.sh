#!/bin/bash

#SBATCH --output=job_sim.out
#SBATCH --time=11:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --mem=0
#SBATCH --job-name=example
#SBATCH --account=def-goluskin


source "params.conf"
source $SETUP_PATH
source $WORKFLOW_PATH

print_params
load_modules
#export HDF5_USE_FILE_LOCKING='False'
#export HDF5_USE_FILE_LOCKING='True'

create_restart
determine_ic
run_sim
post_process
res_check