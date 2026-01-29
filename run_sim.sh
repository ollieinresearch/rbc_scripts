#!/bin/bash

#SBATCH --output=job_sim.out
#SBATCH --time=11:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --mem=0
#SBATCH --job-name=example
#SBATCH --account=def-goluskin
# Tells SLURM to send the SIGUSR1 signal 1800 seconds (30 min, or ~ 3Tb of transfer) before end of the time limit.
#Gabadoo SBATCH --signal=B:SIGUSR1@1800
#Gabadoo SBATCH --signal=B:SIGUSR1@300


source "params.conf"
source $SETUP_PATH
source $WORKFLOW_PATH
# Signal handler is dealt with in workflow

print_params
load_modules
#export HDF5_USE_FILE_LOCKING='False'
#export HDF5_USE_FILE_LOCKING='True'

#move_to_slrmtmp
create_restart
determine_ic
srun_sim
#move_from_slrmtmp
post_process
res_check

exit 0