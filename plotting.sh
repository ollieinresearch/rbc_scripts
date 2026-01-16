#!/bin/bash

#SBATCH --output=test_trill.out
#SBATCH --time=11:59:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=192
#SBATCH --mem=0
#SBATCH --job-name=example
#SBATCH --account=def-goluskin


source "params.conf"
source $SETUP_PATH
source $WORKFLOW_PATH

print_params
ENV_PATH="$SCRATCH/ded3_vis"
load_modules

plot_snapshots