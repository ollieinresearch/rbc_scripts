#!/bin/bash

#SBATCH --output=post_process.out
#SBATCH --time=02:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --mem=0
#SBATCH --job-name=raXprXeX_X_post
#SBATCH --account=def-goluskin


source "params.conf"
source $SETUP_PATH
source $WORKFLOW_PATH

print_params
load_modules

post_process
res_check