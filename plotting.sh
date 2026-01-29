#!/bin/bash

#SBATCH --output=test_trill.out
#SBATCH --time=02:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --mem=0
#SBATCH --job-name=example
#SBATCH --account=def-goluskin


source "params.conf"
source $SETUP_PATH
source $WORKFLOW_PATH

print_params
load_modules_plotting
ffmpeg -y -r 30 -pattern_type glob -i 'res_check_3d/*.png' -threads 32 -pix_fmt yuv420p $PWD/res_check_3d/movie.mp4
plot_snapshots