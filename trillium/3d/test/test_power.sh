#!/bin/bash

PATH_TO_SCRIPTS="/home/ollie/links/scratch/rbc_scripts"

# Load the required modules
module purge
module load StdEnv/2020
module load fftw-mpi mpi4py hdf5-mpi python/3.10.2

# Dedalus performance tip!
export OMP_NUM_THREADS=1
export NUMEXPR_MAX_THREADS=1

# For our virtual environment
env=$SLURM_TMPDIR/env

# Create the virtual environment on each node: 
virtualenv --no-download $env
source $env/bin/activate

pip install --no-index --upgrade pip
pip install --no-index -r $PATH_TO_SCRIPTS/requirements.txt

mkdir res_check

srun python3 power_test.py

ffmpeg -y -r 15 -i res_check/write_%06d.png -threads 192 -pix_fmt yuv420p res_check/movie.mp4