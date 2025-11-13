#!/bin/bash

# Load the required modules
module purge
module load CCEnv arch/avx512 StdEnv/2020
module load fftw-mpi mpi4py hdf5-mpi python/3.10.2

# Dedalus performance tip!
export OMP_NUM_THREADS=1
export NUMEXPR_MAX_THREADS=1
export PIP_CACHE_DIR=$BBUFFER/.cache
export MPLCONFIGDIR=$BBUFFER/.cache/matplotlib

# For our virtual environment
env=$SLURM_TMPDIR/env

# Create the virtual environment on each node: 
virtualenv --no-download $env
source $env/bin/activate

pip install --no-index --upgrade --cache-dir=$BBUFFER/.cache pip
pip install --no-index --cache-dir=$BBUFFER/.cache -r $BBUFFER/requirements.txt

PATH_TO_SCRIPTS=$BBUFFER/scripts/2d/test

mkdir res_check

srun python3 $PATH_TO_SCRIPTS/power_test.py

ffmpeg -y -r 15 -i res_check/write_%06d.png -threads 40 -pix_fmt yuv420p res_check/movie.mp4