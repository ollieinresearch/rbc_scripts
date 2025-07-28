#!/bin/bash

#SBATCH --output=job_movie.out
#SBATCH --time=00:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --mem=0
#SBATCH --account=def-goluskin

# Load the required modules
module purge
module load CCEnv arch/avx512 StdEnv/2020
module load python/3.10.2 mpi4py/3.1.3 fftw-mpi/3.3.8 hdf5-mpi/1.10.6 scipy-stack

# For our virtual environment
env=$SLURM_TMPDIR/env
PATH_TO_SCRIPTS="$BBUFFER/rbc_scripts/3d/visualization"
CACHE_DIR=$BBUFFER/.cache

# Dedalus performance tip!
export OMP_NUM_THREADS=1
export NUMEXPR_MAX_THREADS=1
export PIP_CACHE_DIR=$CACHE_DIRpip
export MPLCONFIGDIR=$CACHE_DIR/matplotlib

# Create the virtual environment on each node: 
virtualenv --no-download $env
source $env/bin/activate

pip install --no-index --upgrade --cache-dir=$CACHE_DIR pip
pip install --no-index --cache-dir=$CACHE_DIR -r $PATH_TO_SCRIPTS/requirements_vis2.txt


# merge snapshot processes into set files - don't need one big file here since plotting slices
# srun python3 -m dedalus merge_procs snapshots --cleanup

mkdir $PWD/visualization
srun python3 $PATH_TO_SCRIPTS/isos_test.py $PWD/snapshots/*.h5 --basepath=$PWD
ffmpeg -y -r 15 -pattern_type glob -i 'visualization/*.png' -threads 40 -pix_fmt yuv420p visualization/movie.mp4
