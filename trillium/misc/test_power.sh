#!/bin/bash
#SBATCH --output=job_post.out
#SBATCH --job-name=r1e7_pr3e-1_power
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --mem=0
#SBATCH --account=def-goluskin

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
PATH_TO_SCRIPTS=$BBUFFER/scripts

# Create the virtual environment on each node: 
virtualenv --no-download $env
source $env/bin/activate

pip install --no-index --upgrade --cache-dir=$BBUFFER/.cache pip
pip install --no-index --cache-dir=$BBUFFER/.cache -r $PATH_TO_SCRIPTS/requirements.txt

srun python3 -m dedalus merge_procs state --cleanup

mkdir $PWD/res_check

srun python3 $PATH_TO_SCRIPTS/power.py $PWD/state/*.h5
ffmpeg -y -r 15 -pattern_type glob -i 'res_check/*.png' -threads 40 -pix_fmt yuv420p res_check/movie.mp4
