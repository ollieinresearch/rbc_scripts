#!/bin/bash
#SBATCH --output=job_post.out
#SBATCH --time=00:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --mem=0
#SBATCH --account=def-goluskin

# Load the required modules
module purge
module load CCEnv arch/avx512 StdEnv/2020
module load fftw-mpi mpi4py hdf5-mpi python/3.10.2

# For our virtual environment
env=$SLURM_TMPDIR/env

# Create the virtual environment on each node: 
virtualenv --no-download $env
source $env/bin/activate

pip install --no-index --upgrade pip
pip install --no-index -r $BBUFFER/requirements.txt

# Dedalus performance tip!
export OMP_NUM_THREADS=1
export NUMEXPR_MAX_THREADS=1
export PIP_CACHE_DIR=$BBUFFER/.cache/pip
export MPLCONFIGDIR=$BBUFFER/.cache/matplotlib

# Required for multiple node simulations
source $env/bin/activate;

#path to all python scripts for simulations; change as needed
PATH_TO_SCRIPTS="$BBUFFER/scripts/3d"

for t in 0 150 200; do
    for im in "simpson13" "trapezoid"; do
        for i in 1 2 5 10; do
            srun python3 $PATH_TO_SCRIPTS/partial_analysis.py $PWD/analysis/analysis.h5 --time=$t --basepath=$PWD --freq=$i --int_method=$im
        done
    done
done