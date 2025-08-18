#!/bin/bash

#SBATCH --output=job.out
#SBATCH --time=00:45:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=64G

#SBATCH --account=def-goluskin

PATH_TO_SCRIPTS="/home/ollie/links/scratch/rbc_scripts"


module purge
module load StdEnv/2020
module load fftw-mpi mpi4py hdf5-mpi python

env=$SLURM_TMPDIR/env

virtualenv --no-download $env
source $env/bin/activate
pip install --no-index --upgrade pip
pip install --no-index -r $PATH_TO_SCRIPTS/requirements.txt


export OMP_NUM_THREADS=1
export NUMEXPR_MAX_THREADS=1


srun python3 $PATH_TO_SCRIPTS/redo_analysis.py