#!/bin/bash

#SBATCH --output=job.out
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --mem=16G

#SBATCH --account=def-goluskin

module purge
module load StdEnv/2020
module load fftw-mpi mpi4py hdf5-mpi python

env=$SLURM_TMPDIR/env


virtualenv --no-download $env
source $env/bin/activate
pip install --no-index --upgrade pip
pip install --no-index -r /home/ollie/requirements.txt


export OMP_NUM_THREADS=1
export NUMEXPR_MAX_THREADS=1


# Filepath for specific resolution
FILEPATH="$(find . -mindepth 1 -maxdepth 1 -type d | sed 's|^\./||')"

#path to all python scripts for simulations; change as needed
PATH_TO_SCRIPTS="/project/def-goluskin/ollie/scripts/2d"

# merge snapshot processes into set files - don't need one big file here since plotting slices
python3 $PATH_TO_SCRIPTS/analysis.py ${FILEPATH}/analysis/analysis.h5 --time=200 --output=$FILEPATH/output_tests
