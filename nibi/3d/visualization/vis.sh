#!/bin/bash
#SBATCH --output=job_movie.out
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=96
#SBATCH --mem=0
#SBATCH --account=def-goluskin

# Load the required modules
module purge
module load StdEnv/2020
module load python/3.10.2 mpi4py/3.1.3 fftw-mpi/3.3.8 hdf5-mpi/1.10.6 scipy-stack

# For our virtual environment
env=$SLURM_TMPDIR/env

#path to all python scripts for simulations; change as needed
PATH_TO_SCRIPTS="/scratch/ollie/rbc_scripts/3d/visualization"

# Dedalus performance tip!
export OMP_NUM_THREADS=1
export NUMEXPR_MAX_THREADS=1

# Create the virtual environment on each node: 
virtualenv --no-download $env
source $env/bin/activate

pip install --no-index --upgrade pip
pip install --no-index -r $PATH_TO_SCRIPTS/requirements_vis.txt

mkdir $PWD/visualization
srun python3 $PATH_TO_SCRIPTS/isos_test.py $PWD/snapshots/*.h5 --basepath=$PWD
ffmpeg -y -r 15 -pattern_type glob -i 'visualization/*.png' -threads 96 -pix_fmt yuv420p visualization/movie.mp4