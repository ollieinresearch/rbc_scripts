#!/bin/bash
#SBATCH --output=job_power.out
#SBATCH --job-name=r1e7_pr3e-1_power
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --mem=0
#SBATCH --account=def-goluskin

# Load the required modules
module purge
module load StdEnv/2020
module load python/3.10.2 mpi4py fftw-mpi hdf5-mpi

# Dedalus performance tip!
export OMP_NUM_THREADS=1
export NUMEXPR_MAX_THREADS=1

# For our virtual environment
env=$SLURM_TMPDIR/env
PATH_TO_SCRIPTS="/project/def-goluskin/ollie/scripts"

# Create the virtual environment on each node: 
virtualenv --no-download $env
source $env/bin/activate

pip install --no-index --upgrade pip
pip install --no-index -r $PATH_TO_SCRIPTS/requirements.txt

srun python3 -m dedalus merge_procs state --cleanup

mkdir $PWD/res_check_3d

srun python3 $PATH_TO_SCRIPTS/3d/power.py $PWD/state/*.h5
ffmpeg -y -r 15 -pattern_type glob -i 'res_check_3d/*.png' -threads 192 -pix_fmt yuv420p res_check_3d/movie.mp4
