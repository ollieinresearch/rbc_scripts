#!/bin/bash
#SBATCH --output=job_movie.out
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem=0
#SBATCH --account=def-goluskin

# Load the required modules
module purge
module load StdEnv/2023
module load python/3.10.2 mpi4py fftw-mpi hdf5-mpi vtk

# For our virtual environment
env=$SLURM_TMPDIR/env

#path to all python scripts for simulations; change as needed
PATH_TO_SCRIPTS="/home/ollie/scratch/ollie/scripts/3d/visualization"

# Create the virtual environment on each node: 
virtualenv --no-download $env
source $env/bin/activate

pip install --no-index --upgrade pip
pip install --no-index -r $PATH_TO_SCRIPTS/requirements_vis.txt

mt="$(cat snapshots/max_temp.txt)"
mv="$(cat snapshots/max_vert.txt)"

srun python3 $PATH_TO_SCRIPTS/isos.py $PWD/snapshots/*.h5 --basepath=$PWD
ffmpeg -y -r 15 -i $PWD/visualization/write_%06d.png -threads 48 -pix_fmt yuv420p visualization/movie.mp4