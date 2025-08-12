#!/bin/bash
#SBATCH --output=job_movie.out
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem=0
#SBATCH --account=def-goluskin

# Path to all python scripts for simulations; change as needed
PATH_TO_SCRIPTS="/home/ollie/scratch/rbc_scripts"
SCRIPTS_3D="/home/ollie/scratch/rbc_scripts/3d"
PATH_TO_ENV="/home/ollie/scratch/dedalus"

################################################################################

module --force purge
module load StdEnv/2020
module load python/3.10.2 mpi4py/3.1.3 fftw-mpi/3.3.8 hdf5-mpi/1.12.1 scipy-stack/2023b

source $PATH_TO_ENV/bin/activate

# Dedalus performance tip!
export OMP_NUM_THREADS=1
export NUMEXPR_MAX_THREADS=1


srun python3 $PATH_TO_SCRIPTS/isos.py $PWD/snapshots/*.h5 --basepath=$PWD
ffmpeg -y -r 15 -i $PWD/visualization/write_%06d.png -threads 48 -pix_fmt yuv420p visualization/movie.mp4