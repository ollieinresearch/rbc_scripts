#!/bin/bash
#SBATCH --output=job_post.out
#SBATCH --job-name=r1e7_pr3e-1_power
#SBATCH --time=00:30:00
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


srun python3 -m dedalus merge_procs state --cleanup

mkdir $PWD/res_check

srun python3 $PATH_TO_SCRIPTS/power.py $PWD/state/*.h5
ffmpeg -y -r 15 -pattern_type glob -i 'res_check/*.png' -threads 48 -pix_fmt yuv420p res_check/movie.mp4
