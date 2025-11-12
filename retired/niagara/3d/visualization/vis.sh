#!/bin/bash

#SBATCH --output=job_movie.out
#SBATCH --time=02:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --mem=0
#SBATCH --account=def-goluskin

# Path to all python scripts for simulations; change as needed
PATH_TO_SCRIPTS="$BBUFFER/rbc_scripts"
SCRIPTS_3D="$BBUFFER/rbc_scripts/3d"
PATH_TO_ENV="$BBUFFER/vis_env"

# Load the required modules
module load CCEnv arch/avx512 StdEnv/2020
module load python/3.10.2 mpi4py/3.1.3 fftw-mpi/3.3.8 hdf5-mpi/1.12.1 scipy-stack/2023b

source $PATH_TO_ENV/bin/activate;

# Dedalus performance tip!
export OMP_NUM_THREADS=1
export NUMEXPR_MAX_THREADS=1

# merge snapshot processes into set files - don't need one big file here since plotting slices
# srun python3 -m dedalus merge_procs snapshots --cleanup

nu=$(grep "Final Nusselt number:" outputs/info.txt | awk -F': ' '{print $2}')

mkdir $PWD/visualization

srun python3 $PATH_TO_SCRIPTS/isos_test.py $PWD/snapshots/*.h5 --basepath=$PWD --nu=$nu
ffmpeg -y -r 15 -pattern_type glob -i 'visualization/*.png' -threads 40 -pix_fmt yuv420p visualization/movie.mp4
