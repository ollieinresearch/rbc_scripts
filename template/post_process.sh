#!/bin/bash
#SBATCH --output=job_post.out
#SBATCH --time=01:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem=0
#SBATCH --account=def-goluskin

################################################################################

# Path to all python scripts for simulations; change as needed
PATH_TO_SCRIPTS="/home/ollie/scratch/rbc_scripts"
SCRIPTS_3D="/home/ollie/scratch/rbc_scripts/3d"
PATH_TO_ENV="/home/ollie/scratch/dedalus"

module --force purge
module load StdEnv/2020
module load python/3.10.2 mpi4py/3.1.3 fftw-mpi/3.3.8 hdf5-mpi/1.12.1 scipy-stack/2023b

source $PATH_TO_ENV/bin/activate

# Dedalus performance tip!
export OMP_NUM_THREADS=1
export NUMEXPR_MAX_THREADS=1

################################################################################


# Post processing
if [ -f "field_analysis/field_analysis.h5" ]; then
    rm field_analysis/field_analysis.h5
    echo "removing existing field_analysis file"
fi

if [ -f "analysis/analysis.h5" ]; then
    rm analysis/analysis.h5
    echo "removing existing analysis file"
fi

# Merge processors into sets
srun python3 -m dedalus merge_procs field_analysis --cleanup
srun python3 -m dedalus merge_procs analysis --cleanup
srun python3 -m dedalus merge_procs snapshots --cleanup
srun python3 -m dedalus merge_procs state --cleanup
srun python3 -m dedalus merge_sets field_analysis/field_analysis.h5 field_analysis/*.h5
srun python3 -m dedalus merge_sets analysis/analysis.h5 analysis/*.h5

# For deciding the restart path
RECENT=$(find state/. -maxdepth 1 -type f -exec basename {} \; | sort -V | tail -n 1)

if [ ! -d "restart" ]; then
  mkdir restart
  echo "creating directory for restart"
fi

rm -rf restart/restart.h5
ln -sv $PWD/state/$RECENT $PWD/restart/restart.h5

srun python3 $PATH_TO_SCRIPTS/analysis.py $PWD/analysis/analysis.h5 --time=0 --basepath=$PWD

srun python3 $PATH_TO_SCRIPTS/power.py $PWD/state/*.h5
ffmpeg -y -r 15 -pattern_type glob -i 'res_check/*.png' -threads 48 -pix_fmt yuv420p res_check/movie.mp4
