#!/bin/bash

#SBATCH --output=job_movie.out
#SBATCH --time=02:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --mem=0
#SBATCH --account=def-goluskin

# Path to all python scripts for simulations; change as needed
PATH_TO_SCRIPTS="$SCRATCH/rbc_scripts"
SCRIPTS_3D="$SCRATCH/rbc_scripts/3d/visualization"
PATH_TO_ENV="$SCRATCH/ded3_vis"

# Load the required modules
ml StdEnv/2023
ml python/3.11.5 mpi4py/4.0.3 fftw-mpi/3.3.10 scipy-stack/2023b vtk/9.3.0 ffmpeg/7.1.1 
ml hdf5-mpi/1.14.4

source /scinet/vast/etc/vastpreload-openmpi.bash

source $PATH_TO_ENV/bin/activate;

# Dedalus performance tip!
export OMP_NUM_THREADS=1
export NUMEXPR_MAX_THREADS=1

# merge snapshot processes into set files - don't need one big file here since plotting slices
# srun python3 -m dedalus merge_procs snapshots --cleanup

nu=$(grep "Final Nusselt number:" outputs/info.txt | awk -F': ' '{print $2}')

mkdir $PWD/visualization

mpirun -n 32 python3 $SCRIPTS_3D/plotting_v3.py $PWD/snapshots/*.h5 --basepath=$PWD --nu=$nu
ffmpeg -y -r 60 -pattern_type glob -i 'visualization/*.png' -threads 32 -pix_fmt yuv420p visualization/movie.mp4
