#!/bin/bash

#SBATCH --output=job_vis.out
#SBATCH --time=00:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --mem=0
#SBATCH --account=def-goluskin

# Path to all python scripts for simulations; change as needed
PATH_TO_SCRIPTS="$SCRATCH/rbc_scripts/2d"
PATH_TO_ENV="$SCRATCH/dedalus"


################################################################################


# Load the required modules
module load StdEnv/2020
module load python/3.10.2 mpi4py/3.1.3 fftw-mpi/3.3.8 hdf5-mpi/1.12.1 scipy-stack/2023b ffmpeg

source $PATH_TO_ENV/bin/activate;

# Dedalus performance tip!
export OMP_NUM_THREADS=1
export NUMEXPR_MAX_THREADS=1

# merge snapshot processes into set files - don't need one big file here since plotting slices
srun python3 -m dedalus merge_procs snapshots --cleanup

srun python3 $PATH_TO_SCRIPTS/max_vort.py snapshots/*.h5
python3 $PATH_TO_SCRIPTS/max_vort_2.py $PWD/snapshots

mv="$(cat snapshots/max_vort.txt)"

srun python3 $PATH_TO_SCRIPTS/plot_slices.py snapshots/*.h5 --output=frames --max_vort=$mv

# Piece frames together!
rm frames/movie.mp4
ffmpeg -y -r 24 -pattern_type glob -i 'frames/*.png' -threads 192 -pix_fmt yuv420p frames/movie.mp4
