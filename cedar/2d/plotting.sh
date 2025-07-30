#!/bin/bash

#SBATCH --output=job_vis.out
#SBATCH --time=00:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem=0
#SBATCH --account=def-goluskin

# Load the required modules
module purge
module load StdEnv/2020
module load fftw-mpi mpi4py hdf5-mpi python ffmpeg

# For our virtual environment
env=$SLURM_TMPDIR/env

#path to all python scripts for simulations; change as needed
PATH_TO_SCRIPTS="/home/ollie/scratch/rbc_scripts/2d"
PATH_TO_GEN_SCRIPTS="/home/ollie/scratch/rbc_scripts"

# Create the virtual environment on each node: 
srun --ntasks $SLURM_NNODES --tasks-per-node=1 bash << EOF
virtualenv --no-download $env
source $env/bin/activate

pip install --no-index --upgrade pip
pip install --no-index -r $PATH_TO_GEN_SCRIPTS/requirements.txt
EOF

# Dedalus performance tip!
export OMP_NUM_THREADS=1
export NUMEXPR_MAX_THREADS=1

# Required for multiple node simulations
source $env/bin/activate;

# merge snapshot processes into set files - don't need one big file here since plotting slices
srun python3 -m dedalus merge_procs snapshots --cleanup

srun python3 $PATH_TO_SCRIPTS/max_vort.py snapshots/*.h5

mv="$(cat snapshots/max_vort.txt)"

srun python3 $PATH_TO_SCRIPTS/plot_slices.py snapshots/*.h5 --output=frames --max_vort=$mv

# Piece frames together!
rm frames/movie.mp4
ffmpeg -r 20 -i frames/write_%06d.png -threads 48 -pix_fmt yuv420p frames/movie.mp4