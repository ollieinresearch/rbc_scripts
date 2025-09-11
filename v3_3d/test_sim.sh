#!/bin/bash

#SBATCH --output=test_v3.out
#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --mem=0
#SBATCH --job-name=example
#SBATCH --account=def-goluskin
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ollietheengineer@uvic.ca


# Path to all python scripts for simulations; change as needed
PATH_TO_SCRIPTS="$SCRATCH/rbc_scripts"
SCRIPTS_3D="$SCRATCH/rbc_scripts/3d"
PATH_TO_ENV="$SCRATCH/ded3test"

# Load the required modules
module load StdEnv/2023
ml python/3.11.5 mpi4py/3.1.4 fftw-mpi/3.3.10 hdf5-mpi/1.14.2

source $PATH_TO_ENV/bin/activate;

# Dedalus performance tip!
export OMP_NUM_THREADS=1
export NUMEXPR_MAX_THREADS=1


srun python3 /home/ollie/links/scratch/rbc_scripts/v3_3d/ded_example.py