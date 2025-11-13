#!/bin/bash
#SBATCH --output=job_post.out
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --mem=0
#SBATCH --account=def-goluskin


# Path to all python scripts for simulations; change as needed
PATH_TO_SCRIPTS="$SCRATCH/rbc_scripts"
SCRIPTS_3D="$SCRATCH/rbc_scripts/3d"
PATH_TO_ENV="$SCRATCH/ded3test"

# Load the required modules
module load StdEnv/2023
ml python/3.11.5 mpi4py/3.1.4 fftw-mpi/3.3.10 hdf5-mpi/1.14.2
source /scinet/vast/etc/vastpreload-openmpi.bash
which mpirun
source $PATH_TO_ENV/bin/activate;
export MPLCONFIGDIR=$SCRATCH/.cache/matplotlib
# Dedalus performance tip!
export OMP_NUM_THREADS=1
export NUMEXPR_MAX_THREADS=1


srun python3 merge.py
