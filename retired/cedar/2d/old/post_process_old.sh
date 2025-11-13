#!/bin/bash

#SBATCH --output=job.out
#SBATCH --time=02:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=32G

#SBATCH --account=def-goluskin

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ollietheengineer@uvic.ca

module purge
module load StdEnv/2020
module load fftw-mpi mpi4py hdf5-mpi python

env=$SLURM_TMPDIR/env


virtualenv --no-download $env
source $env/bin/activate
pip install --no-index --upgrade pip
pip install --no-index -r /home/ollie/requirements.txt


export OMP_NUM_THREADS=1
export NUMEXPR_MAX_THREADS=1



#FILEPATH="$(find . -mindepth 1 -maxdepth 1 -type d | sed 's|^\./||')"

#path to all python scripts for simulations; change as needed
PATH_TO_SCRIPTS="/project/def-goluskin/ollie/scripts/2d"



fp="/project/def-goluskin/RB_data/2D/ra1e05/ra1e05_pr1e-03/320x160_Gam2_ra1e05_pr1e-03_analysis"

# merge analysis data into one file to plot time series data
if [ -f "$fp/analysis.h5" ]; then
    rm $fp/analysis.h5
    echo "removing existing analysis file"
fi

srun python3 -m dedalus merge_procs $fp --cleanup
srun python3 -m dedalus merge_sets /project/def-goluskin/ollie/ollie_rb_data/redoing_test/ra1e5/pr1e-3/160_Gam2/analysis/analysis.h5 $fp/*.h5

fp="/project/def-goluskin/RB_data/2D/ra1e05/ra1e05_pr3e-02/192x96_Gam2_ra1e05_pr3e-02_analysis"

# merge analysis data into one file to plot time series data
if [ -f "$fp/analysis.h5" ]; then
    rm $fp/analysis.h5
    echo "removing existing analysis file"
fi

srun python3 -m dedalus merge_procs $fp --cleanup
srun python3 -m dedalus merge_sets /project/def-goluskin/ollie/ollie_rb_data/redoing_test/ra1e5/pr3e-2/96_Gam2/analysis/analysis.h5 $fp/*.h5
