#!/bin/bash

#SBATCH --output=job_sim.out
#SBATCH --time=11:59:00
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --mem=0
#SBATCH --account=def-goluskin


# Load the required modules
module purge
module load StdEnv/2020
module load python/3.10.2 mpi4py fftw-mpi hdf5-mpi

# For our virtual environment
env=$SLURM_TMPDIR/env
PATH_TO_SCRIPTS="/home/ollie/scratch/scripts/2d"
PATH_TO_GEN_SCRIPTS="/home/ollie/scratch/scripts"

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

# Run Dedalus RB script in 2D with specified parameters.
# Can either use an initial condition or start fresh.
# Merge files for analysis.
# Run analysis script to produce info and plots.

################################################################################
# User specified parameters

# Rayleigh number
RA=1e7
# Exponent of 10 for Pr. (ie if Pr=1, PR_EXP=0, and Pr=0.1 -> PR_EXP=-1)
PR_EXP=-4
# Vertical resolution
RES=200
# Timestep- if using fixed timestep this matters. Otherwise just leave
# sufficiently small that the simulation won't blow up in 25 iterations
DT=0.0002
# Dimensionless time to run the simulation for
SIM_TIME=200
# Method for timestepping. Can be RK222, RK443, CNAB2, MCNAB2, SBDF4
STEPPER=RK222
# Aspect ratio
GAM=2
# Use an initial condition from a previous run? 1=yes, 0=no
IC=1
################################################################################

# specify desired time for initial condition
if [ $IC -eq 1 ]; then
  IC_ARRAY=($(python3 $PATH_TO_GEN_SCRIPTS/initial_condition.py 0 --file=$PWD/restart/restart.h5))
  TOTAL_TIME=$(echo "$SIM_TIME+${IC_ARRAY[1]}" | bc)
  IND=${IC_ARRAY[0]}

else
  echo "Not running with a prior simulation's initial conditions; starting fresh!"
  rm -rf {analysis,state,preliminary_outputs,snapshots,outputs,field_analysis,restart}

  TOTAL_TIME=$SIM_TIME
  IND=-1
fi

srun python3 $PATH_TO_SCRIPTS/rayleigh_benard_script.py --Ra=$RA --Pr_exp=$PR_EXP --res=$RES --dt=$DT --sim_time=$TOTAL_TIME --stepper=$STEPPER --Gamma=$GAM --index=$IND --basepath=$PWD --cfl

# Post processing
if [ -f "analysis/analysis.h5" ]; then
    rm analysis/analysis.h5
    echo "removing existing analysis file"
fi

if [ -f "field_analysis/field_analysis.h5" ]; then
    rm field_analysis/field_analysis.h5
    echo "removing existing field_analysis file"
fi

# merge analysis data into one file to plot time series data
srun python3 -m dedalus merge_procs analysis --cleanup
srun python3 -m dedalus merge_procs field_analysis --cleanup
srun python3 -m dedalus merge_sets analysis/analysis.h5 analysis/*.h5
srun python3 -m dedalus merge_sets field_analysis/field_analysis.h5 field_analysis/*.h5


# merge snapshot processes into set files - don't need one big file here since plotting slices
srun python3 -m dedalus merge_procs snapshots --cleanup

# merge state processes into sets
srun python3 -m dedalus merge_procs state --cleanup


RECENT=$(find state/. -maxdepth 1 -type f -exec basename {} \; | sort -V | tail -n 1)


if [ ! -d "restart" ]; then
  mkdir restart
  echo "creating directory for restart"
fi

rm -rf restart/restart.h5

ln -sv $PWD/state/$RECENT $PWD/restart/restart.h5

srun python3 $PATH_TO_GEN_SCRIPTS/analysis.py $PWD/analysis/analysis.h5 --time=0 --basepath=$PWD

srun python3 $PATH_TO_GEN_SCRIPTS/power.py $PWD/state/*.h5
ffmpeg -y -r 15 -pattern_type glob -i 'res_check/*.png' -threads 48 -pix_fmt yuv420p res_check/movie.mp4
