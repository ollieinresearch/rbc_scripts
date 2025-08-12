#!/bin/bash

#SBATCH --output=job_sim.out
#SBATCH --time=11:59:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=24
#SBATCH --mem=30G

#SBATCH --account=def-goluskin

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ollietheengineer@uvic.ca

# Load the required modules
module purge
module load StdEnv/2020
module load fftw-mpi mpi4py hdf5-mpi python

# For our virtual environment
env=$SLURM_TMPDIR/env

# Create the virtual environment on each node: 
srun --ntasks $SLURM_NNODES --tasks-per-node=1 bash << EOF
virtualenv --no-download $env
source $env/bin/activate

pip install --no-index --upgrade pip
pip install --no-index -r /home/ollie/requirements.txt
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
PR_EXP=-3.5
# Vertical resolution
RES=192
# Timestep- if using fixed timestep this matters. Otherwise just leave
# sufficiently small that the simulation won't blow up in 25 iterations
DT=0.0002
# Dimensionless time to run the simulation for
SIM_TIME=175
# Method for timestepping. Can be RK222, RK443, CNAB2, MCNAB2, SBDF4
STEPPER=RK443
# Aspect ratio
GAM=2
# Use an initial condition from a previous run? 1=yes, 0=no
IC=1
################################################################################

# path to all python scripts for simulations; change as needed
PATH_TO_SCRIPTS="/project/def-goluskin/ollie/scripts/2d"

# specify desired time for initial condition
if [ $IC -eq 1 ]; then
  IC_ARRAY=($(python3 $PATH_TO_SCRIPTS/initial_condition.py 0 --file=$PWD/restart/restart.h5))
  TOTAL_TIME=$(echo "$SIM_TIME+${IC_ARRAY[1]}" | bc)
  IND=${IC_ARRAY[0]}

else
  echo "Not running with a prior simulation's initial conditions; starting fresh!"
  rm -rf {analysis,state,preliminary_outputs,snapshots,outputs,nu_analysis,restart}

  TOTAL_TIME=$SIM_TIME
  IND=-1
fi

srun python3 $PATH_TO_SCRIPTS/rayleigh_benard_script.py --Ra=$RA --Pr_exp=$PR_EXP --res=$RES --dt=$DT --sim_time=$TOTAL_TIME --stepper=$STEPPER --Gamma=$GAM --index=$IND --basepath=$PWD --cfl

# Post processing
if [ -f "analysis/analysis.h5" ]; then
    rm analysis/analysis.h5
    echo "removing existing analysis file"
fi

if [ -f "nu_analysis/nu_analysis.h5" ]; then
    rm nu_analysis/nu_analysis.h5
    echo "removing existing nu_analysis file"
fi

# merge analysis data into one file to plot time series data
srun python3 -m dedalus merge_procs analysis --cleanup
srun python3 -m dedalus merge_procs nu_analysis --cleanup
srun python3 -m dedalus merge_sets analysis/analysis.h5 analysis/*.h5
srun python3 -m dedalus merge_sets nu_analysis/nu_analysis.h5 nu_analysis/*.h5


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

srun python3 $PATH_TO_SCRIPTS/prelim_test.py $PWD/nu_analysis/nu_analysis.h5 --output=$PWD/preliminary_outputs

srun python3 $PATH_TO_SCRIPTS/nu_analysis.py $PWD/nu_analysis/nu_analysis.h5 --time=7481 --output=$PWD/outputs

sbatch multinode_run_sim.sh