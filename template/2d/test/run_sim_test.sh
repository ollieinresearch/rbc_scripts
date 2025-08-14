#!/bin/bash


# Load the required modules
module purge
module load StdEnv/2020
module load fftw-mpi mpi4py hdf5-mpi python

source /home/ollie/ded_env/bin/activate

# Dedalus performance tip
export OMP_NUM_THREADS=1
export NUMEXPR_MAX_THREADS=1


# Run Dedalus RB script in 2D with specified parameters, with initial conditions
# From a previous simulation with a well established flow.
# Merge files for analysis.
# Run analysis script to produce info and plots.

################################################################################
# User specified parameters

# Rayleigh number
RA=1e5
# Exponent of 10 for Pr. (ie if Pr=1, PR_EXP=0, and Pr=0.1 -> PR_EXP=-1)
PR_EXP=-0.5
# Vertical resolution
RES=48
# Timestep- if using fixed timestep this matters. Otherwise just leave
# sufficiently small that the simulation won't blow up in 25 iterations
DT=0.01
# Dimensionless time to run the simulation for
SIM_TIME=0.5
# Method for timestepping. Can be RK222, RK443, CNAB2, MCNAB2, SBDF4
STEPPER=RK222
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
  rm -rf {analysis,state,preliminary_outputs,snapshots,outputs}

  TOTAL_TIME=$SIM_TIME
  IND=-1
fi

srun python3 $PATH_TO_SCRIPTS/rayleigh_benard_script_test.py --Ra=$RA --Pr_exp=$PR_EXP --res=$RES --dt=$DT --sim_time=$TOTAL_TIME --stepper=$STEPPER --Gamma=$GAM --index=$IND --basepath=$PWD