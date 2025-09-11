#!/bin/bash

#SBATCH --output=test.out
#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --mem=0
#SBATCH --job-name=example
#SBATCH --account=def-goluskin
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ollietheengineer@uvic.ca



################################################################################
# Simulation parameters
START_TIME=0
# Rayleigh number
RA=1e5
# Exponent of 10 for Prandtl number. (ie if Pr=1, PR_EXP=0)
PR_EXP=0
# Vertical resolution
RES=96
# Timestep- if using fixed timestep this matters. Otherwise just leave
# sufficiently small that the simulation won't blow up in 25 iterations
DT=0.005
# Dimensionless time to run the simulation for
SIM_TIME=40
# Method for timestepping. Can be RK222, RK443, CNAB2, MCNAB2, SBDF4
STEPPER=RK222
# Aspect ratio in x and y resp.
LX=2
LY=2
# Provide the size of the square 2D process mesh; note that if your resolution is given by R,
# you have a grid of size Lx*R x Ly*R x R, and if your mesh is of size [n,m], then on each core you compute
# (Lx*R/n)x(Ly*R/m) pencils of length R. Most efficient when MESHX=MESHY.
MESHX=16
MESHY=12

# Use initial condition? (1=yes, 0=no) - 0 deletes any data in the folder from a different run!
IC=0
################################################################################


# Path to all python scripts for simulations; change as needed
PATH_TO_SCRIPTS="$SCRATCH/rbc_scripts"
SCRIPTS_3D="$SCRATCH/rbc_scripts/test"
PATH_TO_ENV="$SCRATCH/ded3test"

# Load the required modules
ml StdEnv/2023
ml python/3.11.5 mpi4py/3.1.4 fftw-mpi/3.3.10 hdf5-mpi/1.14.2
ml


source $PATH_TO_ENV/bin/activate;


# Dedalus performance tip!
export OMP_NUM_THREADS=1
export NUMEXPR_MAX_THREADS=1


# Specify desired time for initial condition - 0 implies most recent
if [ $IC -eq 1 ]; then
  IC_ARRAY=($(python3 $PATH_TO_SCRIPTS/initial_condition.py $START_TIME --file=$PWD/restart/restart.h5))
  TOTAL_TIME=$(echo "$SIM_TIME+${IC_ARRAY[1]}" | bc)
  IND=${IC_ARRAY[0]}
else
  echo "Not running with a prior simulation's initial conditions; starting fresh!"
  rm -rf restart
  # rm -rf {analysis,state,preliminary_outputs,snapshots,outputs,field_analysis,restart}

  TOTAL_TIME=$SIM_TIME
  IND=-1
fi

echo "Right before the first srun"
srun python3 $SCRIPTS_3D/gpt_test.py --Ra=$RA --Pr_exp=$PR_EXP --res=$RES --dt=$DT --sim_time=$TOTAL_TIME --index=$IND --basepath=$PWD --stepper=$STEPPER --Lx=$LX --Ly=$LY --meshx=$MESHX --meshy=$MESHY --cfl

exit


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

mkdir res_check
mkdir res_check_3d

srun python3 $PATH_TO_SCRIPTS/power.py $PWD/state/*.h5
srun python3 $SCRIPTS_3D/power.py $PWD/state/*.h5
ffmpeg -y -r 15 -pattern_type glob -i 'res_check/*.png' -threads 96 -pix_fmt yuv420p res_check/movie.mp4
ffmpeg -y -r 15 -pattern_type glob -i 'res_check_3d/*.png' -threads 96 -pix_fmt yuv420p res_check_3d/movie.mp4
