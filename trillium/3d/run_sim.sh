#!/bin/bash

#SBATCH --output=job_sim.out
#SBATCH --time=11:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --mem=0
#SBATCH --job-name=ra6pr1e-2_48
#SBATCH --account=def-goluskin



################################################################################
# Simulation parameters
START_TIME=0
# Rayleigh number
RA=6
# Exponent of 10 for Prandtl number. (ie if Pr=1, PR_EXP=0)
PR=-2
# Vertical resolution
RES=48
# Timestep- if using fixed timestep this matters. Otherwise just leave
# sufficiently small that the simulation won't blow up in 25 iterations
DT=0.005
# Dimensionless time to run the simulation for
SIM_TIME=100
# Method for timestepping. Can be RK222, RK443, CNAB2, MCNAB2, SBDF4
STEPPER=SBDF4
# Aspect ratio in x and y resp.
LX=2
LY=2
# Provide the size of the square 2D process mesh; note that if your resolution is given by R,
# you have a grid of size Lx*R x Ly*R x R, and if your mesh is of size [n,m], then on each core you compute
# (Lx*R/n)x(Ly*R/m) pencils of length R. Most efficient when MESHX=MESHY.
MESHX=16
MESHY=12

# Use initial condition? (1=yes, 0=no) - 0 deletes any data in the folder from a different run!
IC=1



# CFL Stuff
CFL_SAFETY=0.1
CFL_THRESHOLD=0.075
CFL_CADENCE=1


# When to start averaging - post process
AVG_TIME=0
# Exponent of 10 for minimum y axis on power spectra 
YMIN=-18
# Exponent of 10 for maximum y axis on power spectra 
YMAX=-3
# Exponent of 10 for minimum y axis on power spectra 
YMIN_3D=-12
# Exponent of 10 for maximum y axis on power spectra 
YMAX_3D=6
################################################################################


# Path to all python scripts for simulations; change as needed
PATH_TO_SCRIPTS="$SCRATCH/rbc_scripts"
SCRIPTS_3D="$SCRATCH/rbc_scripts/3d"
PATH_TO_ENV="$SCRATCH/ded3"

# Load the required modules
ml StdEnv/2023
ml python/3.11.5 mpi4py/4.0.3 fftw-mpi/3.3.10 scipy-stack/2023b ffmpeg/7.1.1
ml hdf5-mpi/1.14.4

source /scinet/vast/etc/vastpreload-openmpi.bash

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
  # rm -rf {analysis,state,preliminary_outputs,snapshots,outputs,field_analysis,restart}
  rm -rf restart
  TOTAL_TIME=$SIM_TIME
  IND=-1
fi

mpirun python3 $SCRIPTS_3D/rayleigh_benard_script.py --Ra=$RA --Pr=$PR --nz=$RES --dt=$DT --sim_time=$TOTAL_TIME --index=$IND --basepath=$PWD --stepper=$STEPPER --Lx=$LX --Ly=$LY --meshx=$MESHX --meshy=$MESHY --cfl_safety=$CFL_SAFETY --cfl_threshold=$CFL_THRESHOLD --cfl_cadence=$CFL_CADENCE --cfl --snapshots


# For deciding the restart path
RECENT=$(find state/. -maxdepth 1 -type f -exec basename {} \; | sort -V | tail -n 1)

if [ ! -d "restart" ]; then
  mkdir restart
  echo "creating directory for restart"
fi

rm -rf restart/restart.h5
ln -sv $PWD/state/$RECENT $PWD/restart/restart.h5

python3 $PATH_TO_SCRIPTS/analysis_v3.py $PWD --time=$AVG_TIME

mkdir res_check
mkdir res_check_3d

mpirun python3 $PATH_TO_SCRIPTS/power_v3.py $PWD/state/*.h5 --ymin=$YMIN --ymax=$YMAX
mpirun python3 $SCRIPTS_3D/power_v3.py $PWD/state/*.h5 --ymin=$YMIN_3D --ymax=$YMAX_3D
ffmpeg -y -r 60 -pattern_type glob -i 'res_check/*.png' -threads 32 -pix_fmt yuv420p res_check/movie.mp4
ffmpeg -y -r 60 -pattern_type glob -i 'res_check_3d/*.png' -threads 32 -pix_fmt yuv420p res_check_3d/movie.mp4
