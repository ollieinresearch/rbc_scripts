#!/bin/bash
#SBATCH --output=job_post.out
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --mem=0
#SBATCH --account=def-goluskin

# When to start averaging - post process
AVG_TIME=50
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

source $PATH_TO_ENV/bin/activate;

# Dedalus performance tip!
export OMP_NUM_THREADS=1
export NUMEXPR_MAX_THREADS=1

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
