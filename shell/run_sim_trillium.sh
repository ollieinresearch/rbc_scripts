#!/bin/bash
# Trillium (VAST): mpirun launch packed per L3 cache. Lustre clusters: use run_sim.sh.
# Submit from a login node with: bash $SCRIPTS_PATH/shell/submit.sh
#SBATCH --output=job_sim_%j.out
#SBATCH --time=11:59:00
#SBATCH --nodes=27
#SBATCH --ntasks-per-node=192
#SBATCH --cpus-per-task=1
#SBATCH --job-name=ra8pr-2.5
#SBATCH --account=rrg-goluskin

source params.conf
source "$SETUP_PATH"
source "$WORKFLOW_PATH"

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK OMP_PLACES=cores OMP_PROC_BIND=true
export CORES_PER_L3CACHE=8
export RANKS_PER_L3CACHE=$(( CORES_PER_L3CACHE / OMP_NUM_THREADS ))   # valid up to 8 threads
export HDF5_USE_FILE_LOCKING=FALSE

print_params
load_modules
determine_ic
trill_run_sim
echo "Simulation exited with $? ($(cat stop_reason.txt 2>/dev/null))"
submit_post
