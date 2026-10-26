#!/bin/bash
# Trillium (VAST): mpirun launch packed per L3 cache. Lustre clusters: use run_sim.sh.
# Submit from a login node with: bash $SCRIPTS_PATH/shell/submit.sh
#SBATCH --output=jobs/job_sim_%j.out
#SBATCH --time=11:59:00
#SBATCH --nodes=27
#SBATCH --ntasks-per-node=192
#SBATCH --cpus-per-task=1
#SBATCH --job-name=ra8pr-2.5
#SBATCH --account=rrg-goluskin
#SBATCH --mem=0

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
rm -f stop_reason.txt
trill_run_sim
rc=$?
reason=$(cat stop_reason.txt 2>/dev/null) || reason="no stop_reason.txt: crashed or watchdog, see watchdog/"
echo "Simulation exited with $rc ($reason)"
submit_post
