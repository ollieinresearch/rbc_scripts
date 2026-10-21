#!/bin/bash
#SBATCH --output=jobs/job_post_%j.out
#SBATCH --time=02:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=192
#SBATCH --mem=0
#SBATCH --account=rrg-goluskin


# Submitted by run_sim.sh with --dependency=afterany; runs in the same run directory.
source params.conf
source "$SETUP_PATH"
source "$WORKFLOW_PATH"

load_modules
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK NUMEXPR_MAX_THREADS=$SLURM_CPUS_PER_TASK

post_process
if [[ -n $SNAPSHOTS ]]; then
    load_modules_plotting
    export OMP_NUM_THREADS=1
    plot_snapshots
    load_modules
fi
backup_and_prune
collect_results
