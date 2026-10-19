#!/bin/bash
# Lustre clusters (Fir, Rorqual, Nibi): srun launch. Trillium: use run_sim_trillium.sh.
# Submit with: bash $SCRIPTS_PATH/shell/submit.sh   (queues post-processing too)
#SBATCH --output=jobs/job_sim_%j.out
#SBATCH --time=11:59:00
#SBATCH --nodes=27
#SBATCH --ntasks-per-node=192
#SBATCH --mem=0
#SBATCH --job-name=ra8pr-2.5
#SBATCH --account=rrg-goluskin

set -o pipefail
source params.conf
source "$SETUP_PATH"
source "$WORKFLOW_PATH"

print_params
load_modules
export HDF5_USE_FILE_LOCKING=FALSE

determine_ic
srun_sim
echo "Simulation step exited with $? ($(cat stop_reason.txt 2>/dev/null))"

submit_post
