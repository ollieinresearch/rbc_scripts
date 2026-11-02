#!/bin/bash
# Benchmark the configurations in bench.conf (see trill_bench in workflow.sh).
# Run from a folder of its own holding params.conf (for RA, PR, DT, STEPPER, LX, LY and
# the paths) and bench.conf. Set --nodes to the largest NODES in bench.conf.
#SBATCH --output=bench_%j.out
#SBATCH --time=01:30:00
#SBATCH --nodes=9
#SBATCH --ntasks-per-node=192
#SBATCH --cpus-per-task=1
#SBATCH --job-name=bench
#SBATCH --account=rrg-goluskin
#SBATCH --mem=0

source params.conf
source "$SETUP_PATH"
source "$WORKFLOW_PATH"

export OMP_PLACES=cores OMP_PROC_BIND=true
load_modules
trill_bench