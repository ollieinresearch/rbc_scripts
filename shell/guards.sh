# Two additions to shell/workflow.sh and run_sim_trillium.sh.

# 1. First line inside trill_run_sim and srun_sim: refuse to start when the process mesh
#    does not fill the allocation exactly (an unsaved params.conf edit once left 16 of 48
#    nodes idle). SLURM_NTASKS may be unset when only --nodes and --ntasks-per-node are
#    given, hence the fallback.
check_ranks() {
    local per=${SLURM_NTASKS_PER_NODE%%(*}
    local slots=${SLURM_NTASKS:-$(( SLURM_JOB_NUM_NODES * per ))}
    if (( MESHX * MESHY != slots )); then
        echo "MESHX*MESHY = $(( MESHX * MESHY )) ranks but the job has $slots slots" \
             "($SLURM_JOB_NUM_NODES nodes); check params.conf and #SBATCH --nodes. Not starting."
        exit 1
    fi
}

# 2. In run_sim_trillium.sh (and run_sim.sh): make the job's Slurm state follow the
#    simulation. Without the final "exit", the script returns the status of its last
#    command (submit_post), so a run killed for memory still shows COMPLETED.
#
#   determine_ic
#   rm -f stop_reason.txt
#   trill_run_sim
#   rc=$?
#   echo "Simulation exited with $rc ($(cat stop_reason.txt 2>/dev/null || echo 'no stop_reason.txt'))"
#   submit_post
#   exit $rc
