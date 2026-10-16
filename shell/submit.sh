#!/bin/bash
# Run from a login node inside the run directory (any cluster):
#     bash $SCRIPTS_PATH/shell/submit.sh                      # picks the launcher for this cluster
#     bash $SCRIPTS_PATH/shell/submit.sh run_sim.sh [sbatch options]
# Queues the simulation, then post.sh to start when it ends (whatever its exit state).
set -e
if [[ -n $1 && $1 != -* ]]; then script=$1; shift
elif [[ ${CC_CLUSTER:-} == trillium ]]; then script=run_sim_trillium.sh
else script=run_sim.sh; fi
[[ -f $script ]] || { echo "no $script in $PWD (copy it from \$SCRIPTS_PATH/shell/)"; exit 1; }
source params.conf
echo pending > .post_queued
jid=$(sbatch --parsable "$@" "$script")
if pid=$(sbatch --parsable --dependency=afterany:"$jid" --job-name="post_${jid}" \
         "$SCRIPTS_PATH/shell/post.sh"); then
    echo "$pid" > .post_queued
    echo "simulation job $jid ($script), post-processing job $pid (starts after $jid)"
else
    rm -f .post_queued
    echo "simulation job $jid queued; post-processing NOT queued"
fi
