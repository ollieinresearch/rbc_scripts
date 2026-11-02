#!/bin/bash
# One benchmark configuration on exactly the nodes it needs. Submitted by
# "python -m rbc.bench", which sets --nodes, --ntasks-per-node, --cpus-per-task, --time
# and --output. Runs in the benchmark folder (params.conf and bench.conf).
# Arguments: RES_Z NODES C MESHX MESHY RESULTS_DIR
#
# 1. Rehearsal launch: build everything and take BENCH_REHEARSAL steps (warms the nodes).
# 2. Measured launch: time BENCH_SECONDS of stepping, then build the production outputs
#    (params.conf: SNAPSHOTS, A_FREQ, ...) and take BENCH_IO steps to measure memory.
# 3. Delete the output files written in step 2 (they are only there for the memory).
#SBATCH --mem=0

res=$1 nodes=$2 c=$3 mx=$4 my=$5 dir=$6
if [[ -z $SLURM_JOB_ID ]]; then
    echo "bench_one.sh: not inside a Slurm job; run it from a debugjob or submit it with sbatch"
    exit 1
fi
source params.conf
[[ -f bench.conf ]] && source bench.conf
source "$SETUP_PATH"
source "$WORKFLOW_PATH"
load_modules

export OMP_NUM_THREADS=$c OMP_PLACES=cores OMP_PROC_BIND=true
export HDF5_USE_FILE_LOCKING=FALSE RBC_WATCHDOG=${BENCH_WATCHDOG:-900}
run="$dir/${res}_${nodes}n_c${c}_${mx}x${my}"
mkdir -p "$run"

args=(--Ra "$RA" --Pr "$PR" --Lx "$LX" --Ly "$LY" --nz "$res" --nx $(( 2 * res ))
      --meshx "$mx" --meshy "$my" --stepper "$STEPPER" --dt "$DT" --sim_time 1e9
      --basepath "$run" --a_freq "$A_FREQ" --snaps_freq "$SNAPS_FREQ" --spec_freq "$SPEC_FREQ"
      --state_freq 1e9 --state_layout "${STATE_LAYOUT:-c}" --par "${PARA:-virtual}"
      --par_small "${PAR_SMALL:-virtual}")
[[ ${SNAPSHOTS:-0} == 1 ]] && args+=(--snapshots)

launch() {
    mpirun -n $(( mx * my )) --bind-to core --map-by "ppr:$(( 8 / c )):l3cache:pe=$c" \
        -x HDF5_USE_FILE_LOCKING -x OMP_NUM_THREADS -x RBC_WATCHDOG \
        python3 "${SIM_SCRIPT:-$SCRIPTS_PATH/sim/rayleigh_benard.py}" "${args[@]}" "$@" < /dev/null
}

if (( ${BENCH_REHEARSAL:-10} > 0 )); then
    echo "=== $(date +%T) rehearsal (${BENCH_REHEARSAL:-10} steps, not recorded)"
    launch --bench_rehearsal "${BENCH_REHEARSAL:-10}" || echo "rehearsal exited with $?"
fi
echo "=== $(date +%T) measured run"
launch --bench "${BENCH_MAX_ITERS:-100000}" --bench_seconds "${BENCH_SECONDS:-240}" \
    --bench_io "${BENCH_IO:-30}" --bench_csv "$run/bench.csv"
rc=$?
rm -rf "$run"/{state,snapshots,analysis,horizontal_analysis,spectra}
echo "=== $(date +%T) done (exit $rc)"
exit $rc
