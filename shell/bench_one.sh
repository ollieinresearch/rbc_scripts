#!/bin/bash
# One benchmark configuration on exactly the nodes it needs. Submitted by
# "python -m rbc.bench", which sets --nodes, --ntasks-per-node, --cpus-per-task, --time
# and --output. Runs in the benchmark folder (params.conf and bench.conf).
# Arguments: RES_Z NODES C MESHX MESHY RESULTS_DIR [RES_HOR, default 2 RES_Z]
#
# 1. Rehearsal launch: build everything and take BENCH_REHEARSAL steps (warms the nodes).
# 2. Measured launch, a short production job: the production outputs (params.conf:
#    SNAPSHOTS, A_FREQ, ...) are built before the first step, then BENCH_SECONDS of
#    stepping (default 600). The first BENCH_DISCARD seconds (default 60) are not timed.
#    Memory is read on every node throughout (bench_mem.csv) and at the end, after a
#    final write of every output. No mpirun timeout: stepping ends on its own, and ends
#    early enough to write the results if the Slurm time limit would otherwise cut it off.
# 3. Log what step 2 wrote (to show the checkpoint and outputs really were written),
#    then delete it (it is only there for the memory measurement).
#SBATCH --mem=0

res=$1 nodes=$2 c=$3 mx=$4 my=$5 dir=$6 nx=${7:-$(( 2 * $1 ))}
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
if (( nx == 2 * res )); then name=$res; else name="${res}_${nx}"; fi   # as rbc.bench.tag
run="$dir/${name}_${nodes}n_c${c}_${mx}x${my}"
mkdir -p "$run"

args=(--Ra "$RA" --Pr "$PR" --Lx "$LX" --Ly "$LY" --nz "$res" --nx "$nx"
      --meshx "$mx" --meshy "$my" --stepper "$STEPPER" --dt "$DT" --sim_time 1e9
      --basepath "$run" --a_freq "$A_FREQ" --snaps_freq "$SNAPS_FREQ" --spec_freq "$SPEC_FREQ"
      --state_freq 1e9 --state_layout "${STATE_LAYOUT:-c}" --par "${PARA:-virtual}"
      --par_small "${PAR_SMALL:-virtual}")
[[ ${SNAPSHOTS:-0} == 1 ]] && args+=(--snapshots)

# End of this allocation as Unix time, so the measured run can stop in time to save its results.
# 0 (no deadline) if Slurm does not say or the answer is not in the future.
end=$(squeue -h -j "$SLURM_JOB_ID" -o %e 2>/dev/null)
deadline=0
if [[ $end =~ ^[0-9]{4}- ]]; then
    deadline=$(date -d "$end" +%s 2>/dev/null || echo 0)
    (( deadline > $(date +%s) )) || deadline=0
fi
echo "=== allocation ends ${end:-unknown} (deadline $deadline)"

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
launch --bench "${BENCH_MAX_ITERS:-10000000}" --bench_seconds "${BENCH_SECONDS:-600}" \
    --bench_discard "${BENCH_DISCARD:-60}" --bench_outputs --bench_deadline "$deadline" \
    --bench_csv "$run/bench.csv"
rc=$?
echo "=== outputs written by the measured run (deleted next):"
for d in state analysis horizontal_analysis snapshots spectra; do
    [[ -d $run/$d ]] && echo "  $d: $(find "$run/$d" -name '*.h5' | wc -l) .h5 files, $(du -sh "$run/$d" | cut -f1)"
done
rm -rf "$run"/{state,snapshots,analysis,horizontal_analysis,spectra}
echo "=== $(date +%T) done (exit $rc)"
exit $rc
