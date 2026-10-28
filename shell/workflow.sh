#!/bin/bash
# Functions sourced by run_sim.sh / post.sh. Expects params.conf and setup.sh sourced first.
# Python tools come from the rbc package: pip install -e $SCRIPTS_PATH into each venv.

SIM_SCRIPT="$SCRIPTS_PATH/sim/rayleigh_benard.py"

# Trillium's sbatch adds --export=NONE, which sets SLURM_EXPORT_ENV=NONE in the job. Every
# srun inherits it, including the one mpirun uses to start its daemons on the other nodes,
# so ranks off the batch node would start without the modules/venv/variables loaded here.
# (https://support.schedmd.com/show_bug.cgi?id=6772)
export SLURM_EXPORT_ENV=ALL
export RBC_WATCHDOG=${WATCHDOG:-1200} RBC_WATCHDOG_ALL=${WATCHDOG_ALL:-0}

print_params() {
    echo "Time: $(date)  Job: $SLURM_JOB_ID  Cluster: ${CC_CLUSTER:-?}"
    disk_report
    for v in RA PR RES_Z RES_HOR LX LY MESHX MESHY N C DT SIM_TIME STEPPER CFL CFL_SAFETY \
             CFL_THRESHOLD CFL_CADENCE IC START_TIME PARA PAR_SMALL STATE_LAYOUT STATE_FREQ \
             SNAPSHOTS STOP_MARGIN GRACEFUL_STOP VAST_PATCH WATCHDOG FIXED_GET_FILE; do
        printf '  %-14s %s\n' "$v" "${!v}"
    done
}

rb_args() {
    RB_ARGS=(
        --Ra="$RA" --Pr="$PR" --nz="$RES_Z" --nx="$RES_HOR"
        --dt="$DT" --sim_time="$TOTAL_TIME" --index="$IND" --basepath="$PWD"
        --stepper="$STEPPER" --Lx="$LX" --Ly="$LY" --meshx="$MESHX" --meshy="$MESHY"
        --cfl_safety="$CFL_SAFETY" --cfl_threshold="$CFL_THRESHOLD" --cfl_cadence="$CFL_CADENCE"
        --a_freq="$A_FREQ" --snaps_freq="$SNAPS_FREQ" --state_freq="$STATE_FREQ"
        --spec_freq="${SPEC_FREQ:-0.05}" --par="$PARA" --par_small="${PAR_SMALL:-gather}"
        --state_layout="${STATE_LAYOUT:-c}"
        ${CFL:+--cfl} ${SNAPSHOTS:+--snapshots} ${NO_SPEC:+--no_spectra}
    )
    [[ ${VAST_PATCH:-1} == 0 ]] && RB_ARGS+=(--no_vast_patch)
}

# Stop STOP_MARGIN s before the job's end so the final checkpoint fits (rbc.shutdown).
set_deadline() {
    export RBC_GRACEFUL=${GRACEFUL_STOP:-1}
    if [[ $RBC_GRACEFUL == 0 ]]; then
        unset RBC_DEADLINE
        echo "Graceful stop disabled (GRACEFUL_STOP=0)"
        return
    fi
    local end
    end=$(squeue -h -j "$SLURM_JOB_ID" -o %e 2>/dev/null)
    if [[ -z $end || $end == N/A || $end == Unknown ]]; then
        echo "set_deadline: no end time from squeue; relying on signals only"
        return
    fi
    # Margin = STOP_MARGIN, capped at 10% of the time left, so short test jobs still run.
    local end_s now left margin
    end_s=$(date -d "$end" +%s); now=$(date +%s); left=$(( end_s - now ))
    margin=${STOP_MARGIN:-900}
    (( margin > left / 10 )) && margin=$(( left / 10 ))
    export RBC_DEADLINE=$(( end_s - margin ))
    echo "Graceful stop planned at $(date -d @"$RBC_DEADLINE") (${margin} s before the job's end)"
}


# A run restarts from its own state/ once it has checkpoints; RESTART_FROM only seeds
# a new run folder (and may stay set afterwards).
determine_ic() {
    local src=${RESTART_FROM:-$PWD}
    if [[ $IC -eq 1 ]] && { compgen -G "$PWD/state/state_s*.h5" || compgen -G "$src/state/state_s*.h5"; } >/dev/null; then
        rm -f "$PWD"/state/*.lock
        read -r IND T_IC _ < <(python3 -m rbc.restart "$PWD" --from "$src" --time "$START_TIME" \
            --shape "$RES_HOR" "$RES_HOR" "$RES_Z")
        [[ -n $IND ]] || { echo "determine_ic: no usable restart (see the message above)"; exit 1; }
        TOTAL_TIME=$(echo "$SIM_TIME + $T_IC" | bc -l)
        echo "Restarting from write $IND at t = $T_IC; running to t = $TOTAL_TIME"
    else
        echo "Starting fresh"
        rm -rf restart
        TOTAL_TIME=$SIM_TIME
        IND=-1
    fi
}

# Lustre clusters (Rorqual, Fir): Slurm launch.
srun_sim() {
    rb_args
    set_deadline
    srun -n "$N" -c "$C" --cpu-bind=cores python3 "$SIM_SCRIPT" "${RB_ARGS[@]}"
}

# Trillium: Open MPI launch, ranks packed per L3 cache. The VAST fixes are now in
# rbc.vast (applied by the sim script), so fixed_get_file.py is no longer needed;
# set FIXED_GET_FILE=path to wrap with it anyway.
trill_run_sim() {
    rb_args
    set_deadline
    mpirun -n "$N" --bind-to core --map-by "ppr:${RANKS_PER_L3CACHE}:l3cache:pe=${OMP_NUM_THREADS}" \
        -x RBC_DEADLINE -x RBC_GRACEFUL -x RBC_WATCHDOG -x RBC_WATCHDOG_ALL \
        -x HDF5_USE_FILE_LOCKING -x OMP_NUM_THREADS \
        python3 ${FIXED_GET_FILE:+"$FIXED_GET_FILE"} "$SIM_SCRIPT" "${RB_ARGS[@]}"
}

# One-off check of the fast restart path at small scale (compares with Dedalus' loader).
verify_load() {
    rb_args
    srun -n "$N" -c "$C" python3 "$SIM_SCRIPT" "${RB_ARGS[@]}" --verify_load
}

# Trillium version of verify_load (the VAST preload needs mpirun).
trill_verify_load() {
    rb_args
    mpirun -n "$N" --bind-to core --map-by "ppr:${RANKS_PER_L3CACHE}:l3cache:pe=${OMP_NUM_THREADS}" \
        -x HDF5_USE_FILE_LOCKING -x OMP_NUM_THREADS \
        python3 "$SIM_SCRIPT" "${RB_ARGS[@]}" --verify_load
}

# Post-processing runs as its own 1-node job after the sim job ends (any exit state).
# Post-processing job, chained to this one. Preferred: queue both from a login node
# with shell/submit.sh (works on every cluster; required on Trillium, which refuses
# sbatch from compute nodes). submit.sh leaves .post_queued so this does not double-submit.
submit_post() {
    if [[ -f .post_queued ]]; then
        echo "post-processing already queued by submit.sh (job $(cat .post_queued))"
        rm -f .post_queued
        return
    fi
    if [[ ${CC_CLUSTER:-} == trillium ]] || ! sbatch --dependency=afterany:"$SLURM_JOB_ID" \
            --account="$SLURM_JOB_ACCOUNT" --job-name="post_${SLURM_JOB_NAME}" "$SCRIPTS_PATH/shell/post.sh"; then
        echo "could not queue post-processing from here; on a login node run:"
        echo "  cd $PWD && sbatch $SCRIPTS_PATH/shell/post.sh"
    fi
}

# AVG_TIME may be a number or "auto" (MSER-5); profiles and spectra then use the window
# start that rbc.analysis.nusselt actually chose.
post_process() {
    python3 -m rbc.analysis.nusselt "$PWD" --t0 "${AVG_TIME:-0}"
    local t0
    t0=$(python3 -c "import json; print(json.load(open('outputs/nusselt.json'))['t0'])") || t0=${AVG_TIME:-0}
    python3 -m rbc.analysis.profiles "$PWD" --t0 "$t0"
    [[ -d spectra ]] && python3 -m rbc.spectra.cumulative "$PWD" --time="$t0"
    python3 -m rbc.spectra.compare "$(dirname -- "$PWD")" || true
}

# 3D render movie with the spectra panel underneath (needs SNAPSHOTS and the vis env).
plot_snapshots() {
    python3 -m rbc.vis.render3d "$PWD" --mode "${VIS_MODE:-iso}" --spectra --movie \
        --jobs "${VIS_JOBS:-0}" --stride "${VIS_STRIDE:-1}" --fps "$FPS"
}

# Movie of horizontal slices of this run (needs SNAPSHOTS and the vis env).
plot_slices() {
    python3 -m rbc.vis.slices "$PWD" --fields T w --heights mid bl --movie \
        --jobs "${VIS_JOBS:-0}" --fps "${FPS:-30}" --out "$PWD/visualization"
}

# Slices of any number of runs side by side, as a movie: plot_compare_slices RUN_A RUN_B RUN_C ...
# (without arguments, the Ra = 1e6 pair at Pr = 1 and Pr = 10^-3.5).
plot_compare_slices() {
    local runs=("$@")
    (( ${#runs[@]} )) || runs=("$SCRATCH/ollie_rb_data/3d/ra1e6/pr1e0/64_Gam2_2"
                                 "$SCRATCH/ollie_rb_data/3d/ra1e6/pr3e-4/216_Gam2_2")
    python3 -m rbc.vis.slices "${runs[@]}" --fields T w --heights mid bl --movie --jobs "${VIS_JOBS:-0}"
}

# Spectra computed from saved snapshots (runs made before in-sim spectra existed).
snapshot_spectra() {
    srun -n "${SPEC_TASKS:-32}" -c 1 python3 -m rbc.spectra.snapshots snapshots/snapshots_s*.h5 --plot-every=0
}

# Mirror every run under DATA_ROOT into each of BACKUP_DIRS (same path below $SCRATCH):
# analysis, spectra, outputs, logs and the BACKUP_STATES newest readable checkpoints.
# Older checkpoints are then removed from each backup, and on scratch everything but
# the KEEP_STATES newest, but only if every copy succeeded. Snapshots are not copied.
backup_and_prune() {
    local src=${DATA_ROOT:-$SCRATCH/ollie_rb_data/3d} rel_root dst ok=1
    rel_root=${src#"$SCRATCH"/}
    for dst in ${BACKUP_DIRS:-$HOME}; do
        dst="$dst/$rel_root"
        mkdir -p "$dst" || { ok=0; continue; }
        rsync -a -m --include='*/' \
            --include='analysis/***' --include='horizontal_analysis/***' --include='spectra/***' \
            --include='outputs/***' --include='params.conf' --include='stop_reason.txt' \
            --include='job_*.out' --include='run_sim*.sh' --exclude='*' \
            "$src/" "$dst/" || { echo "backup to $dst failed"; ok=0; continue; }
        _backup_states "$src" "$dst" || ok=0
    done
    if [[ $ok -ne 1 ]]; then
        echo "a backup failed; nothing pruned on scratch"
        return 1
    fi
    local run
    while IFS= read -r -d '' run; do
        python3 -m rbc.restart "${run%/state}" --keep "${KEEP_STATES:-4}" --prune >/dev/null
    done < <(find "$src" -type d -name state -prune -print0)
}

# Copy the BACKUP_STATES newest readable checkpoints of each run from $1 to $2,
# then delete older checkpoints in $2.
_backup_states() {
    local src=$1 dst=$2 run rel p ok=1
    while IFS= read -r -d '' run; do
        run=${run%/state}
        rel=${run#"$src"/}
        mapfile -t keep < <(python3 -m rbc.restart "$run" --keep "${BACKUP_STATES:-2}")
        [[ ${#keep[@]} -gt 0 ]] || continue
        mkdir -p "$dst/$rel/state"
        for p in "${keep[@]}"; do
            [[ -e $p ]] && { rsync -a "$p" "$dst/$rel/state/" || ok=0; }
        done
        [[ $ok -eq 1 ]] && python3 -m rbc.restart "$dst/$rel" --keep "${BACKUP_STATES:-2}" --prune >/dev/null
    done < <(find "$src" -type d -name state -prune -print0)
    [[ $ok -eq 1 ]]
}

# results.csv (and results.parquet) -> data/ in the repo, then a copy next to the run
# tree in each backup. Both are git-ignored until the paper is out.
collect_results() {
    local src=${DATA_ROOT:-$SCRATCH/ollie_rb_data/3d} dst f
    python3 -m rbc.collect "$src" --no-plots --jobs "${SLURM_CPUS_PER_TASK:-8}" || return 1
    for dst in ${BACKUP_DIRS:-$HOME}; do
        dst="$dst/${src#"$SCRATCH"/}"
        mkdir -p "$dst"
        for f in "$SCRIPTS_PATH"/data/results.csv "$SCRIPTS_PATH"/data/results.parquet; do
            [[ -f $f ]] && cp -p "$f" "$dst/"
        done
    done
}

# Quota summary at the top of every job log: a job that fills its quota dies mid-write
# with an unhelpful HDF5 error (e.g. about byte counts overflowing).
disk_report() {
    command -v diskusage_report >/dev/null && diskusage_report 2>/dev/null | head -n 20
}

# Sweeps file-handler modes and HDF5 locking, TEST_SECONDS per case.
test_parallel() {
    rb_args
    local lock para
    for lock in FALSE TRUE; do
        export HDF5_USE_FILE_LOCKING=$lock
        for para in virtual mpio gather; do
            reset_test
            echo "=== locking=$lock par=$para"
            timeout "${TEST_SECONDS:-1200}" srun python3 "$SIM_SCRIPT" "${RB_ARGS[@]}" --par="$para"
        done
    done
    export HDF5_USE_FILE_LOCKING=FALSE
}

reset_test() {
    rm -f state/*.loc state/*.lock
    local recent restart
    recent=$(ls state/state_s*.h5 2>/dev/null | sort -V | tail -n 1)
    restart=$(readlink -f restart/restart.h5)
    if [[ -n $recent && $(readlink -f "$recent") != "$restart" ]]; then
        echo "Removing test-written $recent"
        rm -rf "${recent%.h5}" "$recent"
    fi
}
