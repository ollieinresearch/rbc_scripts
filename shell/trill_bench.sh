# Add to shell/workflow.sh.
# Times each line of $BENCH_CONF ("RES_Z NODES C MESHX MESHY", from rbc.plan --emit) for
# $BENCH_ITERS steps at fixed $DT from a random start, with no output. Results are rows
# of $BENCH_DIR/bench.csv. Run inside one allocation as large as the largest line;
# smaller lines use the first nodes and leave the rest idle.
trill_bench() {
    local conf=${BENCH_CONF:-bench.conf} dir=${BENCH_DIR:-$PWD/bench}
    local sim=${SIM_SCRIPT:-$SCRIPTS_PATH/sim/rayleigh_benard.py}
    local res nodes c mx my rpn n run rc
    mkdir -p "$dir"
    export HDF5_USE_FILE_LOCKING=FALSE RBC_WATCHDOG=${BENCH_WATCHDOG:-600}
    while read -r res nodes c mx my <&3 || [[ -n $res ]]; do
        [[ -z $res || $res == \#* ]] && continue
        rpn=$(( 24 * (8 / c) ))
        n=$(( mx * my ))
        if (( n != nodes * rpn )); then
            echo "bench: skipping '$res $nodes $c $mx $my': ${mx}x${my} = $n ranks, $nodes nodes give $(( nodes * rpn ))"
            continue
        fi
        if (( nodes > SLURM_JOB_NUM_NODES )); then
            echo "bench: skipping '$res $nodes $c $mx $my': needs $nodes nodes, job has $SLURM_JOB_NUM_NODES"
            continue
        fi
        run="$dir/${res}_${nodes}n_c${c}_${mx}x${my}"
        mkdir -p "$run"
        export OMP_NUM_THREADS=$c
        echo "=== $(date +%T) bench RES_Z=$res, $nodes nodes, C=$c, mesh ${mx}x${my}"
        timeout "${BENCH_TIMEOUT:-1200}" mpirun -n "$n" --bind-to core \
            --map-by "ppr:$(( 8 / c )):l3cache:pe=$c" \
            -x HDF5_USE_FILE_LOCKING -x OMP_NUM_THREADS -x RBC_WATCHDOG \
            python3 "$sim" --Ra "$RA" --Pr "$PR" --Lx "$LX" --Ly "$LY" --nz "$res" \
            --meshx "$mx" --meshy "$my" --stepper "$STEPPER" --dt "$DT" --sim_time 1e9 \
            --basepath "$run" --bench "${BENCH_ITERS:-100}" --bench_csv "$dir/bench.csv" < /dev/null
        rc=$?
        (( rc == 0 )) || echo "bench: RES_Z=$res ${mx}x${my} on $nodes nodes failed (exit $rc; 124 = timeout)"
    done 3< "$conf"
    [[ -f $dir/bench.csv ]] && { column -s, -t "$dir/bench.csv" 2>/dev/null || cat "$dir/bench.csv"; }
}
