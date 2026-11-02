#!/bin/bash
# One benchmark configuration on exactly the nodes it needs. Submitted by
# "python -m rbc.bench submit", which sets --nodes, --ntasks-per-node, --cpus-per-task,
# --time and --output. Run from a folder holding params.conf (RA, PR, DT, STEPPER, LX,
# LY and the paths). Arguments: RES_Z NODES C MESHX MESHY RESULTS_DIR
#SBATCH --mem=0

res=$1 nodes=$2 c=$3 mx=$4 my=$5 dir=$6
source params.conf
source "$SETUP_PATH"
source "$WORKFLOW_PATH"
load_modules

export OMP_NUM_THREADS=$c OMP_PLACES=cores OMP_PROC_BIND=true
export HDF5_USE_FILE_LOCKING=FALSE RBC_WATCHDOG=${BENCH_WATCHDOG:-600}
run="$dir/${res}_${nodes}n_c${c}_${mx}x${my}"
mkdir -p "$run"

mpirun -n $(( mx * my )) --bind-to core --map-by "ppr:$(( 8 / c )):l3cache:pe=$c" \
    -x HDF5_USE_FILE_LOCKING -x OMP_NUM_THREADS -x RBC_WATCHDOG \
    python3 "${SIM_SCRIPT:-$SCRIPTS_PATH/sim/rayleigh_benard.py}" \
    --Ra "$RA" --Pr "$PR" --Lx "$LX" --Ly "$LY" --nz "$res" --meshx "$mx" --meshy "$my" \
    --stepper "$STEPPER" --dt "$DT" --sim_time 1e9 --basepath "$run" \
    --bench "${BENCH_ITERS:-100}" --bench_csv "$run/bench.csv" < /dev/null
exit $?
