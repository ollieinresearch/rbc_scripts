#!/bin/bash


print_params() {
    echo "Start time: $START_TIME"
    echo "Rayleigh number:$RA"
    echo "Prandtl number: $PR"
    echo "Vertical resolution: $RES"
    echo "Factor for horizontal resolution: $GAM"
    echo "Initial timestep: $DT"
    echo "Length of simulation (dimensionless units): $SIM_TIME"
    echo "Timestepper: $STEPPER"
    echo "Aspect ratio for (x,y): ($LX, $LY)"
    echo "Mesh shape (x, y): ($MESHX, $MESHY)"
    echo "Start from a previous run? $IC"
    echo "Safety factor for CFL: $CFL_SAFETY"
    echo "Threshold for CFL: $CFL_THRESHOLD"
    echo "Cadence for CFL: $CFL_CADENCE"
    echo "Touch tmp_file? $TMP"
    echo "Parallel filehandler mode: $PARA"
    echo "Use CFL? $CFL"
    echo "Save snapshots? $SNAPSHOTS"

}



determine_ic() {

    # Specify desired time for initial condition - 0 implies most recent
    if [ $IC -eq 1 ]; then
        IC_ARRAY=($(python3 $SCRIPTS_PATH/initial_condition.py $START_TIME --file=$PWD/restart/restart.h5))
        TOTAL_TIME=$(echo "$SIM_TIME+${IC_ARRAY[1]}" | bc)
        IND=${IC_ARRAY[0]}
    else
        echo "Not running with a prior simulation's initial conditions; starting fresh!"
        # rm -rf {analysis,state,preliminary_outputs,snapshots,outputs,field_analysis,restart}
        rm -rf restart
        TOTAL_TIME=$SIM_TIME
        IND=-1
    fi

}



run_sim() {

    mpirun python3 "$SCRIPTS_3D/rayleigh_benard_script.py" \
    --Ra="$RA" \
    --Pr="$PR" \
    --nz="$RES" \
    --gamma=$GAM \
    --dt="$DT" \
    --sim_time="$TOTAL_TIME" \
    --index="$IND" \
    --basepath="$PWD" \
    --stepper="$STEPPER" \
    --Lx="$LX" \
    --Ly="$LY" \
    --meshx="$MESHX" \
    --meshy="$MESHY" \
    --cfl_safety="$CFL_SAFETY" \
    --cfl_threshold="$CFL_THRESHOLD" \
    --cfl_cadence="$CFL_CADENCE" \
    ${CFL:+--cfl} \
    ${SNAPSHOTS:+--snapshots}



}


srun_sim() {

    srun python3 "$SCRIPTS_3D/rayleigh_benard_script.py" \
    --Ra="$RA" \
    --Pr="$PR" \
    --nz="$RES" \
    --gamma=$GAM \
    --dt="$DT" \
    --sim_time="$TOTAL_TIME" \
    --index="$IND" \
    --basepath="$PWD" \
    --stepper="$STEPPER" \
    --Lx="$LX" \
    --Ly="$LY" \
    --meshx="$MESHX" \
    --meshy="$MESHY" \
    --cfl_safety="$CFL_SAFETY" \
    --cfl_threshold="$CFL_THRESHOLD" \
    --cfl_cadence="$CFL_CADENCE" \
    ${CFL:+--cfl} \
    ${SNAPSHOTS:+--snapshots}



}


test_parallel() {
    mpirun --timeout 300 python3 "$SCRIPTS_3D/rayleigh_benard_script.py" \
                --Ra="$RA" \
                --Pr="$PR" \
                --nz="$RES" \
                --gamma=$GAM \
                --dt="$DT" \
                --sim_time="$TOTAL_TIME" \
                --index="$IND" \
                --basepath="$PWD" \
                --stepper="$STEPPER" \
                --Lx="$LX" \
                --Ly="$LY" \
                --meshx="$MESHX" \
                --meshy="$MESHY" \
                --cfl_safety="$CFL_SAFETY" \
                --cfl_threshold="$CFL_THRESHOLD" \
                --cfl_cadence="$CFL_CADENCE" \
                ${CFL:+--cfl} \
                ${SNAPSHOTS:+--snapshots}

            
    for PARA in virtual mpio gather; do 
        for TMP in False True; do
            echo $PARA $TMP
            echo $PARA $TMP
            echo $PARA $TMP
            mpirun --timeout 1200 python3 "$SCRIPTS_3D/rayleigh_benard_script.py" \
                --Ra="$RA" \
                --Pr="$PR" \
                --nz="$RES" \
                --gamma=$GAM \
                --dt="$DT" \
                --sim_time="$TOTAL_TIME" \
                --index="$IND" \
                --basepath="$PWD" \
                --stepper="$STEPPER" \
                --Lx="$LX" \
                --Ly="$LY" \
                --meshx="$MESHX" \
                --meshy="$MESHY" \
                --cfl_safety="$CFL_SAFETY" \
                --cfl_threshold="$CFL_THRESHOLD" \
                --cfl_cadence="$CFL_CADENCE" \
                ${CFL:+--cfl} \
                ${SNAPSHOTS:+--snapshots}
            echo $PARA $TMP
            echo $PARA $TMP
            echo $PARA $TMP
            srun --time 20 python3 "$SCRIPTS_3D/rayleigh_benard_script.py" \
                --Ra="$RA" \
                --Pr="$PR" \
                --nz="$RES" \
                --gamma=$GAM \
                --dt="$DT" \
                --sim_time="$TOTAL_TIME" \
                --index="$IND" \
                --basepath="$PWD" \
                --stepper="$STEPPER" \
                --Lx="$LX" \
                --Ly="$LY" \
                --meshx="$MESHX" \
                --meshy="$MESHY" \
                --cfl_safety="$CFL_SAFETY" \
                --cfl_threshold="$CFL_THRESHOLD" \
                --cfl_cadence="$CFL_CADENCE" \
                ${CFL:+--cfl} \
                ${SNAPSHOTS:+--snapshots}
        done
    done



}


post_process() {

    # For deciding the restart path
    RECENT=$(find state/. -maxdepth 1 -type f -exec basename {} \; | sort -V | tail -n 1)

    if [ ! -d "restart" ]; then
    mkdir restart
    echo "creating directory for restart"
    fi

    # Set the new restart path
    rm -rf restart/restart.h5
    ln -sv $PWD/state/$RECENT $PWD/restart/restart.h5

    python3 $SCRIPTS_PATH/analysis_v3.py $PWD --time=$AVG_TIME
    python3 $SCRIPTS_PATH/h_analysis.py $PWD --time=$AVG_TIME
    
}

res_check() {
    #mkdir res_check
    #srun -c 3 python3 $PATH_TO_SCRIPTS/power_v3.py $PWD/state/*.h5 --ymin=$YMIN --ymax=$YMAX
    #ffmpeg -y -r 30 -pattern_type glob -i 'res_check/*.png' -threads 32 -pix_fmt yuv420p res_check/movie.mp4

    mkdir $PWD/res_check_3d
    export OMP_NUM_THREADS=3
    export NUMEXPR_MAX_THREADS=3
    srun -n 32 --cpus-per-task=6 python3 $SCRIPTS_3D/power_v3.py $PWD/state/*.h5 --mins=$MINS --maxs=$MAXS
    ffmpeg -y -r 30 -pattern_type glob -i 'res_check_3d/*.png' -threads 32 -pix_fmt yuv420p $PWD/res_check_3d/movie.mp4
  
}



plot_snapshots() {

    nu=$(grep "Final Nusselt number:" outputs/info.txt | awk -F': ' '{print $2}')

    mkdir $PWD/visualization

    srun -n 32 --cpus-per-task=6 python3 $SCRIPTS_3D/plotting_v3.py $PWD/snapshots/*.h5 --basepath=$PWD --nu=$nu
    ffmpeg -y -r 30 -pattern_type glob -i 'visualization/*.png' -threads 32 -pix_fmt yuv420p visualization/movie.mp4

}