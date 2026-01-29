#!/bin/bash


load_modules() {

    # Load the required modules
    ml StdEnv/2023
    ml python/3.11.5 mpi4py/4.0.3 fftw-mpi/3.3.10 scipy-stack/2023b ffmpeg/7.1.1 hdf5-mpi/1.14.4
    ml
    source /scinet/vast/etc/vastpreload-openmpi.bash
    which mpirun

    source $ENV_PATH/bin/activate;

    # Dedalus performance tip!
    export OMP_NUM_THREADS=1
    export NUMEXPR_MAX_THREADS=1
    
}



load_modules_test() {

    # Load the required modules
    ml StdEnv/2023
    ml python/3.11.5 mpi4py/4.0.3 fftw-mpi/3.3.10 scipy-stack/2023b ffmpeg/7.1.1 hdf5-mpi/1.14.4
    ml

    source $ENV_PATH/bin/activate;

    # Dedalus performance tip!
    export OMP_NUM_THREADS=1
    export NUMEXPR_MAX_THREADS=1
    
}


load_modules_plotting() {

    # Load the required modules
    ml StdEnv/2023
    ml python/3.11.5 mpi4py/4.0.3 fftw-mpi/3.3.10 scipy-stack/2023b ffmpeg/7.1.1 vtk/9.3.0
    ml

    source /scinet/vast/etc/vastpreload-openmpi.bash

    source $VIS_ENV_PATH/bin/activate;

    # Dedalus performance tip!
    export OMP_NUM_THREADS=1
    export NUMEXPR_MAX_THREADS=1
    
}

create_restart() {
    rm -rf $PWD/state/*.lock
    # For deciding the restart path
    RECENT=$(find state/. -maxdepth 1 -type f -exec basename {} \; | sort -V | tail -n 1)

    if [ ! -d "restart" ]; then
    mkdir restart
    echo "creating directory for restart"
    fi

    # Set the new restart path
    rm -rf restart/restart.h5
    ln -sv $PWD/state/$RECENT $PWD/restart/restart.h5

}

