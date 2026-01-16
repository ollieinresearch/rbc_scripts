#!/bin/bash


load_modules() {

    # Load the required modules
    ml StdEnv/2023
    ml python/3.11.5 mpi4py/4.0.3 fftw-mpi/3.3.10 scipy-stack/2023b ffmpeg/7.1.1 hdf5-mpi/1.14.4

    source /scinet/vast/etc/vastpreload-openmpi.bash

    source $ENV_PATH/bin/activate;

    # Dedalus performance tip!
    export OMP_NUM_THREADS=1
    export NUMEXPR_MAX_THREADS=1
    export HDF5_USE_FILE_LOCKING='False'
    
}

load_modules_plotting() {

    # Load the required modules
    ml StdEnv/2023
    ml python/3.11.5 mpi4py/4.0.3 fftw-mpi/3.3.10 scipy-stack/2023b ffmpeg/7.1.1 vtk/9.3.0
    ml hdf5-mpi/1.14.4

    source /scinet/vast/etc/vastpreload-openmpi.bash

    source $VIS_ENV_PATH/bin/activate;

    # Dedalus performance tip!
    export OMP_NUM_THREADS=1
    export NUMEXPR_MAX_THREADS=1
    export HDF5_USE_FILE_LOCKING='False'
    
}

