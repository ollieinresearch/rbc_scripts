#!/bin/bash
# Module environments, sourced by every job script. The two venvs are created once by
# shell/install_envs.sh; both get the rbc package installed in editable mode.
#
# Two environments because pyvista needs the vtk module, and the plotting stack it pulls
# in differs slightly from the simulation one (hdf5-mpi is not loaded with vtk).

MODULES_SIM="StdEnv/2023 python/3.11.5 mpi4py/4.0.3 fftw-mpi/3.3.10 scipy-stack/2023b ffmpeg/7.1.1 hdf5-mpi/1.14.4"
MODULES_VIS="StdEnv/2023 python/3.11.5 mpi4py/4.0.3 fftw-mpi/3.3.10 scipy-stack/2023b ffmpeg/7.1.1 vtk/9.3.0"

# Trillium only: SciNet's VAST preload for Open MPI (absent elsewhere, so skipped).
_vast_preload() {
    [[ -f /scinet/vast/etc/vastpreload-openmpi.bash ]] && source /scinet/vast/etc/vastpreload-openmpi.bash
}

# Reinstall the rbc package (editable, no dependencies) only if it cannot be imported. An
# editable install records the repository's absolute path, so moving or renaming the
# folder breaks it while pip list still shows it.
_check_rbc() {
    python3 -c "import rbc" 2>/dev/null && return
    echo "rbc not importable in $(which python3); reinstalling from ${SCRIPTS_PATH:?}"
    pip install --no-index --no-deps -q -e "$SCRIPTS_PATH" \
        && python3 -c "import rbc" 2>/dev/null \
        || echo "WARNING: rbc package not installed in $(which python3)"
}

load_modules() {
    module load $MODULES_SIM
    _vast_preload
    source "$ENV_PATH/bin/activate"
    export OMP_NUM_THREADS=${OMP_NUM_THREADS:-$C}
    export NUMEXPR_MAX_THREADS=$OMP_NUM_THREADS
    _check_rbc
}

load_modules_trill_test() {
    module load StdEnv/2023 python/3.13.2 mpi4py/4.1.0 fftw-mpi/3.3.10 scipy-stack/2026a ffmpeg/7.1.1 hdf5-mpi/1.14.6
    _vast_preload
    source "$ENV_PATH/bin/activate"
    export OMP_NUM_THREADS=1 NUMEXPR_MAX_THREADS=1
    _check_rbc
}

load_modules_plotting() {
    module purge
    module load $MODULES_VIS
    _vast_preload
    source "$VIS_ENV_PATH/bin/activate"
    export OMP_NUM_THREADS=1 NUMEXPR_MAX_THREADS=1
    _check_rbc
}
