#!/bin/bash
# Create (or recreate) both virtual environments from the Alliance wheelhouse and
# install the rbc package into each. Run once per cluster, on a login node:
#     cd $SCRATCH/rbc_scripts && bash shell/install_envs.sh
# Uses ENV_PATH / VIS_ENV_PATH from the environment if set, else $SCRATCH/ded3(_vis).
set -euo pipefail
SCRIPTS_PATH=${SCRIPTS_PATH:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
ENV_PATH=${ENV_PATH:-$SCRATCH/ded3}
VIS_ENV_PATH=${VIS_ENV_PATH:-$SCRATCH/ded3_vis}
source "$SCRIPTS_PATH/shell/setup.sh"

make_env() {   # make_env PATH MODULES REQUIREMENTS EXTRA
    local path=$1 modules=$2 req=$3 extra=$4
    module purge
    module load $modules
    if [[ -d $path ]]; then
        echo "removing old environment $path"
        rm -rf -- "$path"
    fi
    virtualenv --no-download "$path"
    source "$path/bin/activate"
    pip install --no-index --upgrade pip
    pip install --no-index -r "$req"
    pip install --no-index -e "$SCRIPTS_PATH$extra"
    python -c "import dedalus, rbc; print('ok:', dedalus.__version__, rbc.__file__)"
    deactivate
}

make_env "$ENV_PATH" "$MODULES_SIM" "$SCRIPTS_PATH/requirements.txt" ""
make_env "$VIS_ENV_PATH" "$MODULES_VIS" "$SCRIPTS_PATH/requirements_vis.txt" "[vis]"
echo "done: $ENV_PATH and $VIS_ENV_PATH"
