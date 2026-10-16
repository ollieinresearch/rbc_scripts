# Installing on a cluster

Do this once on each cluster you use. Everything goes on {term}`scratch`, where jobs read
and write fastest.

## 1. Get the repository

```bash
cd $SCRATCH
git clone https://github.com/OWNER/rbc_scripts.git    # the address on the repository's GitHub page
```

`$SCRATCH/rbc_scripts` is the location every job script expects (`SCRIPTS_PATH` in
{doc}`../user-guide/parameters`).

## 2. Create the two environments

```bash
cd $SCRATCH/rbc_scripts
bash shell/install_envs.sh
```

This takes a few minutes and creates two Python virtual environments, each from the
Alliance {term}`wheelhouse` (`pip install --no-index`, no internet downloads):

| Environment | Path | Used for | Modules |
|---|---|---|---|
| `ded3` | `$SCRATCH/ded3` (`ENV_PATH`) | simulations, analysis, figures | `MODULES_SIM` in `shell/setup.sh` (includes `hdf5-mpi`) |
| `ded3_vis` | `$SCRATCH/ded3_vis` (`VIS_ENV_PATH`) | 3D rendering with pyvista | `MODULES_VIS` (includes `vtk`, not `hdf5-mpi`) |

There are two because pyvista needs the `vtk` module, whose dependencies differ slightly
from the simulation stack. Job scripts switch between them with `load_modules` and
`load_modules_plotting` ({doc}`../user-guide/workflow`).

Each environment gets the packages in `requirements.txt` (or `requirements_vis.txt`) and
then the `rbc` package itself, installed in *editable* mode: Python imports it straight
from `$SCRATCH/rbc_scripts`, so `git pull` updates it without reinstalling.

```{dropdown} What is in requirements.txt, and why
`requirements.txt` is Dedalus 3.0.4 with its own declared dependencies (`docopt`,
`xarray`, `pytest` and its plugins, which Dedalus uses for its self-tests), plus `pandas`
for the analysis. NumPy, SciPy, Matplotlib, mpi4py and h5py come from the `scipy-stack`,
`mpi4py` and `hdf5-mpi` modules, built by the Alliance against the cluster's MPI and HDF5,
which is why they are not pip-installed. Versions are pinned so every cluster runs the
same code.
```

## 3. Check

```bash
source $SCRATCH/ded3/bin/activate
python -c "import dedalus, rbc; print(dedalus.__version__, rbc.__version__)"
python -m rbc.collect --help
```

To check that Dedalus itself works, run its tests on a login node (a few minutes):
`python -m dedalus test`.

## 4. Optional: your laptop

The analysis and plotting run anywhere with Python 3.10 or newer; only the simulation
needs Dedalus and MPI.

```bash
git clone https://github.com/OWNER/rbc_scripts.git && cd rbc_scripts
python -m venv .venv && source .venv/bin/activate
pip install -e .            # add [vis] for 3D rendering, [docs] to build these pages
```

Copy `data/results.csv` from the cluster (it is not in git) and every
`python -m rbc.plots.*` command works locally, with LaTeX labels if LaTeX is installed.

## Updating

```bash
cd $SCRATCH/rbc_scripts && git pull
```

Rerun `install_envs.sh` only when `requirements*.txt` or the module versions in
`shell/setup.sh` change. It deletes and recreates both environments.

Next: {doc}`first-run`.
