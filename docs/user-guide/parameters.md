# `params.conf` reference

Each run directory has its own `params.conf`, copied from `shell/params.conf` and edited.
It is a shell file: every line is `export NAME=value`, and lines starting with `#` are off.
Job scripts `source` it, so it may use shell arithmetic such as `$((RES_Z * 2))`.

The `#SBATCH` lines at the top of the job script (nodes, time, job name, account) are the
other half of a run's settings; they must agree with `N` and `C` below
({ref}`sbatch-lines`).

## Simulation

| Variable | Template | Meaning |
|---|---|---|
| `RES_Z` | 192 | Chebyshev modes $N_z$ in the vertical. |
| `RES_HOR` | `2 * RES_Z` | Fourier modes $N_x = N_y$; twice `RES_Z` keeps the same resolution per unit length in a Γ = 2 box. Must be even. |
| `MESHX`, `MESHY` | 24, 24 | {term}`Process mesh`. Each should divide `RES_HOR`. |
| `N` | `MESHX * MESHY` | Number of MPI ranks. |
| `C` | 1 | Cores (OpenMP threads) per rank. |
| `RA` | 7 | $\log_{10} Ra$. |
| `PR` | -2 | $\log_{10} Pr$, on half-decades. |
| `SIM_TIME` | 55 | Free-fall times to run: from 0 for a fresh start, added to the restart time otherwise. |
| `STEPPER` | `RK222` | Dedalus timestepper class ({ref}`timestepping`). |
| `CFL` | 1 | Adaptive timestep; comment out for a fixed `DT`. |
| `CFL_SAFETY` | 1/3 | Safety factor $C$ in $\Delta t \le C\,\Delta x/\lvert u\rvert$. |
| `CFL_THRESHOLD` | 0.1 | Relative change below which the timestep is not updated. |
| `CFL_CADENCE` | 2 | Iterations between CFL updates. |
| `DT` | 0.0005 | Fixed timestep, or the first one with `CFL` on a fresh start (a restart uses the checkpoint's). |
| `LX`, `LY` | 2 | Box size; the height is 1. |

## Output

| Variable | Template | Meaning |
|---|---|---|
| `A_FREQ` | 0.005 | Simulation time between `analysis` and `horizontal_analysis` writes. |
| `SPEC_FREQ` | 0.05 | Simulation time between spectra writes. |
| `SNAPSHOTS` | off | Set to 1 to write full 3D snapshots (large). |
| `SNAPS_FREQ` | 1/30 | Simulation time between snapshots. |
| `STATE_FREQ` | 3600 | Wall-clock seconds between checkpoints. A final one is always written when the job ends cleanly. |
| `STATE_LAYOUT` | `c` | Checkpoints in coefficient space (`c`: exact, and 3.4 = (3/2)³ times smaller than the dealiased grid) or grid space (`g`). |
| `PARA` | `virtual` | Dedalus file-handler mode for checkpoints and snapshots: `virtual`, `mpio` or `gather` ({ref}`handler-modes`). |
| `PAR_SMALL` | `virtual` | Mode for the small `analysis` and `horizontal_analysis` output. |
| `NO_SPEC` | off | Set to 1 to skip the in-simulation spectra. |

## Restarts

| Variable | Template | Meaning |
|---|---|---|
| `IC` | 1 | 1: restart from a checkpoint if there is one; 0: start fresh (from rest with small noise). |
| `START_TIME` | 0 | Restart from the checkpoint nearest this time; 0 means the newest. |
| `RESTART_FROM` | off | Another run directory whose checkpoint seeds this one, e.g. a lower resolution of the same point; it is regridded if the sizes differ. Used only while this run has no checkpoints of its own, so it can stay set ({doc}`restarts`). |

## Robustness

| Variable | Template | Meaning |
|---|---|---|
| `STOP_MARGIN` | 900 | Seconds before the job's time limit at which to stop and write the final checkpoint (at most a tenth of the job). |
| `GRACEFUL_STOP` | 1 | 0 turns off the deadline and signal handling ({ref}`graceful-stop`). |
| `WATCHDOG` | 1200 | Seconds without progress before the stacks are dumped to `watchdog/` and the job ends; 0 is off ({ref}`watchdog`). |
| `WATCHDOG_ALL` | off | 1 dumps every rank's stack, not one per node. |
| `VAST_PATCH` | on | 0 uses Dedalus' own file-handler code instead of {mod}`rbc.vast`. |
| `FIXED_GET_FILE` | off | Path of the older `fixed_get_file.py` wrapper, to run the simulation through it (needs `VAST_PATCH=0`). |

## Post-processing and visualization

| Variable | Template | Meaning |
|---|---|---|
| `AVG_TIME` | `auto` | Start of the averaging window, or `auto` for the MSER-5 estimate ({ref}`transient`). |
| `FPS` | 30 | Frame rate of movies. |
| `VIS_MODE` | `iso` | 3D render style: `iso` (isosurfaces) or `volume`. |
| `VIS_JOBS` | 0 | Frames rendered in parallel; 0 chooses from the cores and memory ({ref}`vis-workers`). |
| `VIS_STRIDE` | 1 | Horizontal subsampling for renders (2 uses a quarter of the memory). |

## Storage

| Variable | Template | Meaning |
|---|---|---|
| `DATA_ROOT` | `$SCRATCH/ollie_rb_data/3d` | Top of the run tree. Backed up, pruned and collected as a whole. |
| `BACKUP_DIRS` | `$HOME` and `$PROJECT` | Where the tree is mirrored (same path below `$SCRATCH`). Fir has no `$PROJECT`, so `/project/def-goluskin/$USER` is used there. |
| `BACKUP_STATES` | 2 | Newest checkpoints copied to each backup; older ones are removed from the backups. |
| `KEEP_STATES` | 4 | Newest checkpoints kept on scratch after a successful backup. |

## Paths

| Variable | Template | Meaning |
|---|---|---|
| `SCRIPTS_PATH` | `$SCRATCH/rbc_scripts` | The repository. |
| `SETUP_PATH`, `WORKFLOW_PATH` | in `shell/` | Sourced by every job script. |
| `ENV_PATH` | `$SCRATCH/ded3` | Simulation environment. |
| `VIS_ENV_PATH` | `$SCRATCH/ded3_vis` | Visualization environment. |

## Used without being in the template

| Variable | Default | Meaning |
|---|---|---|
| `TEST_SECONDS` | 1200 | Length of each case of `test_parallel`. |
| `SPEC_TASKS` | 32 | Ranks for `snapshot_spectra`. |
| `CORES_PER_L3CACHE` | 8 | Trillium: cores sharing an L3 cache; set in `run_sim_trillium.sh`. |
