# The workflow scripts

Everything a job does is a shell function in `shell/workflow.sh`, called by short job
scripts. This page follows a run through them; {doc}`parameters` lists the variables they
read.

## The files

| File | Role |
|---|---|
| `shell/sim_ready.sh` | Run in a new run directory: copies the job script for this cluster, `post.sh`, `submit.sh` and (if absent) `params.conf` into it. |
| `shell/submit.sh` | Run on a login node in the run directory. Makes `jobs/` for the logs, submits the job script for this cluster and queues `post.sh` after it. |
| `shell/run_sim.sh` | Job script for Fir, Rorqual and Nibi ({term}`Lustre`): launches with `srun`. |
| `shell/run_sim_trillium.sh` | Job script for Trillium ({term}`VAST`): launches with `mpirun`, ranks placed per {term}`L3 cache`. |
| `shell/post.sh` | One-node post-processing job. |
| `shell/workflow.sh` | The functions below. |
| `shell/setup.sh` | Module lists and the `load_modules*` functions. |
| `shell/params.conf` | Template of the per-run settings; copied into each run directory. |
| `shell/install_envs.sh` | Creates the two Python environments ({doc}`../getting-started/install`). |
| `sim/rayleigh_benard.py` | The Dedalus simulation ({doc}`../theory/equations`). |

Every job script starts the same way: `source params.conf`, then `setup.sh` and
`workflow.sh` from the paths it sets.

## A simulation job

`run_sim.sh` (or `run_sim_trillium.sh`) calls, in order:

`print_params`
: Writes the date, job, cluster, every setting and a disk-quota summary
  (`disk_report`, using the Alliance `diskusage_report` command) to the log.

`load_modules`
: Loads `MODULES_SIM`, the VAST preload on Trillium, and activates the simulation
  environment. Sets `OMP_NUM_THREADS` to `C`.

`determine_ic`
: Decides between a fresh start and a restart ({doc}`restarts`): with `IC=1` it restarts
  if this run, or the run named by `RESTART_FROM`, has a checkpoint. It runs
  {mod}`rbc.restart`, which takes this run's own newest readable checkpoint (or the one
  nearest `START_TIME`), or the other run's if this one has none yet, links it as
  `restart/restart.h5` (or regrids it there), and prints its index and time; the run then
  goes to that time plus `SIM_TIME`.

`set_deadline`
: Asks Slurm when the job ends and exports `RBC_DEADLINE`, `STOP_MARGIN` seconds earlier
  (at most a tenth of the time left), for the {ref}`graceful stop <graceful-stop>`.

`srun_sim` or `trill_run_sim`
: Builds the command line from `params.conf` (`rb_args`) and launches the simulation.
  The two launchers are compared in {doc}`clusters`.

`submit_post`
: Queues `post.sh` if `submit.sh` did not already (it leaves a `.post_queued` marker).
  Trillium refuses `sbatch` from compute nodes, so there it only prints the command.

The simulation writes `stop_reason.txt` (`finished`, `deadline`, `SIGUSR1`, ...) and a
row of `outputs/jobs.csv` when it ends.

## A post-processing job

`post.sh` runs on one full node after the simulation job, whatever its exit state:

`post_process`
: {mod}`rbc.analysis.nusselt` (Nu, Re, convergence), then {mod}`rbc.analysis.profiles`
  and {mod}`rbc.spectra.cumulative` from the averaging start that the first one chose,
  then {mod}`rbc.spectra.compare` over every resolution of the same (Ra, Pr).

`plot_snapshots`
: Only with `SNAPSHOTS=1`: switches to the visualization environment and renders a movie
  with {mod}`rbc.vis.render3d`.

`backup_and_prune`
: Mirrors every run under `DATA_ROOT` to each of `BACKUP_DIRS`, then deletes old
  checkpoints ({doc}`storage`).

`collect_results`
: Runs {mod}`rbc.collect` over `DATA_ROOT` and copies `results.csv` (and
  `results.parquet`) into the backups.

## Other functions

`verify_load`, `trill_verify_load`
: Loads the newest checkpoint with both the fast loader and Dedalus' own, and reports the
  largest difference (should be 0). Run once at small scale after changing
  {mod}`rbc.vast`. The Trillium version launches with `mpirun`, which the VAST preload
  needs.

`plot_slices`
: A movie of horizontal slices of this run (temperature and vertical velocity at mid-height
  and in the thermal boundary layer) into `visualization/`, with {mod}`rbc.vis.slices`.
  Needs `SNAPSHOTS=1` and the visualization environment (`load_modules_plotting`).

`plot_compare_slices RUN_A RUN_B RUN_C ...`
: The same for any number of runs side by side, one row each; without arguments, the
  $Ra = 10^6$ runs at $Pr = 1$ and $Pr = 10^{-3.5}$.

`test_parallel`, `reset_test`
: Run the simulation for `TEST_SECONDS` with each file-handler mode and HDF5 locking
  setting, deleting the test output between cases. For diagnosing filesystem problems.

`snapshot_spectra`
: Spectra of runs made before spectra were computed in the simulation
  ({mod}`rbc.spectra.snapshots`).

## Running pieces by hand

Every Python step is an ordinary command, so any of them can be rerun on a login node
or in an interactive job, for example after changing the averaging window:

```bash
source $SCRATCH/ded3/bin/activate
cd $RUN
python -m rbc.analysis.nusselt . --t0 120
python -m rbc.analysis.profiles . --t0 120
```

The full list is in {doc}`../reference/cli`.
