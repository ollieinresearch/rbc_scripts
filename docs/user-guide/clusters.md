---
html_theme.sidebar_secondary.remove: true
---

# Fir and Trillium

The same repository, environments and commands work on every Alliance cluster with the
`StdEnv/2023` modules. The project mostly ran on **Fir** and **Trillium**; Rorqual and
Nibi behave like Fir. The differences are few, and each is shown below side by side, with
the lines that differ highlighted and numbered. The cluster pages are the reference for
anything not covered here: [Fir](https://docs.alliancecan.ca/wiki/Fir),
[Trillium](https://docs.alliancecan.ca/wiki/Trillium).

## At a glance

| | Fir (also Rorqual, Nibi) | Trillium |
|---|---|---|
| Filesystem | {term}`Lustre` | {term}`VAST` (NFS) |
| Job script | `shell/run_sim.sh` | `shell/run_sim_trillium.sh` |
| Launcher | `srun` | `mpirun`, ranks packed per {term}`L3 cache` |
| `sbatch` from a compute node | allowed | refused: always use `submit.sh` from a login node |
| `$PROJECT` | not defined (`/project/def-goluskin/<user>` used) | defined |
| VAST preload | absent, skipped | sourced by `setup.sh` |
| Cores per node | 192 | 192 |

## Job script header

::::{grid} 1 1 2 2
:gutter: 2

:::{grid-item-card} Fir: `run_sim.sh`
:class-card: cluster-card cluster-fir

```{code-block} bash
:emphasize-lines: 5,14

#SBATCH --output=job_sim_%j.out
#SBATCH --time=11:59:00
#SBATCH --nodes=27
#SBATCH --ntasks-per-node=192
#SBATCH --mem=0                       # (1)
#SBATCH --job-name=ra8pr-2.5
#SBATCH --account=rrg-goluskin

source params.conf
source "$SETUP_PATH"
source "$WORKFLOW_PATH"
load_modules
determine_ic
srun_sim                              # (2)
```
:::

:::{grid-item-card} Trillium: `run_sim_trillium.sh`
:class-card: cluster-card cluster-trillium

```{code-block} bash
:emphasize-lines: 5,10,11,12,15

#SBATCH --output=job_sim_%j.out
#SBATCH --time=11:59:00
#SBATCH --nodes=27
#SBATCH --ntasks-per-node=192
#SBATCH --cpus-per-task=1             # (1)
#SBATCH --job-name=ra8pr-2.5
#SBATCH --account=rrg-goluskin

source params.conf ...                # as on Fir
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores OMP_PROC_BIND=true  # (3)
export RANKS_PER_L3CACHE=$((8 / OMP_NUM_THREADS))
load_modules
determine_ic
trill_run_sim                         # (2)
```
:::
::::

{bdg-primary}`1` **Memory and cores.** On Fir the job asks for all of each node's memory
(`--mem=0`) rather than the default amount per core. Trillium schedules whole nodes (see
its page), and the cores per task are stated so that OpenMP and the binding agree.

{bdg-danger}`2` **Launcher.** Fir uses `srun`, which Slurm's own binding handles.
Trillium's MPI is launched with `mpirun`, which needs the placement spelled out (next
section).

{bdg-success}`3` **Thread placement (Trillium).** Each rank's threads stay on its own
cores. With one thread per rank (`C=1`), 8 ranks share each L3 cache.

## Launch line

::::{grid} 1 1 2 2
:gutter: 2

:::{grid-item-card} Fir: `srun_sim`
:class-card: cluster-card cluster-fir

```{code-block} bash
:emphasize-lines: 3,4

rb_args
set_deadline
srun -n "$N" -c "$C" \
    --cpu-bind=cores \
    python3 "$SIM_SCRIPT" \
    "${RB_ARGS[@]}"
```
:::

:::{grid-item-card} Trillium: `trill_run_sim`
:class-card: cluster-card cluster-trillium

```{code-block} bash
:emphasize-lines: 3,4,5,6,7,8

rb_args
set_deadline
mpirun -n "$N" --bind-to core \
    --map-by "ppr:8:l3cache:pe=1" \
    -x RBC_DEADLINE -x RBC_GRACEFUL \
    -x RBC_WATCHDOG -x RBC_WATCHDOG_ALL \
    -x HDF5_USE_FILE_LOCKING \
    -x OMP_NUM_THREADS \
    python3 "$SIM_SCRIPT" \
    "${RB_ARGS[@]}"
```
:::
::::

{bdg-danger}`2` On Trillium (shown with `C=1`; the script uses `$RANKS_PER_L3CACHE` and `$OMP_NUM_THREADS`), `--map-by ppr:8:l3cache:pe=1` places 8 ranks per L3
cache, each bound to one core; FFT-heavy spectral codes are limited by memory bandwidth,
and spreading ranks evenly over the caches uses all of it. The `-x` options forward the
variables that the simulation reads to ranks on other nodes, which `mpirun` does not do
by default. `workflow.sh` also sets `SLURM_EXPORT_ENV=ALL`, because Trillium's `sbatch`
implies `--export=NONE`, which would otherwise stop the module environment reaching
the ranks on other nodes ({ref}`export-env`).

## Submitting

`shell/submit.sh` from a login node works on both and is the only way on Trillium, whose
compute nodes cannot submit jobs. On Fir, a job script submitted with plain `sbatch` can
queue its own post-processing; `submit.sh` leaves a marker so it is never queued twice.

## Filesystem

| | Lustre (Fir) | VAST (Trillium) |
|---|---|---|
| What goes wrong | Metadata server overload when thousands of ranks create files at once | A node can cache "file does not exist" for a file another node just created |
| Symptom | slow or stuck first write with `PARA=virtual` at large rank counts | hang at the first write or at a new set ({ref}`vast-hang`) |
| Remedy | `PARA=mpio` (one shared file) if it occurs | {mod}`rbc.vast` (on by default, harmless on Lustre) |
| Restart loading | either loader | the fast loader matters most: see {ref}`n2-opens` |

## Paths

`BACKUP_DIRS` uses `$PROJECT` where it exists and `/project/def-goluskin/$USER` where it
does not (Fir). Everything else is relative to `$SCRATCH` and `$HOME`, which every
cluster defines.
