# Your first run

A complete small simulation: $Ra = 10^5$, $Pr = 1$, 48 Chebyshev modes, one node for an
hour. It is cheap, and it exercises every part of the workflow: submission, a fresh start,
the final checkpoint, post-processing, backups, a restart and the figures.

## 1. Make the run directory

Runs live under `$SCRATCH/ollie_rb_data/3d/ra<Ra>/pr<Pr>/<label>` (naming rules in
{doc}`../user-guide/directories`):

```bash
RUN=$SCRATCH/ollie_rb_data/3d/ra1e5/pr1e0/48_Gam2_2
mkdir -p $RUN && cd $RUN
cp $SCRATCH/rbc_scripts/shell/params.conf .
cp $SCRATCH/rbc_scripts/shell/run_sim.sh .            # on Trillium: run_sim_trillium.sh
```

## 2. Set the parameters

Edit `params.conf`; the lines to change are highlighted (every variable is described in
{doc}`../user-guide/parameters`):

```{code-block} bash
:emphasize-lines: 2,4,5,7,8,9,11
# Simulation
export RES_Z=48                  # Chebyshev modes in z
export RES_HOR=$((RES_Z * 2))    # Fourier modes in x and y: 96 (aspect ratio 2)
export MESHX=16                  # 16 x 12 = 192 ranks = one full node
export MESHY=12
export N=$((MESHX * MESHY))
export RA=5                      # log10 Ra
export PR=0                      # log10 Pr
export SIM_TIME=50               # free-fall times to simulate
...
export IC=0                      # start from rest (small random temperature noise)
```

and the `#SBATCH` lines at the top of the job script:

```{code-block} bash
:emphasize-lines: 2,3,5
#SBATCH --output=job_sim_%j.out
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --job-name=ra5pr0
#SBATCH --account=rrg-goluskin
```

```{dropdown} Why these numbers
:class-container: derivation
- `RES_Z=48` is ample at $Ra = 10^5$: the thermal boundary layer is about $1/(2Nu)$ thick,
  roughly 0.12, and 48 Gauss–Chebyshev nodes put about ten of them inside it
  ({ref}`chebyshev-nodes`).
- The process mesh must multiply to the number of ranks, and each direction of the mesh
  should divide the number of Fourier modes (96 / 16 = 6, 96 / 12 = 8) so every rank holds
  the same amount of data.
- One hour is far more than needed; with Slurm, short jobs start sooner.
- `SIM_TIME=50` is too short for converged statistics (in practice nothing converges in
  less than about 150 free-fall times), which is why step 6 restarts the run.
```

## 3. Submit

From a {term}`login node`, inside the run directory:

```bash
bash $SCRATCH/rbc_scripts/shell/submit.sh
```

It picks the right job script for the cluster, submits it, and queues the
post-processing job to start when the simulation ends, whatever its exit state:

```text
simulation job 4123456 (run_sim.sh), post-processing job 4123457 (starts after 4123456)
```

## 4. Watch it

```bash
squeue -u $USER                  # PD = pending, R = running
tail -f job_sim_4123456.out      # the simulation's log
```

A healthy log shows, in order: the parameters, a disk-quota summary, `Starting
fresh`, the {term}`graceful stop` time, Dedalus building its matrices, `Starting loop`,
`First step done`, and then every 500 iterations a line such as

```text
Iteration=1001, t=4.12000, dt=4.1e-03, Re=58.31
```

$Re$ grows from zero as convection sets in, then fluctuates about a steady value. If the
log stops changing for a long time, see {doc}`../troubleshooting`.

## 5. Read the results

When the post-processing job has run, `outputs/` holds:

| File | What it is |
|---|---|
| `info.txt`, `info.pdf` | Nu from three identities, Re, and the {doc}`convergence diagnostics <../theory/convergence>` |
| `profiles.pdf`, `profiles.json` | mean temperature and velocity profiles, boundary-layer thicknesses |
| `cumulative_spectra_*.pdf`, `spec_range.txt` | time-averaged spectra and {term}`spectral decades` |
| `jobs.csv` | ranks, wall time and core-hours of each job |

The run has also been copied to `$HOME/ollie_rb_data/3d/ra1e5/pr1e0/48_Gam2_2` and to
project space ({doc}`../user-guide/storage`), and `data/results.csv` in the repository
has a row for it.

Open `info.pdf`. The top-left panel is the instantaneous Nu: it rises from 1 (pure
conduction) through an overshoot to a statistically steady state. With 50 time units the
three Nu values in `info.txt` will not yet agree to the precision needed, and
`converged_stat` will be `False`.

## 6. Restart for longer

Edit `params.conf`:

```bash
export IC=1             # restart from the newest checkpoint
export SIM_TIME=200     # 200 more time units
export AVG_TIME=auto    # let the analysis find where the transient ends
```

and submit again. The log now says `Restarting from write 0 at t = 50.0...; running to t =
250.0...`, and the run continues exactly where it stopped. See
{doc}`../user-guide/restarts` for what else a restart can do (change resolution, start from
another run).

## 7. Make figures

```bash
cd $SCRATCH/rbc_scripts
python -m rbc.plots.all --png
```

writes every paper figure to `figures/`. With one run most are empty; the
{doc}`../gallery/index` shows what they look like with a full set of runs.

## Where next

- {doc}`../user-guide/workflow`: what each script did in this run.
- {doc}`../theory/convergence`: how to decide a run is long enough.
- {doc}`../theory/resolution`: how to decide it is well enough resolved.
