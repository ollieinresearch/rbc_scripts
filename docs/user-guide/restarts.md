# Restarts and resolution changes

Runs take many jobs. Each job ends with a {term}`checkpoint`, and the next starts from it.

## Continuing a run

With `IC=1` (the template's setting) and a checkpoint to start from, `determine_ic`
calls {mod}`rbc.restart`, which:

1. lists `state/state_s*.h5` newest first;
2. skips any set it cannot read: a job killed while writing leaves a set that exists but
   is incomplete, which the old `ls | tail -1` would have chosen. For {term}`virtual
   files <virtual file>` it checks that every process file is present, because a missing
   one reads back as zeros rather than as an error;
3. takes the newest readable write (or, with `START_TIME`, the write nearest that time);
4. links it as `restart/restart.h5` and prints its index and time.

The run then goes on to `restart time + SIM_TIME`. Continuity is exact: the checkpoint
holds the full state in coefficient space, and the first timestep is the one in use when
it was written.

The {term}`process mesh` may change between jobs. The loader in {mod}`rbc.vast` reads the
checkpoint's map of which process file holds which block, and each rank opens only the
files that overlap its own part of the domain: one when the mesh is unchanged, a few when
it changed. The log says `Loaded state from process files (at most k per rank)`. It falls
back to Dedalus' own loader, saying why, only for a plain (non-virtual) file such as a
regridded restart, a change of resolution, or a gap in the map ({ref}`n2-opens`).

## Starting over

`IC=0` starts from rest with small random temperature noise
($10^{-5} z(1-z)$, low-pass filtered), deleting `restart/`. Existing output is not
deleted; new sets are numbered after it, so remove old output first if the run should
really start again.

## Starting from another run

`RESTART_FROM` seeds a new run directory with another run's checkpoint: usually a lower
resolution of the same (Ra, Pr), or a nearby (Ra, Pr) whose flow is a better initial
condition than rest:

```bash
export IC=1
export RESTART_FROM="$PWD/../64_Gam2_2"
```

The newest checkpoint of that run is used (or the one nearest `START_TIME`). If its grid
differs from this run's (`RES_HOR`, `RES_Z`), it is regridded into `restart/restart.h5`
first.

`RESTART_FROM` only matters while this run has no checkpoints of its own. From the second
job on, the run continues from its own `state/` and the setting is ignored (the log says
so), so it can stay in `params.conf`.

## Changing resolution

A resolution change always makes a new run directory, named after its new $N_z$, seeded
from the old run. To go from 64 to 128 modes:

```bash
mkdir ../128_Gam2_2 && cd ../128_Gam2_2
bash $SCRATCH/rbc_scripts/shell/sim_ready.sh
# in params.conf:
export RES_Z=128          # RES_HOR follows as 2 * RES_Z
export IC=1
export RESTART_FROM="$PWD/../64_Gam2_2"
```

`determine_ic` passes the new shape to {mod}`rbc.restart`, which sees that the
checkpoint is 128 × 128 × 64 and regrids it with {mod}`rbc.regrid` to 256 × 256 × 128.
Refining is exact (the same function, written with more coefficients), and coarsening
removes only the highest modes; the proof is in {doc}`../theory/regridding`. Without
this step Dedalus refuses to load a coefficient-space checkpoint of another size
("Cannot change global shape when loading coeff data").

Changing `RES_Z` inside an existing run is refused, because the run's statistics and
checkpoints would mix two resolutions:

```text
restart: checkpoint state_s12.h5 has shape (128, 128, 64) but RES_HOR/RES_Z give
(256, 256, 128). Change resolution in a new run folder with RESTART_FROM=<this run>
```

After a resolution change the statistics need a short new transient (the small scales
adjust within a few time units), so start the averaging window after it.

By hand:

```bash
python -m rbc.regrid state/state_s12.h5 restart/regridded.h5 --nx 256 --nz 128
```

## Keeping and pruning checkpoints

```bash
python -m rbc.restart RUN --keep 4            # list the 4 newest readable checkpoints
python -m rbc.restart RUN --keep 4 --prune    # and delete older ones
```

Pruning deletes nothing unless that many readable checkpoints exist, and never deletes
the one `restart/restart.h5` points to. `backup_and_prune` does this automatically after
a successful backup ({doc}`storage`).
