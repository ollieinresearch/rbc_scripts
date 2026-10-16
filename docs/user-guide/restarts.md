# Restarts and resolution changes

Runs take many jobs. Each job ends with a {term}`checkpoint`, and the next starts from it.

## Continuing a run

With `IC=1` (the template's setting), `determine_ic` calls {mod}`rbc.restart`, which:

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

The {term}`process mesh` may change between jobs. If it does, the fast loader in
{mod}`rbc.vast` finds that the process files do not match the new decomposition and
falls back to Dedalus' own loader, which redistributes the data; the log says so.

## Starting over

`IC=0` starts from rest with small random temperature noise
($10^{-5} z(1-z)$, low-pass filtered), deleting `restart/`. Existing output is not
deleted; new sets are numbered after it, so remove old output first if the run should
really start again.

## Starting from another run

Set `RESTART_FROM` to another run directory, usually a lower resolution of the same (Ra,
Pr), or a nearby (Ra, Pr) whose flow is a better initial condition than rest:

```bash
export IC=1
export RESTART_FROM="$PWD/../64_Gam2_2"
```

The newest checkpoint of that run is used. If its grid differs from this run's
(`RES_HOR`, `RES_Z`), it is regridded into `restart/restart.h5` first.

## Changing resolution

To raise a run from 64 to 128 modes in place:

```bash
export RES_Z=128          # RES_HOR follows as 2 * RES_Z
export IC=1
```

`determine_ic` passes the new shape to {mod}`rbc.restart`, which sees that the
checkpoint is 128 × 128 × 64 and regrids it with {mod}`rbc.regrid` to 256 × 256 × 128.
Refining is exact (the same function, written with more coefficients), and coarsening
removes only the highest modes; the proof is in {doc}`../theory/regridding`. Without
this step Dedalus refuses to load a coefficient-space checkpoint of another size
("Cannot change global shape when loading coeff data").

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
