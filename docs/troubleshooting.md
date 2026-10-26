# Troubleshooting

Problems met so far, by symptom. Click one to see the cause and the fix.

## The job hangs

:::{dropdown} The log stops after "Starting loop" or at "First step done", on Trillium
**Cause.** Most likely a rank waiting in a collective that the others never reach. Two
known cases: the VAST file-creation race ({ref}`vast-hang`) and, in an earlier version of
the restart loader, mismatched field layouts ({ref}`fftw-deadlock`).

**What to do.** With the {ref}`watchdog <watchdog>` on (default), the job ends by itself after
`WATCHDOG` seconds and leaves stacks in `watchdog/`. The last lines of each stack show
where it was stuck: in `get_file` or `create_current_file` means file creation (check that
`VAST_PATCH` is not 0); in `distributor.py` `_plan` means the layouts disagree (run
`verify_load` at small scale). Send the stacks with any report.
:::

:::{dropdown} `sbatch` in a Trillium job does nothing, or the post-processing never starts
**Cause.** Trillium does not allow submitting jobs from compute nodes.

**Fix.** Submit with `bash $SCRIPTS_PATH/shell/submit.sh` from a login node, which queues the
post-processing as a dependent job. For a run already finished:
`cd RUN && sbatch $SCRIPTS_PATH/shell/post.sh`.
:::

:::{dropdown} A multi-node job starts, but ranks on other nodes fail to import modules
**Cause.** `--export=NONE` keeps the job's environment from the other nodes
({ref}`export-env`).

**Fix.** Already handled by `workflow.sh` (`SLURM_EXPORT_ENV=ALL`) and the `-x` options of
`trill_run_sim`. If a new environment variable must reach the simulation on Trillium, add
`-x NAME` there.
:::

:::{dropdown} The restart takes very long before the first step
**Cause.** Loading a `virtual` checkpoint through the joint file opens every process file
on every rank ({ref}`n2-opens`).

**Fix.** The fast loader avoids it when the process mesh is unchanged. If the log says
"decomposition changed; using Dedalus load_state", restore the previous `MESHX`/`MESHY`, or
accept one slow restart: the next checkpoint matches the new mesh.
:::

## Errors at start-up

:::{dropdown} `TypeError: 'float' object is not callable` in `add_equation`
**Cause.** A variable named `dt` in the simulation script. Dedalus parses the equations in
the script's local namespace, where `dt(...)` must be the time-derivative operator; a
number called `dt` replaces it.

**Fix.** The timestep variable is `step_dt`. Never name a variable after an operator used
in the equations (`dt`, `lap`, `grad`, `div`, `lift`, ...).
:::

:::{dropdown} `Cannot change global shape when loading coeff data`
**Cause.** Restarting from a coefficient-space checkpoint at a different resolution with
Dedalus' own loader.

**Fix.** Set the new `RES_Z` in `params.conf` and submit as usual: `determine_ic` regrids the
checkpoint first ({doc}`user-guide/restarts`).
:::

:::{dropdown} "no usable restart in state/" with `IC=1`
**Cause.** No readable checkpoint: none written yet, or the only ones were cut off
mid-write.

**Fix.** `python -m rbc.restart RUN` reports which sets it skipped. Start fresh (`IC=0`) or
restart from another run (`RESTART_FROM`).
:::

:::{dropdown} "WARNING: rbc package not installed"
**Fix.** `source $ENV_PATH/bin/activate && pip install --no-index -e $SCRIPTS_PATH`, or rerun
`shell/install_envs.sh`.
:::

## Errors during a run

:::{dropdown} The job dies mid-write with an HDF5 error about sizes or byte counts
**Cause.** Usually a full quota: HDF5 reports the failed write, not the reason.

**Fix.** `diskusage_report`. Delete snapshots or old checkpoints
(`python -m rbc.restart RUN --keep 4 --prune`), or move data to project space. The newest
checkpoint may be the corrupt one; the restart skips it automatically.
:::

:::{dropdown} "Re is not finite; stopping."
**Cause.** The solution blew up: the timestep was too large (CFL safety too high, or a
fixed `DT`), or the resolution far too low.

**Fix.** No final checkpoint is written, so the run can restart from the last good one.
Lower `CFL_SAFETY`, or increase the resolution.
:::

:::{dropdown} The timestep keeps shrinking
**Cause.** Velocities growing, normally during the initial transient; persistently, a sign of
under-resolution (grid-scale noise raises the velocity maximum).

**Fix.** Check the spectra for a pile-up ({doc}`theory/resolution`).
:::

:::{dropdown} "Graceful stop planned at ..." appears at the start of every job
That is the announcement of when the run will stop, not a stop. If the job is short, the
margin is capped at a tenth of it. Set `GRACEFUL_STOP=0` to turn it off for testing.
:::

## Results look wrong

:::{dropdown} The three Nu values disagree
Over a short window they differ by storage terms that shrink like $1/T$
({doc}`theory/nusselt`); if the window is long and they still differ, the vorticity form is
most sensitive to under-resolution. Check `converged_stat` and the spectra.
:::

:::{dropdown} `converged` is False although Nu hardly changes when the run continues
The thirds rule is noisy ({doc}`theory/convergence`): the spread of the thirds is about
three standard errors and varies a lot from run to run. Check the standard error and
`converged_stat` instead.
:::

:::{dropdown} Spectra of different resolutions disagree at low mode numbers
Those modes are under-resolved in the coarser run ({doc}`theory/resolution`). Spectra
written by the old FFT-based script also differed by a constant factor because they were
not divided by the number of grid points; recompute them with
{mod}`rbc.spectra.snapshots`.
:::

:::{dropdown} Negative values of $\lvert\omega\rvert^2$ in snapshots
It is written as a product evaluated on the output grid, which can undershoot slightly
where it is small; widespread negative values indicate aliasing from under-resolution. The
visualizations clip them unless `--keep-negative`.
:::

## Visualization

:::{dropdown} `XIO: fatal IO error 22 (Invalid argument) on X server ":99"` from `plot_snapshots`
A render process lost its connection to the virtual X server and was killed by Xlib; in a
process pool that stops every frame still pending. Display `:99` is the one
`pyvista.start_xvfb` always uses: older versions of {mod}`rbc.vis.render3d` called it from
every worker, so each of the `VIS_JOBS` workers raced to start `Xvfb :99`, and all of them
then depended on one server that no process watched.

The current version starts a single Xvfb in the parent on a free display (its messages go
to `visualization/xvfb.log`), checks that rendering works before starting workers, and
replaces each worker every few frames. If a worker still dies, the pool is rebuilt, Xvfb
restarted if it died, and the missing frames tried again (`--retries`). Frames already
written are kept, so running the same command again, or the next post-processing job,
renders only what is missing. If the server keeps dying, look in `xvfb.log`, and try
fewer `--jobs` or `--stride 2` (memory).
:::

:::{dropdown} `vtkOpenGLVolume...Table: ... does not support the required texture size of 32768`
From `--mode volume`. VTK sizes the colour and opacity lookup textures from the range
the transfer functions and the data cover, relative to the spacing of their points; when
they cover very different ranges it asks for more than OpenGL allows and falls back to a
coarser table. One version of {mod}`rbc.vis.render3d` did exactly that: pyvista 0.44
rescales volume scalars to 0–255 before VTK sees them, and that version then set an
opacity function in data units (−3 to 3) on VTK directly, so the two did not line up,
nearly every voxel got the largest opacity, and the boxes rendered solid. The volume
mode now passes its transfer functions through pyvista, as the original
`plotting_v3.py` did, and neither problem occurs.
:::

## MPI binding

:::{dropdown} Runs much slower than expected with `mpirun`
Check the binding with `--report-bindings`. With Open MPI, `--map-by :PE=$C` silently drops
the PE modifier; write the mapping in full (`ppr:8:l3cache:pe=1` on Trillium,
`slot:PE=$C --bind-to core` elsewhere), and make the `#SBATCH` lines agree:
`--ntasks=N --cpus-per-task=C`, not `--ntasks=N*C`.
:::
