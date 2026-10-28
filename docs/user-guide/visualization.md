# Visualizations

Pictures and movies of the fields themselves, from snapshots. They need
`SNAPSHOTS=1` in `params.conf` for at least part of a run (snapshots are large; a short
stretch is usually enough).

## Fixed colour scales

Every field is shown in units that mean the same thing at every Pr and in every frame
({mod}`rbc.vis.fields`):

| Field | Shown as | Why |
|---|---|---|
| temperature | $\theta/\theta_\text{rms}(z)$ with $\theta = T - \overline T(z)$ | removes the mean profile, so plumes stand out at every height |
| vertical velocity | $w/w_\text{rms}(z)$ | |
| vorticity | $\lvert\omega\rvert/\omega_\text{rms}$ on a log scale | $\omega_\text{rms} = \sqrt{(Nu-1)/Pr}$ exactly ({ref}`derivation <nu-vorticity>`), so no frame-by-frame scaling is needed |

The profiles come from `outputs/profiles.npz`, so run the post-processing first.
$\lvert\omega\rvert^2$ is written as a product on the output grid and can dip slightly
below zero where it is small; those points are clipped, or kept and marked with
`--keep-negative`.

## Horizontal slices

```bash
python -m rbc.vis.slices RUN_A RUN_B RUN_C --fields T w --heights mid bl
python -m rbc.vis.slices RUN_A RUN_B RUN_C RUN_D --movie
python -m rbc.vis.slices RUN_A RUN_B RUN_C --write -1 -1 40            # a chosen snapshot per run
python -m rbc.vis.slices RUN_A RUN_B RUN_C --labels "Pr = 1" "Pr = 10^{-2}" "Pr = 10^{-3.5}"
```

Any number of runs, one row each in the order given; one column per (field, height).
Heights: `mid` ($z = 1/2$), `bl` ($z = \delta_T/2 = 1/(4Nu)$, inside the thermal boundary
layer), `bl-edge` ($\delta_T$), or a number. Side by side, $Pr = 1$ and $Pr = 0.01$ show how
much finer the velocity is than the temperature at small Pr.

Each row is labelled with its run and the time of its snapshot. Rows at the same (Ra, Pr)
also get their $N_z$, so resolutions can be compared; `--labels` replaces the labels
(one per run, TeX allowed). `--write` takes one snapshot index for every run, or one per
run (negative counts from the end).

Runs seldom have the same number of snapshots, so a movie pairs them by position:
`--align end` (the default) uses the last $n$ snapshots of every run, $n$ being the
shortest run's count, and `--align start` the first $n$. The log prints how many each run
has. Frames already written are kept, so an interrupted movie resumes; `--redo` starts
over, and frames made with other runs or options are cleared automatically (the settings
are kept in `settings.json` in the frames folder).

Each plane is Fourier-interpolated to at least `--points` points per side (512 by
default; 0 turns it off). The fields are Fourier series in $x$ and $y$, so zero-padding
their coefficients gives the same function on a finer grid: nothing is invented, but a
coarse run is drawn smoothly instead of as blocks (a run with $N_z = 64$ has only 128
points across). Planes are drawn with antialiasing, so a fine run with more points than
the panel has pixels is averaged rather than subsampled, which would make it look grainy.
Movie frames are written at `--dpi` (300); at the default size a panel is then about 350
pixels wide. In `workflow.sh`, `plot_slices` makes this movie for one run and
`plot_compare_slices RUN_A RUN_B RUN_C ...` for any number.

## Joint PDF of w and θ

```bash
python -m rbc.vis.jointpdf RUN_A RUN_B RUN_C --height mid --last 20
```

Any number of runs, one panel each, three to a row (`--columns`), with the same labels as
the slices (`--labels`).

The joint probability density of $(w/w_\text{rms}, \theta/\theta_\text{rms})$ on a plane,
accumulated over snapshots, and the correlation coefficient
$C(z) = \langle w\theta\rangle_h/\sqrt{\langle w^2\rangle_h\langle\theta^2\rangle_h}$.
Mass along the diagonal (warm fluid rising, cool fluid sinking) carries the heat; a round
blob means $w$ and $\theta$ are unrelated. Since
$\langle w\theta\rangle_h = (Nu-1)/\sqrt{Ra\,Pr}$ away from the walls, a smaller $C$ at
small Pr with the same flux means $w$ has variance that $\theta$ does not share.

## 3D renders and movies

Uses the visualization environment (pyvista):

```bash
source $SCRATCH/ded3_vis/bin/activate
python -m rbc.vis.render3d RUN --mode iso --spectra --movie
```

`iso`: temperature isosurfaces at $\theta = \pm1.5\,\theta_\text{rms}$ (warm red, cool
blue: plumes) beside vorticity isosurfaces at $3\,\omega_\text{rms}$ coloured by $w$
(vortex tubes, and whether they sit in rising or sinking fluid).

`volume`: the volume renderings of the original `plotting_v3.py`, with its transfer
functions. Temperature $T$ on $[0, 1]$ is transparent within about $0.2$ of $1/2$ and
opaque beyond about $0.25$, so only the hot and cold extremes (plumes and boundary layers)
show; $|\omega|^2$, on an asinh scale over its range in each frame, is opaque only in the
top fifth of that range, the strongest filaments. Both use `jet` (`--cmap-T`,
`--cmap-omega`), and each colour bar has its opacity curve drawn across it. The Chebyshev
points are resampled to a uniform $z$ grid; `--stretch-z` places them evenly as the
original did, which makes the boundary layers look thicker. Because the vorticity range is
taken frame by frame, its colours are relative to that frame. `--spectra` adds the spectra at the same instant
below, so a burst of small-scale energy can be matched to what the flow is doing.
`post.sh` makes this movie automatically when `SNAPSHOTS=1`.

Movies are assembled with ffmpeg (loaded with the modules). Frames render in parallel
({ref}`vis-workers`), in worker processes that are replaced every `--frames-per-process`
frames; `--stride 2` quarters the memory of a frame if a large run needs it. On a compute node
without a display, one virtual X server (Xvfb) is started for the whole command and
stopped at the end.

Frames already in `visualization/frames_<mode>/` are kept, so an interrupted movie, or the
next post-processing job after more snapshots, renders only the new ones; a frame that
fails is retried (`--retries`), and the MP4 is made only once every frame exists. Frames
made with other options (mode, thresholds, stride, ...) are cleared automatically. The
spectra panel's time mean is the one at the moment each frame was rendered, so for a final
movie re-render everything with `--redo`.

(vis-workers)=
## How many frames at once

Every visualization renders frames (or, for the joint PDF, snapshots) in separate worker
processes. With `--jobs 0`, the default and `VIS_JOBS=0` in `params.conf`, the number of
workers is the smallest of three limits ({mod}`rbc.vis.parallel`), and the log says which
one applied:

| Limit | Value |
|---|---|
| frames | no more workers than frames to render |
| cores | the cores the job may use (`os.sched_getaffinity`, which follows the Slurm allocation: 192 in `post.sh`) |
| memory | 70% of the memory the job may still use (the cgroup limit, or what the node has free), divided by one frame's estimated peak |

The memory of one frame was measured on two grid sizes and rounded up:

| Command | One frame | Example: $432^2 \times 216$ grid |
|---|---|---|
| `rbc.vis.slices` | 0.4 GB + 64 B per horizontal point and per plane, for the read and the image | 0.7 GB for 3 runs × 4 panels |
| `rbc.vis.jointpdf` | 0.4 GB + 1 kB per horizontal point | 0.6 GB |
| `rbc.vis.render3d --mode volume` | 0.6 GB + 100 B per grid point (after `--stride`) | 4.4 GB |
| `rbc.vis.render3d --mode iso` | 0.6 GB + 150 B per grid point | 6.2 GB |

So slices and joint PDFs use every core of a node, while 3D renders of large runs are
limited by memory: a $1024^2 \times 512$ run needs about 76 GB per isosurface frame, so
about seven at once on a node with 768 GB, or about 28 with `--stride 2` (19 GB each). More workers
than that would only be killed for memory. The cores left over go to each render worker's
software OpenGL (`LP_NUM_THREADS`, the cores divided by the workers), which otherwise
starts one thread per core in every worker.

A positive `--jobs` (or `VIS_JOBS`) is used as given, up to the number of frames. Lower it
when the job shares a node, or when the run's files are on a filesystem that struggles
with many readers (each frame of a `virtual` snapshot set opens every process file that
holds a part of the planes it reads).

The parent process reads the list of snapshots, and for `render3d` the spectra, once, and
passes each worker only what its frame needs.

Finished pictures and movies worth keeping go in the {doc}`../gallery/media`.
