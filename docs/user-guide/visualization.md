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
python -m rbc.vis.slices RUN_PR1 RUN_PR001 --fields T w --heights mid bl
python -m rbc.vis.slices RUN_PR1 RUN_PR001 --movie --jobs 16
```

One row per run, one column per (field, height). Heights: `mid` ($z = 1/2$), `bl`
($z = \delta_T/2 = 1/(4Nu)$, inside the thermal boundary layer), `bl-edge` ($\delta_T$),
or a number. Side by side, $Pr = 1$ and $Pr = 0.01$ show how much finer the velocity is
than the temperature at small Pr.

## Joint PDF of w and θ

```bash
python -m rbc.vis.jointpdf RUN_PR1 RUN_PR001 --height mid --last 20 --jobs 8
```

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
python -m rbc.vis.render3d RUN --mode iso --spectra --movie --jobs 16
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
(`--jobs`), in worker processes that are replaced every `--frames-per-process` frames;
`--stride 2` halves the horizontal resolution if memory runs short. On a compute node
without a display, one virtual X server (Xvfb) is started for the whole command and
stopped at the end.

Frames already in `visualization/frames_<mode>/` are kept, so an interrupted movie, or the
next post-processing job after more snapshots, renders only the new ones; a frame that
fails is retried (`--retries`), and the MP4 is made only once every frame exists. The
spectra panel's time mean is the one at the moment each frame was rendered, so for a final
movie re-render everything with `--redo`.

Finished pictures and movies worth keeping go in the {doc}`../gallery/media`.
