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
(vortex tubes, and whether they sit in rising or sinking fluid). `volume`: volume
renderings with fixed transfer functions. `--spectra` adds the spectra at the same instant
below, so a burst of small-scale energy can be matched to what the flow is doing.
`post.sh` makes this movie automatically when `SNAPSHOTS=1`.

Movies are assembled with ffmpeg (loaded with the modules). Frames render in parallel
(`--jobs`), each in its own process; `--stride 2` halves the horizontal resolution if
memory runs short.

Finished pictures and movies worth keeping go in the {doc}`../gallery/media`.
