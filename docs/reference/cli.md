# Command-line reference

Every command's options, generated from the scripts themselves, so they are always current.
Each command also prints this with `--help`. Shell functions that call them are in
{doc}`../user-guide/workflow`.

## Simulation

(cli-rayleigh-benard)=
### `sim/rayleigh_benard.py`

The Dedalus simulation; launched by the job scripts, which build these options from `params.conf` (`rb_args` in `shell/workflow.sh`).

```{argparse}
:module: rayleigh_benard
:func: get_parser
:prog: python sim/rayleigh_benard.py
```

## Run management

(cli-rbc-restart)=
### `rbc.restart`

Choose the restart checkpoint and link it; keep or prune old checkpoints. See {doc}`../user-guide/restarts`.

```{argparse}
:module: rbc.restart
:func: get_parser
:prog: python -m rbc.restart
```

(cli-rbc-regrid)=
### `rbc.regrid`

Change the resolution of a coefficient-space checkpoint. See {doc}`../theory/regridding`.

```{argparse}
:module: rbc.regrid
:func: get_parser
:prog: python -m rbc.regrid
```

(cli-rbc-resources)=
### `rbc.resources`

Job history and core-hours of runs. See {ref}`post-processing <post-processing-cost>`.

```{argparse}
:module: rbc.resources
:func: get_parser
:prog: python -m rbc.resources
```

## Analysis of one run

(cli-rbc-analysis-nusselt)=
### `rbc.analysis.nusselt`

Nu, Re and convergence diagnostics. See {doc}`../theory/nusselt` and {doc}`../theory/convergence`.

```{argparse}
:module: rbc.analysis.nusselt
:func: get_parser
:prog: python -m rbc.analysis.nusselt
```

(cli-rbc-analysis-profiles)=
### `rbc.analysis.profiles`

Time-averaged profiles and boundary layers. See {doc}`../theory/boundary-layers`.

```{argparse}
:module: rbc.analysis.profiles
:func: get_parser
:prog: python -m rbc.analysis.profiles
```

(cli-rbc-spectra-cumulative)=
### `rbc.spectra.cumulative`

Time-averaged spectra and decades. See {doc}`../theory/resolution`.

```{argparse}
:module: rbc.spectra.cumulative
:func: get_parser
:prog: python -m rbc.spectra.cumulative
```

(cli-rbc-spectra-compare)=
### `rbc.spectra.compare`

Spectra of every resolution of one (Ra, Pr) overlaid.

```{argparse}
:module: rbc.spectra.compare
:func: get_parser
:prog: python -m rbc.spectra.compare
```

(cli-rbc-spectra-snapshots)=
### `rbc.spectra.snapshots`

Spectra from snapshots (older runs) and per-write spectra plots.

```{argparse}
:module: rbc.spectra.snapshots
:func: get_parser
:prog: python -m rbc.spectra.snapshots
```

## All runs

(cli-rbc-collect)=
### `rbc.collect`

Collect every run into `results.csv`. See {doc}`../user-guide/cases-and-results`.

```{argparse}
:module: rbc.collect
:func: get_parser
:prog: python -m rbc.collect
```

## Paper figures

(cli-rbc-plots-all)=
### `rbc.plots.all`

Every figure below.

```{argparse}
:module: rbc.plots.all
:func: get_parser
:prog: python -m rbc.plots.all
```

(cli-rbc-plots-param-grid)=
### `rbc.plots.param_grid`

Map of the parameter points and slices.

```{argparse}
:module: rbc.plots.param_grid
:func: get_parser
:prog: python -m rbc.plots.param_grid
```

(cli-rbc-plots-scaling)=
### `rbc.plots.scaling`

Compensated Nu, Nu − 1 or Re along slices. See {doc}`../theory/normalization`.

```{argparse}
:module: rbc.plots.scaling
:func: get_parser
:prog: python -m rbc.plots.scaling
```

(cli-rbc-plots-fixed)=
### `rbc.plots.fixed`

Nu − 1 against Pr at fixed Ra and against Ra at fixed Pr.

```{argparse}
:module: rbc.plots.fixed
:func: get_parser
:prog: python -m rbc.plots.fixed
```

(cli-rbc-plots-exponents)=
### `rbc.plots.exponents`

Local exponents α and β. See {doc}`../theory/exponents`.

```{argparse}
:module: rbc.plots.exponents
:func: get_parser
:prog: python -m rbc.plots.exponents
```

(cli-rbc-plots-boundary-layers)=
### `rbc.plots.boundary_layers`

Boundary-layer ratio and nodes per layer.

```{argparse}
:module: rbc.plots.boundary_layers
:func: get_parser
:prog: python -m rbc.plots.boundary_layers
```

(cli-rbc-plots-resolution)=
### `rbc.plots.resolution`

Resolution studies and the decades framework.

```{argparse}
:module: rbc.plots.resolution
:func: get_parser
:prog: python -m rbc.plots.resolution
```

(cli-rbc-plots-contour)=
### `rbc.plots.contour`

Contours of Nu or Re over the (Ra, Pr) plane.

```{argparse}
:module: rbc.plots.contour
:func: get_parser
:prog: python -m rbc.plots.contour
```

## Visualization

(cli-rbc-vis-slices)=
### `rbc.vis.slices`

Horizontal slices of several runs, and movies.

```{argparse}
:module: rbc.vis.slices
:func: get_parser
:prog: python -m rbc.vis.slices
```

(cli-rbc-vis-jointpdf)=
### `rbc.vis.jointpdf`

Joint PDF of w and θ, and their correlation profile.

```{argparse}
:module: rbc.vis.jointpdf
:func: get_parser
:prog: python -m rbc.vis.jointpdf
```

(cli-rbc-vis-render3d)=
### `rbc.vis.render3d`

3D isosurface or volume renders, and movies (visualization environment).

```{argparse}
:module: rbc.vis.render3d
:func: get_parser
:prog: python -m rbc.vis.render3d
```
