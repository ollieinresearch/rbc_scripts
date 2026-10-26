# History and design choices

## The project

The project began in 2020 with David Goluskin and an earlier student, with the aim of
building a dataset of two- and three-dimensional Rayleigh–Bénard simulations at small
Prandtl number, and of using it to test the Shishkina–Lohse theory of heat transport in
that regime.

```{admonition} To fill in
:class: todo-note

The name of the first student, and the main milestones since 2020.
```

The work now concentrates on three dimensions. Alongside the main question, whether Nu
follows the predicted scaling along slices $Pr \propto Ra^{-1/2}$, two subprojects grew out
of it:

- **how to judge resolution**: what spectral diagnostics say about the accuracy of Nu
  ({doc}`../theory/resolution`);
- **how to estimate and test exponents** from a grid of runs, and so what it takes to
  confirm or refute a predicted scaling ({doc}`../theory/exponents`,
  {doc}`../theory/normalization`).

In 2026 the scripts were reorganized into the `rbc` package documented here.

## Why things are the way they are

| Choice | Reason | More |
|---|---|---|
| Dedalus | spectral accuracy, implicit–explicit timestepping, and equations written as text | {doc}`../theory/equations` |
| Horizontally periodic box | no side-wall layers to resolve; efficient Fourier bases | {doc}`../theory/equations` |
| Γ = 2 | a common choice for periodic boxes (to fill in: the original reason) | |
| Free-fall units | velocities and timesteps of order one at every (Ra, Pr) | {ref}`free-fall-units` |
| Slices $Pr \propto Ra^{-1/2}$ | the Shishkina–Lohse prediction is a single power of Ra along them | {doc}`../theory/exponents` |
| RK222 | SBDF2 gave oscillatory spectra and needed smaller timesteps | {ref}`timestepping` |
| Checkpoints in coefficient space | exact restarts, (3/2)³ times smaller, resolution changes by zero-padding | {doc}`../theory/regridding` |
| Spectra computed during the run | far cheaper than saving snapshots to compute them later | {doc}`../theory/resolution` |
| One checkpoint per set | a job killed mid-write damages one set only, which the restart skips | {doc}`../user-guide/restarts` |
| `Nu − 1` for exponents | the convective part, which the theories describe | {doc}`../theory/nusselt` |
| Logs of Ra and Pr on half-decades as keys | exact comparisons; no `3.16e-3` versus `3e-3` mismatches | {mod}`rbc.cases` |

(old-to-new)=
## From the old scripts

Run directories made before 2026 used separate scripts. Where their code went:

| Old | Now |
|---|---|
| `rayleigh_benard_script.py`, `fixed_get_file.py` | `sim/rayleigh_benard.py`, {mod}`rbc.vast` |
| `analysis_v3.py` | {mod}`rbc.analysis.nusselt` |
| `h_analysis.py` | {mod}`rbc.analysis.profiles` |
| `coeff_spectra.py`, `spectra.py`, `spectra_cumulative.py`, `spectra_compare.py` | {mod}`rbc.spectra.core`, {mod}`~rbc.spectra.snapshots`, {mod}`~rbc.spectra.cumulative`, {mod}`~rbc.spectra.compare` |
| `initial_condition.py`, `create_restart` | {mod}`rbc.restart` |
| `overall_analysis.py`, `redo_analysis.py` | {mod}`rbc.collect` |
| `param_grid.py`, `rbc_slice.py`, `contour.py`, `res_scaling.py` | {mod}`rbc.plots.param_grid`, {mod}`~rbc.plots.scaling`, {mod}`~rbc.plots.contour`, {mod}`~rbc.plots.resolution` |
| `plotting_v3.py` | {mod}`rbc.vis.render3d` |
| `power*.py` | retired |
| post-processing inside `run_sim.sh`; `post_process.sh`, `plotting.sh` | `shell/post.sh`, a dependent one-node job |
