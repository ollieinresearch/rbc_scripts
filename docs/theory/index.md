# Theory and methods

The physics and numerics behind the scripts, with the derivations. Each page states its
results in the text; the working is in {fas}`chevron-down` **Derivation** boxes.

| Page | What it covers |
|---|---|
| {doc}`equations` | Boussinesq equations, free-fall units, boundary conditions, the Dedalus formulation and timestepping |
| {doc}`nusselt` | The three identities for Nu, the plate fluxes, Re and the rms vorticity |
| {doc}`boundary-layers` | Boundary-layer thicknesses, Prandtl–Blasius ratios, nodes inside the layers |
| {doc}`convergence` | When a time average is good enough: standard errors, the transient, the convergence rule |
| {doc}`resolution` | Spectral conventions, decades, what a resolution study measures |
| {doc}`exponents` | Local scaling exponents α and β from the grid of runs |
| {doc}`normalization` | How to place slices on a compensated plot |
| {doc}`regridding` | Why changing the resolution of a checkpoint is exact |
| {doc}`parallel` | Graceful stopping, the watchdog, and the filesystem and MPI problems solved along the way |

```{toctree}
:hidden:

equations
nusselt
boundary-layers
convergence
resolution
exponents
normalization
regridding
parallel
```
