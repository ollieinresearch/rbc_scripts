# Glossary

Terms used throughout these pages. Anywhere a term is underlined, hovering over it shows
its entry here.

```{glossary}
:sorted:

Rayleigh–Bénard convection
  Flow of a fluid layer between a heated bottom plate and a cooled top plate, driven by
  buoyancy. Here: three-dimensional, horizontally periodic, with no-slip plates held at
  fixed temperatures. See {doc}`theory/equations`.

Rayleigh number
Ra
  Strength of the thermal driving relative to diffusion,
  $Ra = g\alpha_T\Delta T H^3/(\nu\kappa)$. Runs are labelled by $\log_{10} Ra$
  (`RA=8` in `params.conf` means $Ra = 10^8$).

Prandtl number
Pr
  Ratio of viscosity to thermal diffusivity, $Pr = \nu/\kappa$. Liquid metals have
  $Pr \sim 10^{-2}$, air about 0.7, water about 7. Labelled by $\log_{10} Pr$ on
  half-decades (`PR=-2.5`).

Nusselt number
Nu
  Total vertical heat flux divided by the flux that conduction alone would carry; 1 for a
  fluid at rest. Three exact identities give it from volume averages; see
  {doc}`theory/nusselt`. The scaling analysis uses $Nu - 1$, the convective part.

Reynolds number
Re
  $Re = U_\text{rms} H/\nu$ with the rms velocity over the volume; in free-fall units
  $Re = \sqrt{Ra/Pr}\,\langle|\boldsymbol u|^2\rangle^{1/2}$.

free-fall units
  Units with length $H$, temperature $\Delta T$ and velocity
  $U_f = \sqrt{g\alpha_T\Delta T H}$, so that time is in free-fall times $H/U_f$. Then
  $\nu = \sqrt{Pr/Ra}$ and $\kappa = 1/\sqrt{Ra\,Pr}$. All simulation times in this
  project are in these units. See {ref}`free-fall-units`.

aspect ratio
Γ
  Width over height of the box, $\Gamma = L_x/H = 2$ in every run here.

slice
  A line $Pr = Pr_0\,(Ra/10^5)^{a}$ through parameter space, labelled by its Prandtl number
  $Pr_0$ at $Ra = 10^5$. The main slices have $a = -1/2$; slices with $a = -1$ are drawn
  with `--ra-inv`. See {func}`rbc.cases.slice_c`.

Pr_0
c
  Label of a slice: its Prandtl number at $Ra = 10^5$, stored as $c = \log_{10} Pr_0$.

family
  All slices of one slope $a$. Family $-1/2$ is drawn solid with filled markers, family
  $-1$ dashed with open markers.

cases.csv
  The hand-edited registry of parameter points and their status (finished, running,
  planned). See {doc}`user-guide/cases-and-results`.

results.csv
  The table written by {mod}`rbc.collect`, one row per run directory. Kept out of git
  until the paper is published. Columns: {ref}`results-columns`.

production run
  The run used for a parameter point's results when there are several resolutions: the
  one pinned by `nz` in `cases.csv`, else the finest converged one.

res_check
  Column of `cases.csv`: 1 where a resolution study was done (drawn as a star).

resolution study
  Several runs of the same (Ra, Pr) at increasing resolution, to measure how Nu converges.
  See {doc}`theory/resolution`.

run directory
  One simulation at one resolution, e.g. `$SCRATCH/ollie_rb_data/3d/ra1e8/pr1e0/108_Gam2_2`.
  Holds `params.conf`, the job script and every output. See
  {doc}`user-guide/directories`.

set
  One file (or folder of process files) written by a Dedalus file handler, named
  `<handler>_s<N>.h5`. Each job starts new sets.

checkpoint
state
  A `state` set holding the full solution, written every `STATE_FREQ` wall-clock seconds
  and at the end of every job. Restarts start from the newest readable one.

virtual file
  Dedalus' default parallel output: every rank writes its own process file, and a small
  joint file presents them as one dataset. Fast to write; reading it naively opens every
  process file. See {doc}`theory/parallel`.

process mesh
  The 2D layout `MESHX` × `MESHY` of MPI ranks over the horizontal directions. Their
  product is the number of ranks `N`.

rank
  One MPI process. A full node of Fir or Trillium runs 192 ranks.

Chebyshev polynomials
  The vertical basis. The grid is the Gauss–Chebyshev points, clustered near the plates
  where the boundary layers are. See {ref}`chebyshev-nodes`.

dealiasing
  Evaluating products on a grid 3/2 times finer than the number of modes, so that the
  quadratic nonlinearity does not alias onto resolved modes.

tau method
  How Dedalus imposes boundary conditions on a spectral basis: extra unknowns (tau
  variables) multiply polynomials added to the equations. See {doc}`theory/equations`.

CFL condition
  The adaptive timestep limit $\Delta t \le C\,\Delta x/|u|$; `CFL_SAFETY` is $C$.

RK222
  The two-stage, second-order implicit–explicit Runge–Kutta timestepper used for all runs.

boundary layer
  Thin regions next to the plates where temperature (thermal boundary layer, thickness
  $\delta_T$) or horizontal velocity (viscous boundary layer, $\delta_u$) change rapidly.
  See {doc}`theory/boundary-layers`.

δ_T
  Thermal boundary-layer thickness, $1/(2Nu)$: where the wall tangent of the mean
  temperature reaches the bulk value 1/2.

δ_u
  Viscous boundary-layer thickness, from the peak of the rms horizontal velocity.

Prandtl–Blasius
  Laminar boundary-layer theory used for the reference curves of $\delta_u/\delta_T$.

statistical convergence
  Whether a time average has been taken over long enough, with the initial transient
  removed, for its sampling error to be small. See {doc}`theory/convergence`.

standard error
SE
  Standard deviation of a time average over repeated, independent records of the same
  length: the sampling error of a mean.

integrated autocorrelation time
τ_int
  $\int_0^\infty \rho(s)\,ds$ for the autocorrelation function $\rho$ of a time series;
  samples further apart than a few $\tau_\text{int}$ are effectively independent.

effective samples
n_eff
  $T/(2\tau_\text{int})$: the number of independent samples a record of length $T$ is
  worth.

MSER-5
  Rule for choosing where the initial transient ends: the start time that minimizes the
  marginal standard error of the remaining batch means. See {ref}`transient`.

section test
  The group's original convergence check: the averaging window is split into thirds, and
  the Nu values of the thirds must agree.

spectral decades
D
  How many decades the energy spectrum drops from its peak to its last modes. With $D$
  decades, Nu is expected to have about $D/2$ significant figures. See
  {doc}`theory/resolution`.

marginal spectrum
  Energy summed over all modes with a given index in one direction (Chebyshev degree, or
  horizontal wavenumber). The resolution rule uses the marginals.

shell spectrum
  Energy summed over shells of constant index radius; the 3D shells mix Fourier
  wavenumbers with Chebyshev degrees.

pile-up
  Energy accumulating in the last resolved modes; a sign of under-resolution.

local exponent
α
β
  Local slopes of $\log Q$ against $\log Pr$ (α) and $\log Ra$ (β) for $Q = Nu - 1$,
  $Nu$ or $Re$, written $Q \sim Pr^{\alpha} Ra^{\beta}$. See {doc}`theory/exponents`.

compensated plot
  $Nu$ divided by the power of $Ra$ that a theory predicts along a slice, so that the
  prediction is a horizontal line. See {doc}`theory/normalization`.

AICc
  Akaike's information criterion with the small-sample correction, used to choose the
  degree of a fitted surface.

Shishkina–Lohse
  The theory of the ultimate regime at small Prandtl number tested here, predicting
  $Nu \sim (Ra\,Pr)^{1/2}$ up to logarithms, which is $Nu \sim Ra^{1/4}$ along
  $Pr \propto Ra^{-1/2}$.

graceful stop
  Ending a run early in a way every rank agrees on, so a final checkpoint can be written.
  Triggered by the job's deadline or by `SIGUSR1`. See {ref}`graceful-stop`.

watchdog
  A timer that dumps the Python stack of every node's first rank if the run makes no
  progress for `WATCHDOG` seconds, then ends the job. See {ref}`watchdog`.

regrid
  Changing the resolution of a coefficient-space checkpoint by zero-padding or
  truncation, which is exact. See {doc}`theory/regridding`.

scratch
  Large, fast, unbacked-up storage (`$SCRATCH`) where runs live. Old files may be purged.

project space
  Backed-up group storage (`$PROJECT`, or `/project/def-goluskin/<user>` on Fir).

home
  Small, backed-up personal storage (`$HOME`).

Lustre
  The parallel filesystem of Fir, Rorqual and Nibi.

VAST
  The NFS-based filesystem of Trillium. Its client-side caching is the reason for
  {mod}`rbc.vast`.

L3 cache
  Shared cache of a group of 8 cores on Trillium's CPUs; ranks are placed per L3 cache.

login node
  Where you land with `ssh`; for editing and submitting jobs, not computing. Trillium
  only accepts `sbatch` here.

compute node
  Where jobs run.

Slurm
  The job scheduler on every Alliance cluster (`sbatch`, `squeue`, `scancel`, `sacct`).

module
  The Alliance software environment system (`module load python/3.11.5 ...`).

wheelhouse
  The Alliance's own collection of Python wheels, used by `pip install --no-index`.
```
