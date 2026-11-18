# Spectral resolution

A spectral method represents each field by a finite set of coefficients. When the flow is
resolved, the coefficients fall off rapidly with mode number, and the error of anything
computed from the field is set by the size of the coefficients that were left out. The
diagnostics here measure that fall-off, and resolution studies measure its effect on Nu
directly.

## What the spectra are

{mod}`rbc.spectra.core` computes spectra from Dedalus' coefficients without any
interpolation or windowing; its docstring gives the conventions in full. In short:

- each coefficient carries the mean square it contributes to the field ("energy per
  coefficient"), so that summing all of them gives the mean of $f^2$ over the grid
  exactly (a Parseval identity);
- Dedalus' coefficients do not depend on the grid size, so spectra of runs at different
  resolutions coincide at low mode numbers. (An earlier FFT-based script did not divide by
  the number of grid points, which shifted spectra of different resolutions by
  $\log_{10}(N_\text{hi}/N_\text{lo})$ decades.)
- the spectra are sums over bins, never averages, so slopes have their usual meaning.

The spectra saved every `SPEC_FREQ`, for velocity and for the temperature perturbation $b$:

| Name | Binned by | Use |
|---|---|---|
| vertical marginal `nz` | Chebyshev degree $n$, summed over $k_x, k_y$ | truncation in $z$ |
| horizontal marginal `h` | $k_x$ (and $k_y$), summed over the rest; mean of x and y | truncation in $x, y$ |
| horizontal shells `perp` | $\mathrm{rint}\sqrt{k_x^2 + k_y^2}$ | the horizontal cascade |
| 3D index shells `3d` | $\mathrm{rint}\sqrt{k_x^2 + k_y^2 + n^2}$ | an overall picture |
| mean profile `meanz` | Chebyshev degree of the horizontal mean | the mean boundary layers |

The horizontal mean ($k_x = k_y = 0$) is left out of the first four and given its own
spectrum, because it is dominated by the mean temperature profile, which is not a
fluctuation.

## Decades

The **decades** $D$ of a spectrum are $\log_{10}$ of its peak minus the mean of
$\log_{10}E$ over its last three modes ({func}`rbc.spectra.core.decades`). Using the end of
the spectrum, rather than its minimum, matters when energy piles up at the grid scale:
the minimum then lies before the cutoff, and "maximum minus minimum" would report more
decades for a worse-resolved run. `pileup` measures how far the end sits above the
spectrum's lowest stretch; above about 0.1 decades it signals under-resolution.

`D_<field>_<panel>` in `results.csv` is the decades of the time-averaged spectrum, the
conservative number, since the time average of the tail is dominated by the moments when
the small scales are most energetic. `Dmin_*` is the smallest value of any single write,
which shows whether a shortfall is steady or comes in bursts.

### D decades, D/2 significant figures

Energy is the square of amplitude, so $D$ decades of energy is $D/2$ decades of amplitude.
The truncation error of a quantity that depends on the field is of the order of the
neglected amplitude relative to the largest, $10^{-D/2}$, so Nu has about $D/2$
significant figures. This rule of thumb came from Keaton Burns (Dedalus).

It is an estimate, not a bound: Nu is an integral dominated by the large scales and the
boundary layers, and often converges faster than the rule suggests. The resolution study
checks it directly: {mod}`rbc.plots.resolution` `decades` plots the actual Nu error of each
run against its $D$ with the line $10^{-D/2}$.

## Screening, not a pass mark

Decades are a **screening** diagnostic: they say where truncation is most likely to
matter and which direction limits a run, but they do not measure the error in Nu. Only
comparing runs that differ in resolution (or timestep) does that, and that comparison is
what the error bars on Nu are built from ({ref}`resolution-error-budget`). A run at 3.95
decades is not meaningfully worse than one at 4.00: by the rule of thumb the difference
is $10^{-1.975}$ against $10^{-2}$, 1.06% against 1.00%.

{mod}`rbc.analysis.resolution` runs the screening tests of one run and writes
`outputs/resolution.txt`, `resolution.pdf` and `resolution.json`:

Decades per direction
: The x, y and z marginals of velocity and temperature are tested separately (the
  horizontal one is not averaged over x and y), so the direction that limits the run is
  named (`D_run_at`). Each $D$ comes with a standard error from a block bootstrap over
  time, the 5th percentile of the single writes, and the pile-up. A direction is `pass`
  if $D$ reaches the target (4 by default, `--target`), `within 2 SE` if $D + 2\,\mathrm{SE}$
  does, `below` otherwise. `D_run` is the smallest over the six. The mean-profile spectrum
  (`meanz`) is reported with a status but does not enter `D_run`.

Derivative-weighted decades `Dk2`
: The decades of $n^2E(n)$. $\langle|\nabla T|^2\rangle$ and $\langle|\omega|^2\rangle$ are
  sums of squared derivatives, whose truncation follows $k^2E$ rather than $E$; this is
  why those identities converge after $\langle wT\rangle$. Exact in $x$ and $y$; a proxy in
  $z$ (with $\zeta=\cos\theta$, $n^2E(n)$ is the spectrum of $\sin\theta\,\partial_\zeta f$).

Dissipation scales
: $\eta = (Pr^2/(Ra(Nu-1)))^{1/4}$ from the exact global dissipation, and the thermal
  scale $\eta\,Pr^{-3/4}$ for $Pr<1$ (Batchelor $\eta\,Pr^{-1/2}$ for $Pr>1$), against the
  largest spacing in each direction: $L/N$ horizontally, the mid-plane spacing of the
  Gauss–Chebyshev grid ($\approx \pi H/(2N_z)$) vertically. $k_\text{max}\eta = \pi\eta/\Delta
  \ge 1$ is the criterion of Grötzbach (1983), $\ge 1.5$ the stricter spectral-DNS one.
  These are volume averages; the local scale in plumes and boundary layers is smaller,
  which is why the spectra, not $\eta$, are the primary test.

Boundary-layer nodes (vertical only)
: The number of Chebyshev nodes inside the thermal and viscous layers. This says
  nothing about the horizontal directions: they are Fourier, with uniform spacing, and a
  Fourier mode is global, so there is no "node inside the layer" to count. Horizontal
  resolution is judged from the x and y spectra and $\Delta x/\eta$. Even vertically,
  the decades of the mean-profile spectrum test the layer more directly than a node count,
  so the counts are reported, not tested.

Exact identities
: The spread between Nu from $\langle wT\rangle$ and from $\langle|\nabla T|^2\rangle$. In a
  dealiased Galerkin code these identities hold for the truncated system too, so their
  agreement is weak evidence of resolution: an under-resolved run can satisfy them and
  still be wrong.

`screen` in `results.csv` is the status of `D_run`, or `pile-up` if any direction's
pile-up exceeds 0.1 decades.

### Why the marginals and not the 3D shells

The marginals are the right test, and a requirement on the 3D shell spectrum would add
nothing to it:

- **Each marginal measures one truncation.** The vertical marginal's last modes are
  exactly the energy at the highest Chebyshev degrees kept, summed over everything else;
  if it has fallen by 4 decades, the vertical truncation error is small whatever happens
  horizontally. The same holds for the horizontal marginal. Together they bound the
  truncation error in every direction.
- **The 3D shell radius is not a physical scale.** It adds a Fourier wavenumber to a
  Chebyshev degree. A Fourier mode $k$ is a wavelength $L_x/k$ everywhere; Chebyshev degree
  $n$ resolves scales of about $1/n$ in the middle of the layer but about $1/n^2$ near the
  plates. Equal radius therefore mixes very different physical scales, and the tail of the
  shell spectrum is not the energy at "the smallest scale".
- **Shells dilute a slow direction.** At a given radius the shell sums many modes. If one
  direction is under-resolved, its highest modes share their shells with many
  well-resolved modes from the other directions, so the shell spectrum can still look
  steep while one marginal is flat.
- **The last shells are incomplete.** Shells beyond the smallest cutoff along any axis
  contain only some of their modes and drop off artificially, so the spectrum is cut at
  the last complete shell, which is also where the axis cutoffs are.

The 3D shells have more points, so they look smoother, but they carry no information
about truncation that the marginals do not. They remain useful as a picture of the
cascade.

The horizontal mean ($k_x=k_y=0$) is left out of the marginals, and it is where the mean
thermal boundary layer lives. Its spectrum is reported separately (`D_<field>_meanz` in
`resolution.txt` and `results.csv`).

## Resolution studies

At selected points (`res_check` in `cases.csv`, stars on the parameter map) the same flow
is run at several $N_z$. {mod}`rbc.plots.resolution`:

`study`
: Nu from each identity relative to the finest run, and the decades, against $N_z$. Nu
  from $\langle wT\rangle$ usually converges first; the vorticity identity last.

`decades`
: Spectra of the runs nearest 4, 5 and 6 decades, and the Nu error against $D$. It also
  fits $D$ linearly against $N_z$ and prints the $N_z$ needed for each target.

```{dropdown} Why D is linear in N_z
:class-container: derivation

For a field analytic in a neighbourhood of the domain, Chebyshev coefficients decay
geometrically, $|c_n| \sim e^{-an}$, so $\log_{10}E$ at the cutoff falls like
$-2aN_z/\ln 10$ while the peak stays put: $D \approx D_0 + (2a/\ln10)\,N_z$. Fourier
coefficients of a smooth periodic field behave the same way. The fit therefore
extrapolates to the resolution needed for a target, as long as the runs are already in
this exponential regime.
```

(resolution-error-budget)=
## Error bars on Nu

Nu has two kinds of error, kept separate:

- **statistical**: the standard error of the time mean (`Nu_se`, {doc}`convergence`).
  It is random and shrinks with longer averaging;
- **systematic**: the bias from truncation and timestepping. It is not random: rerunning
  longer does not reduce it, and it may have the same sign at neighbouring points.

{mod}`rbc.analysis.variants` measures the systematic part by comparing the runs of one
(Ra, Pr) folder, whatever they are: different resolutions, timesteppers or CFL safety
factors ({doc}`../user-guide/directories` has the naming). For a pair,

$$
\delta = Nu_b - Nu_a, \qquad \sigma_\delta = \sqrt{SE_a^2 + SE_b^2}, \qquad
B = |\delta| + 2\sigma_\delta,
$$

and $B$ is a 95% upper bound on the bias of the coarser run. A pair agrees when
$|\delta| \le \max(2\sigma_\delta,\ 0.5\%\,(Nu-1))$ (`--tol`).

Resolution pairs
: Runs with the same timestepper and safety factor. Each is compared with the finest
  of its group, which takes the bound from its comparison with the next coarser run. For
  the finest run this is conservative: with exponential convergence its own error is a
  small fraction of that difference.

Timestep pairs
: Runs on the same grid with a different timestepper or safety factor, compared with
  the baseline (the label without a suffix). A timestep variant inherits the baseline's
  resolution bound; a run without a timestep pair inherits the largest timestep bound
  measured with its stepper and safety factor on the same or a coarser grid (with a
  CFL-limited step, dt shrinks as the grid is refined).

The bounds add linearly, $B_\text{sys} = B_\text{res} + B_{dt}$, because they are biases,
not independent random errors. Note the floor: when two runs agree, $B$ is still about
$2\sqrt2\,SE$, so a tight systematic bound needs both runs averaged long enough.

Points without any comparison get an **estimated** bound ({func}`rbc.analysis.variants.fill_bounds`):
the largest measured relative bound, scaled by $10^{(D_\text{ref}-D_\text{run})/2}$. This
leans on the rule of thumb above, so `Nu_sys_source` in `results.csv` records whether
each bound was `measured` or `estimated`.

### Propagation to exponents

Every exponent estimate (two-point secants, cell centres, row fits, surface fits, slopes
of the compensated plots) is linear in $L = \log_{10}Q$, $\hat\gamma = \sum_i g_iL_i$. A bias
of at most $b_i = B_i/(Q_i\ln10)$ in each $L_i$ then moves it by at most

$$
\Delta\gamma \le \sum_i |g_i|\,b_i,
$$

whatever the signs and correlations of the biases. That is the systematic error bar on
exponents ({mod}`rbc.analysis.exponents`). It is a worst case: biases of one sign at every
point leave a slope unchanged, and the bound assumes the opposite. For a two-point secant
over $h$ decades between points with relative bound $r = B/Q$ it is $2r/(h\ln 10)$: with
$r = 1\%$ and $h = 1/2$, 0.017. Short steps make large systematic bars. Figures draw the statistical error
as a capped bar and add the systematic bound as a thin whisker; `--errorbars stat` or
`none` turns the whiskers or all bars off ({doc}`../user-guide/figures`).

### What not to do about a small shortfall

Filters, hyperviscosity or a larger dealiasing factor all change the equations being
solved or only hide the tail; none reduces the error of Nu. Changing the timestepper or
CFL safety factor changes only the time error, which the timestep pairs measure. If a run
is short of decades, the remedy is resolution, or a resolution pair that shows it does not
matter.

## Spectra at a glance

- `outputs/cumulative_spectra_velocity.pdf` and `..._temperature.pdf`: the four panels of
  one run, time-averaged, with the 5–95% range over time;
- `outputs/decades_timeseries.pdf`: $D(t)$ of every write;
- `outputs/spec_range.txt`: the table of all decades;
- `outputs/resolution.txt`, `resolution.pdf`: the screening tests of one run;
- `<Ra>/<Pr>/outputs/variants.txt`, `variants.pdf`: every run of one point compared, with
  the bounds;
- `<Ra>/<Pr>/outputs/spectra_compare.pdf`: every resolution of one point overlaid. The
  curves should lie on top of each other at low mode numbers and peel off one by one at
  their cutoffs; if they disagree at low mode numbers, those modes are under-resolved in the
  coarser runs.
