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

## The resolution rule

A run counts as resolved when **both the horizontal and the vertical marginal spectra, of
both velocity and temperature, span at least 4 decades**. The figure `resolution_study`
uses the smallest of those four, `D_run`.

By the rule of thumb, 4 decades corresponds to an error of about 1% in Nu, larger than the
0.5% sampling precision recommended in {doc}`convergence`. Whether the resolution error
really is that large is exactly what the resolution studies show; if their errors sit
well below $10^{-D/2}$, 4 decades is enough, otherwise the requirement should rise to
about 5.

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

```{admonition} Suggestion
:class: tip

The marginals leave out the horizontal mean, which is where the mean thermal boundary
layer lives, so the rule never looks at it. Its spectrum is in `results.csv`
(`D_temp_meanz`). Since Nu is set by the boundary-layer gradient, requiring 4 decades there
too would close the gap.
```

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

## Spectra at a glance

- `outputs/cumulative_spectra_velocity.pdf` and `..._temperature.pdf`: the four panels of
  one run, time-averaged, with the 5–95% range over time;
- `outputs/decades_timeseries.pdf`: $D(t)$ of every write;
- `outputs/spec_range.txt`: the table of all decades;
- `<Ra>/<Pr>/outputs/spectra_compare.pdf`: every resolution of one point overlaid. The
  curves should lie on top of each other at low mode numbers and peel off one by one at
  their cutoffs; if they disagree at low mode numbers, those modes are under-resolved in the
  coarser runs.
