# Statistical convergence

A time-averaged Nu has three kinds of error: **bias from the initial transient** if the
averaging window starts too early, **sampling error** because the window is finite, and
**discretization error** from the grid ({doc}`resolution`). This page is about the first
two: how to tell that a window is long enough and starts late enough, and what the rule
for calling a run converged should be.

{mod}`rbc.analysis.nusselt` computes every quantity below and writes them to `info.txt`
and `results.csv`.

## Sampling error of a time average

Let $x(t)$ be a statistically steady signal with variance $\sigma^2$ and autocorrelation
function $\rho(s)$. Its mean over a window of length $T$ has variance

$$
\operatorname{Var}(\bar x) \approx \frac{2\tau_\text{int}\,\sigma^2}{T},\qquad
\tau_\text{int} = \int_0^\infty \rho(s)\,ds,
$$

for $T \gg \tau_\text{int}$. So the window is worth $n_\text{eff} = T/(2\tau_\text{int})$
independent samples, and the {term}`standard error` is $\sigma/\sqrt{n_\text{eff}}$. Saving
more often does not help once the spacing is below $\tau_\text{int}$; only a longer window
does.

```{dropdown} Derivation
:class-container: derivation

$\operatorname{Var}(\bar x) = \frac{1}{T^2}\int_0^T\!\!\int_0^T \sigma^2\rho(t - t')\,dt\,dt'
= \frac{2\sigma^2}{T}\int_0^T\big(1 - \tfrac{s}{T}\big)\rho(s)\,ds$, by changing variables
to $s = t - t'$. For $T \gg \tau_\text{int}$ the factor $1 - s/T$ is 1 wherever $\rho$ is not
negligible, which leaves $2\sigma^2\tau_\text{int}/T$.
```

Two estimates of the standard error are computed. **Batch means is the one every test
uses**; the autocorrelation estimate sets the block length and is reported for comparison.

**Batch means** (`Nu_se`, `se` in `nusselt.json`). Cut the window into $n$ equal blocks of
length $L = T/n$ and take the standard deviation of the block means divided by $\sqrt n$.
If $L \gg \tau_\text{int}$ the block means are nearly independent, each with variance
$S(0)/L$, where $S(0) = 2\tau_\text{int}\sigma^2$ is the spectral density at zero
frequency, so this estimates $\sqrt{S(0)/T}$, the standard error of the whole mean. The
blocks are made at least $10\,\tau_\text{int}$ long, with between 10 and 20 of them
(`BLOCK_TAU`, `N_BLOCKS_MIN`, `N_BLOCKS_MAX`); the estimate is then itself uncertain by
$1/\sqrt{2(n-1)}$, 16% with 20 blocks and 24% with 10, which `info.txt` prints.

**Autocorrelation** (`Nu_se_tau`). Estimate $\rho$ from the series and sum it up to the
smallest window $W$ with $W \ge 5\,\tau_\text{int}(W)$, the automatic windowing of
[Sokal (1997)](#refs-convergence): stopping early would miss correlation, summing the noisy
tail would add variance. This gives $\tau_\text{int}$ and $n_\text{eff}$, which set the
block length and the record-length test. On a test signal with a known autocorrelation time
of 0.5 and record length 2000, the code returns $\tau_\text{int} = 0.495$ and the exact
standard error to three digits.

The same reasoning for turbulence statistics is set out in detail by
[Oliver et al. (2014)](#refs-convergence).

(why-batch-means)=
### Why not the autocorrelation estimate

For $\langle wT\rangle$ the autocorrelation estimate can be far too large. The kinetic
energy budget (free-fall units, no-slip plates, periodic sides) gives, at every instant,

$$
\frac{dK}{dt} = \langle wT\rangle - \sqrt{\frac{Pr}{Ra}}\,\langle|\boldsymbol\omega|^2\rangle,
\qquad K = \tfrac12\langle|\boldsymbol u|^2\rangle,
\qquad\text{so}\qquad
Nu_{wT}(t) = Nu_\omega(t) + \sqrt{Ra\,Pr}\;\frac{dK}{dt}.
$$

The prefactor is large ($\sqrt{Ra\,Pr} \approx 3200$ at $Ra = 10^8$, $Pr = 0.1$) and $K$ is
small, so $Nu_{wT}(t)$ fluctuates much more than the other identities, but the extra term is
a time derivative: over a window it averages to $\sqrt{Ra\,Pr}\,\Delta K/T$, which falls
like $1/T$, not $1/\sqrt T$. In spectral terms its density $(2\pi f)^2 S_K(f)$ vanishes at
$f = 0$, so it adds nothing to $S(0)$; equivalently its autocorrelation integrates to zero,
through negative lobes at longer lags that cancel the positive peak. Sokal's window stops
summing about $5\,\tau_\text{int}$ in, includes the peak, cuts off the cancellation, and so
overestimates $\tau_\text{int}$ and the standard error. A block mean contains only
$\sqrt{Ra\,Pr}\,\Delta K/L$ of it, which is small for long blocks.

```{dropdown} Evidence and a test
:class-container: derivation

All the identities, and the plate flux, differ from one another by time derivatives of
bounded quantities ($K$ above, and stored heat for the others), so they share $S(0)$ and
should have the same standard error. In a run at $Ra = 10^8$, $Pr = 0.1$ over 503 time
units, batch means (10 blocks) gave 0.064 for both $\langle wT\rangle$ and
$\langle|\nabla T|^2\rangle$; the autocorrelation estimate gave 0.061 for
$\langle|\nabla T|^2\rangle$ but 0.146 for $\langle wT\rangle$, and the two means agreed to
0.003, far closer than an independent error of 0.13 would allow. The running $\tau(W)$ in
`info.pdf` shows the mechanism directly: for a series with a derivative term it rises to a
peak and falls back as the negative lobes enter; for a plain correlated series it levels off.

On synthetic series (a correlated base plus a large derivative term, 500 time units, 150
realizations), the true scatter of the window mean was 0.056; batch means with 30 blocks gave
a median of 0.071, and the autocorrelation estimate 0.101. Without the derivative term,
batch means was within 3% of the truth for 10 to 30 blocks; with it, it overestimated by 7%,
16% and 22% for 10, 20 and 30 blocks (the $\Delta K/L$ term), while its own scatter fell from
26% to 15%. Twenty blocks is the compromise. The overestimate is conservative and shrinks as
$1/L^2$.
```

## How precise does Nu need to be?

The precision should be set by what the numbers are used for: the exponents. With
$L = \log_{10}(Nu - 1)$ and a relative standard error $r$ of $Nu - 1$,
$\sigma_L = r/\ln 10$. Propagating through the estimators of {doc}`exponents` (runs one
decade apart in Ra and half a decade apart in Pr):

| Estimate | Standard error | $r = 0.5\%$ | $r = 1\%$ |
|---|---|---|---|
| $\sigma_L$ | $r/\ln10$ | 0.0022 | 0.0043 |
| β from one two-point secant | $\sqrt2\,\sigma_L$ | 0.0031 | 0.0061 |
| α from one two-point secant | $\sqrt2\,\sigma_L/0.5$ | 0.0061 | 0.012 |
| β at a cell centre | $\sigma_L$ | 0.0022 | 0.0043 |
| α at a cell centre | $2\sigma_L$ | 0.0043 | 0.0087 |
| slope of a 4-point slice fit | $\sigma_L/\sqrt5$ | 0.0010 | 0.0019 |

The effects being measured are of order 0.05 (the gap between a measured slice exponent and
1/4) and 0.1 per decade (drifts of α). A relative standard error of **0.5% of Nu − 1** keeps
every estimate in the table at least ten times smaller than those effects; 1% is still
adequate for the slice slopes but leaves single-cell values of α noisy.

```{dropdown} Derivation of the table
:class-container: derivation

A secant over a step $h$ is $(L_2 - L_1)/h$, with variance $2\sigma_L^2/h^2$; $h = 1$ for
β and $h = 1/2$ for α. A cell estimate averages two secants with four distinct corners:
$\beta = (L_{10} - L_{00} + L_{11} - L_{01})/(2h)$ has variance $4\sigma_L^2/(2h)^2 =
\sigma_L^2/h^2$, and likewise α with $k = 1/2$. A least-squares slope through $n$ points
at $x_i$ has variance $\sigma_L^2/\sum(x_i - \bar x)^2$, and for $x = 0, 1, 2, 3$ the sum is 5.
```

## The original rule: thirds within 1%

The group's rule so far: split the window into thirds; the run is converged if, for
$\langle wT\rangle$ and $\langle|\nabla T|^2\rangle$, the largest difference of $Nu - 1$
between the thirds is under 1% (`converged`, `section_spread_pct`).

(convergence-identities)=
:::{note}
Only those two identities are checked, in both rules. The vorticity form
$1 + Pr\,\langle|\boldsymbol\omega|^2\rangle$ is still computed and reported, but it
weights the smallest scales, and a run can be under-resolved in vorticity while the
velocity and temperature spectra fall off far enough that the heat transport is not
affected; that is one of the points of this project. Its difference from
$\langle wT\rangle$, `vort_diff_pct`, is a resolution diagnostic
({doc}`resolution`), not a convergence test.
:::

It answers the right question (would a different stretch of the run give a different
answer?) but with a fixed threshold rather than one tied to the noise, and with only three
values, so the outcome is largely chance. Each third has standard error $\sqrt3$ times the
whole window's, and the expected range of three independent normal values is
$3/\sqrt\pi \approx 1.69$ standard deviations, so

$$
\mathbb E[\text{spread}] \approx 1.69\times\sqrt3\times r \approx 2.9\,r .
$$

The 1% threshold therefore corresponds on average to $r \approx 0.34\%$, but any one run
can land far from average:

| Relative standard error $r$ of the full window | 0.2% | 0.25% | 0.3% | 0.4% | 0.5% |
|---|---|---|---|---|---|
| Probability that the thirds differ by more than 1% | 10% | 23% | 36% | 56% | 69% |

```{dropdown} How the probabilities were computed
:class-container: derivation

For $n$ independent standard normal values the range $R$ has
$P(R \le w) = n\int\phi(x)\,[\Phi(x + w) - \Phi(x)]^{n-1}dx$. With $n = 3$ and
$w = 1\% / (\sqrt3\,r)$, the threshold in units of one third's standard error, the integral
was evaluated numerically (SciPy `quad`). This assumes the three section means are
independent, which holds when each third is much longer than $\tau_\text{int}$.
```

So a run with an excellent $r = 0.25\%$ fails a quarter of the time, and a run four times
noisier than needed passes a third of the time. That matches experience: runs that "won't
quite agree to 1%" while Nu barely changes with another 200 time units are usually already
precise; the spread is about three standard errors and fluctuates. Using quarters, or
0.75%, moves the threshold but keeps both problems.

### Running averages from different start times

The proposal to compute the running average from several start times, and to call the run
converged when it changes by less than 0.1% over the last 100 time units, contains two
separate ideas.

*Varying the start time* tests for the transient: if dropping more of the beginning
changes the mean by more than its standard error, the transient is not over. This is now
drawn in `info.pdf` (mean over $[t_0, t_1]$ against $t_0$, with a $\pm2$ standard-error
band), and it is the visual counterpart of the MSER rule below.

*A small change over the last 100 units* is weaker than it looks. Adding a stretch
$\Delta T$ to a window $T$ moves the running mean by
$(\Delta T/T)(\bar x_\text{new} - \bar x)$, which shrinks like $1/T$ however poorly
converged the mean is: after 1000 time units it says only that the last 100 units agree
with the rest to about 1%, a single noisy block. It confirms that the mean has stopped
moving, not how accurate it is. Comparing the change with the standard error keeps the
idea and fixes the scale.

## Proposed rule

`converged_stat` is true when, for $\langle wT\rangle$ and $\langle|\nabla T|^2\rangle$:

| Test | Requirement | Purpose |
|---|---|---|
| precision | batch-means standard error ≤ 0.5% of $Nu - 1$ | precise enough for the exponents (table above) |
| record length | $n_\text{eff} \ge 50$, i.e. $T \ge 100\,\tau_\text{int}$ | at least 10 blocks of $10\,\tau_\text{int}$, so the standard error is itself trustworthy |
| sections | the thirds agree within their errors: $p \ge 0.01$ in a $\chi^2$ test | no step change within the window |
| trend | linear trend of the block means below 3 standard errors | no slow drift |

All four use the batch-means standard error. The thresholds are in
`rbc.analysis.nusselt.CRITERIA`. The plate imbalance $\overline{Nu_\text{bot} -
Nu_\text{top}}$ is still computed and printed (in standard errors), but it is not part of
the rule.

```{dropdown} Why 50 effective samples
:class-container: derivation

The windowed estimate of $\tau_\text{int}$ has relative variance about $2(2W+1)/N$ for a
window of $W$ samples out of $N$ (Sokal 1997). With $W = 5\tau_\text{int}$ in samples this
is $20\,\tau_\text{int}/T$, which is 0.2 at $T = 100\,\tau_\text{int}$: $\tau_\text{int}$ is
then known to about 45%, good enough to set the block length. And $T \ge 100\,\tau_\text{int}$
is exactly what lets the batch-means estimate have at least 10 blocks of at least
$10\,\tau_\text{int}$, so that the block means are nearly independent and the standard error
is itself known to about 24% or better.
```


```{admonition} Status
:class: important

This rule is computed for every run but does not yet replace the original one:
`converged` (thirds within 1%) still decides which run is the production run. Switch
once it has been checked against existing runs and agreed with David.
```

(transient)=
## When is the transient over?

Three checks, all automatic (the third is reported only):

**MSER-5** ([White 1997](#refs-convergence)). Average the series in batches of five
samples, and choose the start $d$ that minimizes
$\sum_{i \ge d}(x_i - \bar x_d)^2/(n - d)^2$, the squared standard error of the
remaining mean: discarding transient data lowers it, discarding steady data raises it.
`t0_mser` is the later of the two checked identities. MSER is a minimum: it stops removing
data once the remaining bias is below the noise, so the window start can sensibly be later.
`AVG_TIME=auto` uses it.

**Start-time scan.** The mean against the window start, in `info.pdf`; it should be flat
within its band from the chosen start on.

**Plate balance.** While the mean temperature is still adjusting, the bottom and top plates
carry different fluxes ({doc}`nusselt`). The imbalance is printed in `info.txt` for
inspection; it is not part of `converged_stat`.

The run's own signals add a further, visual check: the kinetic energy and Nu should
fluctuate about a level, with no long climb or decay, before the window starts.

### A minimum averaging time?

A fixed number of free-fall times is a poor universal rule, because the correlation time
varies across parameter space; $n_\text{eff} \ge 50$ is the same requirement measured in
the right unit. In practice no run has met the precision target in under about 150
free-fall times, which is a useful planning figure for how long to ask for.

## Recommended practice

1. Set `AVG_TIME=auto` (or a time chosen from the start-time scan) and run until
   `converged_stat` is true; `info.txt` says which test fails.
2. If only the precision test fails, run longer: the standard error falls like
   $1/\sqrt T$, so halving it takes four times as long.
3. If the section or trend test fails with good precision, look at the series: something
   changed during the run (a resolution change, a restart from another run, a slow
   large-scale reorganization).
4. Quote $Nu \pm$ standard error, never a spread.

(refs-convergence)=
## References

- A. D. Sokal (1997), Monte Carlo methods in statistical mechanics: foundations and new
  algorithms, in *Functional Integration*, NATO ASI Series B 361, Springer.
- K. P. White (1997), An effective truncation heuristic for bias reduction in simulation
  output, *Simulation* 69(6), 323–334.
- T. A. Oliver, N. Malaya, R. Ulerich and R. D. Moser (2014), Estimating uncertainties in
  statistics computed from direct numerical simulation, *Physics of Fluids* 26, 035101
  ([arXiv:1311.0828](https://arxiv.org/abs/1311.0828)).
