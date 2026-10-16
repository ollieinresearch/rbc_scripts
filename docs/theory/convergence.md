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

Two estimates of the standard error are computed:

**Batch means** (`Nu_se`). Split the window into 10 equal blocks; the standard deviation of
the block means divided by $\sqrt{10}$. Simple and assumption-free, provided each block is
much longer than $\tau_\text{int}$. With only 10 blocks the estimate itself is uncertain by
about $1/\sqrt{2\cdot 9} \approx 24\%$.

**Autocorrelation** (`Nu_se_tau`). Estimate $\rho$ from the series and sum it up to the
smallest window $W$ with $W \ge 5\,\tau_\text{int}(W)$, the automatic windowing of
[Sokal (1997)](#refs-convergence): stopping early would miss correlation, summing the noisy
tail would add variance. This gives $\tau_\text{int}$ and $n_\text{eff}$ as well. On a test
signal with a known autocorrelation time of 0.5 and record length 2000, the code returns
$\tau_\text{int} = 0.495$ and the exact standard error to three digits.

The two agree when the window is long; when they disagree, the window is too short for
either to be trusted. The same reasoning for turbulence statistics is set out in detail by
[Oliver et al. (2014)](#refs-convergence).

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

The group's rule so far: split the window into thirds; the run is converged if, for all
three identities, the largest difference of $Nu - 1$ between the thirds is under 1%
(`converged`, `section_spread_pct`).

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

`converged_stat` is true when, for each of the three identities:

| Test | Requirement | Purpose |
|---|---|---|
| precision | standard error ≤ 0.5% of $Nu - 1$ | precise enough for the exponents (table above) |
| record length | $n_\text{eff} \ge 50$, i.e. $T \ge 100\,\tau_\text{int}$ | the standard error itself is trustworthy (to about 20%) |
| sections | the thirds agree within their errors: $p \ge 0.01$ in a $\chi^2$ test | no step change within the window |
| trend | linear trend over the window below 3 standard errors | no slow drift |

and the plate fluxes balance: $\overline{Nu_\text{bot} - Nu_\text{top}}$ below 3 standard
errors. The thresholds are in `rbc.analysis.nusselt.CRITERIA`.

```{dropdown} Why 50 effective samples
:class-container: derivation

The windowed estimate of $\tau_\text{int}$ has relative variance about $2(2W+1)/N$ for a
window of $W$ samples out of $N$ (Sokal 1997). With $W = 5\tau_\text{int}$ in samples this
is $20\,\tau_\text{int}/T$, which is 0.2 at $T = 100\,\tau_\text{int}$: $\tau_\text{int}$ is
then known to about 45%, and the standard error, proportional to $\sqrt{\tau_\text{int}}$,
to about 22%.
```


```{admonition} Status
:class: important

This rule is computed for every run but does not yet replace the original one:
`converged` (thirds within 1%) still decides which run is the production run. Switch
once it has been checked against existing runs and agreed with David.
```

(transient)=
## When is the transient over?

Three checks, all automatic:

**MSER-5** ([White 1997](#refs-convergence)). Average the series in batches of five
samples, and choose the start $d$ that minimizes
$\sum_{i \ge d}(x_i - \bar x_d)^2/(n - d)^2$, the squared standard error of the
remaining mean: discarding transient data lowers it, discarding steady data raises it.
`t0_mser` is the latest over the three identities. MSER is a minimum: it stops removing
data once the remaining bias is below the noise, so the window start can sensibly be later.
`AVG_TIME=auto` uses it.

**Start-time scan.** The mean against the window start, in `info.pdf`; it should be flat
within its band from the chosen start on.

**Plate balance.** While the mean temperature is still adjusting, the bottom and top plates
carry different fluxes ({doc}`nusselt`); the imbalance must be insignificant.

The run's own signals add a fourth, visual check: the kinetic energy and Nu should
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
