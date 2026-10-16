# Local scaling exponents

Theories of convection predict power laws, $Q \sim Pr^{\alpha}Ra^{\beta}$, for
$Q = Nu - 1$, $Nu$ or $Re$, with exponents that change between regimes. Throughout the
code and the figures, **α is the exponent of Pr and β the exponent of Ra**. A single fit
over all runs gives one pair; the question here is whether α and β are constant over the
parameter range, or drift, as they would through a transition (for example from 1/3
towards 1/2).

{mod}`rbc.analysis.exponents` estimates them; {mod}`rbc.plots.exponents` plots them
({ref}`example <sphx_glr_gallery_auto_examples_plot_exponents.py>`).

## Working in logarithms

Let $x = \log_{10}Ra$, $y = \log_{10}Pr$ and $L = \log_{10}Q$. Then

$$
\alpha = \frac{\partial L}{\partial y},\qquad \beta = \frac{\partial L}{\partial x},
$$

the local slopes of the surface $L(x, y)$. The runs sit on a grid with spacing $h = 1$
(a decade) in $x$ and $k = 1/2$ in $y$. The standard error of $L$ is
$\sigma_L = \sigma_Q/(Q\ln10)$, from the time-averaging error of each run
({doc}`convergence`).

## Two-point estimates at midpoints

Along a row of fixed Pr, the secant between neighbours estimates β at their midpoint;
along a column of fixed Ra, it estimates α:

$$
\beta\big(x + \tfrac h2, y\big) \approx \frac{L(x + h, y) - L(x, y)}{h}.
$$

This is a central difference about the midpoint, accurate to $O(h^2)$, and **exact
whenever $L$ is quadratic**, which is the simplest surface whose exponents drift. Error:
$\sqrt{\sigma_{L,1}^2 + \sigma_{L,2}^2}/h$.

## Both exponents at cell centres

Take a grid cell with corners $L_{00} = L(x, y)$, $L_{10} = L(x+h, y)$,
$L_{01} = L(x, y+k)$, $L_{11} = L(x+h, y+k)$. Averaging the two secants on opposite edges
puts both exponents at the centre:

$$
\beta = \frac{(L_{10} - L_{00}) + (L_{11} - L_{01})}{2h},\qquad
\alpha = \frac{(L_{01} - L_{00}) + (L_{11} - L_{10})}{2k}.
$$

With $h = 1$ and $k = 1/2$ the cell's diagonals are slices. The diagonal from $L_{01}$ to
$L_{10}$ goes one decade up in Ra and half a decade down in Pr: one step along
$Pr \propto Ra^{-1/2}$. The other diagonal is a step along $Pr \propto Ra^{+1/2}$. Their
secants are

$$
s_- = \frac{L_{10} - L_{01}}{h} = \beta - \frac{\alpha}{2},\qquad
s_+ = \frac{L_{11} - L_{00}}{h} = \beta + \frac{\alpha}{2},
$$

so the cell centre lies on a slice, and the slope measured along that slice is
$s_- = \beta - \alpha/2$: exactly the quantity the compensated slice plots measure
({doc}`normalization`). In general, along $Pr \propto Ra^{a}$ the exponent is
$\beta + a\alpha$; the [Shishkina–Lohse](https://arxiv.org/abs/2407.16573) values
$\alpha = \beta = 1/2$ give $1/4$ for $a = -1/2$.

```{dropdown} Derivation: the diagonal secants
:class-container: derivation

$L_{10} - L_{01} = (L_{10} - L_{00}) - (L_{01} - L_{00})$, a step $+h$ in $x$ minus a step
$+k$ in $y$. For locally linear $L$ that is $\beta h - \alpha k$, and with $k = h/2$,
$s_- = \beta - \alpha/2$. The same for $s_+$ with a plus sign. Adding and subtracting:
$\beta = (s_+ + s_-)/2$ and $\alpha = s_+ - s_-$, which reproduces the cell formulas above.
```

Each cell estimate comes with its standard errors and the covariance of α and β (they share
the four corners), so the error of $s_-$ is computed correctly.

## Row and column regressions

`--method lines` fits one weighted straight line of $L$ against $x$ along each complete row
(giving β) and against $y$ along each column (α). Fewer, smoother estimates, each an
average over the row; useful to see whether the two-point values scatter about a common
line.

## Fitting a surface

The most economical description is a polynomial $L(x, y)$ of total degree $n$, fitted to all
runs by weighted least squares, whose derivatives give α and β everywhere with errors from
the fit's covariance:

| Degree | Parameters | Exponents |
|---|---|---|
| 1 | 3 | constant: a single power law |
| 2 | 6 | each changes linearly across the plane: the simplest crossover |
| 3 | 10 | quadratic changes |

Degrees are compared with AICc, Akaike's criterion with the small-sample correction
$2p(p+1)/(n - p - 1)$ for $p$ parameters and $n$ runs; with 15–25 runs the correction
matters (without it, the plain criterion picks an overfitted cubic). `--surface 2` prints
the table and overlays the quadratic surface's exponents, drawn only across the range of
the data.

### Other ways to interpolate exponents over the plane

Ordered from least to most machinery:

1. **Piecewise-linear on triangles.** Split each cell along a diagonal and take each
   triangle's gradient. Model-free; equivalent to the cell estimates.
2. **Polynomial surfaces** (implemented, above).
3. **A two-regime model**, $Nu - 1 = A\,Pr^{\alpha_1}Ra^{\beta_1} + B\,Pr^{\alpha_2}Ra^{\beta_2}$,
   in the spirit of Grossmann–Lohse: fits where a transition happens and its two limiting
   exponents directly. Worth trying if the exponents trend.
4. **Gaussian-process regression** of $L(x, y)$, using each run's standard error as its
   noise: its derivatives come with credible bands, the most principled way to claim a
   transition with stated confidence.

Smoothing splines are not recommended at this density of runs: their derivatives are
noisy.

## Testing for constant exponents

In each panel of the figure, the estimates of one exponent are fitted two ways against the
log of one parameter ({func}`rbc.analysis.exponents.flatness`):

- a constant, giving the weighted mean and $\chi^2_\nu$; $\chi^2_\nu \approx 1$ means the
  scatter is what the errors predict;
- a line, giving the trend per decade with its error.

A trend several standard errors from zero, together with $\chi^2_\nu \gg 1$ for the
constant, is evidence that the exponent changes with that parameter. Because the errors of
neighbouring two-point estimates share runs, they are not independent, so treat
$\chi^2_\nu$ as a guide rather than an exact test.
