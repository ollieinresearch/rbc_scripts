# Placing slices on a compensated plot

The compensated plots ({mod}`rbc.plots.scaling`, figures A and B) divide $Q$ along each
slice $Pr \propto Ra^{a}$ by the power of Ra a reference law predicts there,
$Ra^{\beta + a\alpha}$, so that the law is a horizontal line and the fitted slope $m$ of
each slice is its departure from the law. With $Q \sim Pr^{\alpha}Ra^{\beta}$ and the
Shishkina–Lohse values $\alpha = \beta = 1/2$, that is $Nu/Ra^{1/4}$ along
$Pr \propto Ra^{-1/2}$.

Where each slice sits vertically is a separate choice, `--norm`.

## What normalization can and cannot change

Dividing a slice by a constant shifts it on a logarithmic axis without changing its shape.
So **every option gives the same slopes**; the fits and the conclusions about exponents do
not depend on it. The choice only decides what a reader compares by eye: slopes, or
offsets between slices.

## The options

`none`
: The compensated values as they are. Offsets between slices then carry the Pr dependence
  that the reference law does not remove.

`first`
: Each slice divided by its first (lowest-Ra) point, so all start at 1. Emphasizes how
  far each slice drifts, but the lowest-Ra point is often the least asymptotic one, and
  its error moves the whole curve.

`ref`
: Each slice divided by its value at one Ra, `--ref-ra` (default $10^7$): all slices
  pass through 1 there.

`common`
: Each slice divided by its geometric mean over the Ra values that every slice shares.
  This is the least-squares vertical alignment of the slices on their common range: no
  single point is privileged, and the slices pivot about the middle of that range.

`law`
: No per-slice constant at all: $Q/(Pr^\alpha Ra^\beta)$ for every point. If the
  reference law held, every slice would collapse onto one horizontal line. Offsets between
  slices then test the Pr exponent, and the slopes test the slice exponent.

## Is normalizing at $Ra = 10^7$ a good idea?

The suggestion to make all slices equal at $Ra = 10^7$ (`ref`) is less arbitrary than it
looks: $10^7$ is the highest Ra that every slice reaches, so the slices are aligned where
they are closest to their asymptotic behaviour, and their different slopes fan out
towards lower Ra, where a reader sees them most easily. Its weaknesses are that a single
run fixes each whole curve, so that run's sampling error shifts all of it, and that the
common crossing point at $10^7$ suggests a convergence that is only a choice of
normalization.

Recommended:

- **`law` as the main figure.** It has no free constants, and the collapse (or not) is
  itself the test of the theory: every departure visible on it is a departure from the
  law.
- **`common` as the companion** if a pure comparison of slopes is wanted.
- If the reference-point version is kept, justify it in the caption as "normalized at the
  highest Ra common to all slices", and prefer the fitted value there to the data point,
  which is what `ref` does when the point is missing. Doing it always would remove the
  dependence on one run.

```bash
python -m rbc.plots.scaling --q Nu --norm law
python -m rbc.plots.scaling --q Nu --norm common
python -m rbc.plots.scaling --q Nu --norm ref --ref-ra 7
```
