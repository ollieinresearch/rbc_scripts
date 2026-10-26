# Paper figures

Every figure is one module of {mod}`rbc.plots`, reading `data/results.csv` (and
`data/cases.csv`) and writing a PDF to `figures/` (and a PNG with `--png`).

```bash
python -m rbc.plots.all --png        # all of them
python -m rbc.plots.all --results nusselt.csv   # Nu figures from a simple Ra,Pr,Nusselt file
```

Every figure takes `--results`, which may also be a simple three-column file
({ref}`simple-results`).

Figures are made at the size they will be printed, for inclusion at 100% in LaTeX: widths
are fractions of `rbc.style.TEXTWIDTH_IN`, the one number to match to the journal's text
width. Text is typeset with LaTeX where it is installed (`--no-tex` turns this off).
Colours and markers of slices are fixed in {mod}`rbc.style`, so a slice looks the same
in every figure. Each figure below links to its example in the {doc}`../gallery/index`.

## Parameter map

{mod}`rbc.plots.param_grid` ({ref}`example <sphx_glr_gallery_auto_examples_plot_param_grid.py>`)

Every point of `cases.csv`: filled circles finished, open circles running, stars
resolution studies, and the slices as coloured lines. `--ra-inv` adds the
$Pr \propto Ra^{-1}$ slices (dashed); `--span finished` stops each line at its last
finished point.

## Compensated Nu and Re along slices (figures A and B)

{mod}`rbc.plots.scaling` ({ref}`example <sphx_glr_gallery_auto_examples_plot_scaling.py>`)

Nu (or $Nu-1$, or Re) divided by the power of Ra that the reference law predicts along
each slice, so the prediction is a horizontal line and the fitted slope $m$ is the
departure from it. The reference law $Q \sim Pr^\alpha Ra^\beta$ defaults to the
Shishkina–Lohse values; `--alpha` and `--beta` change it. `--fit quadratic` adds a
drifting-exponent fit; `--family -0.5 -1` draws both families. How the slices are placed
vertically is `--norm`, discussed in {doc}`../theory/normalization`.

```bash
python -m rbc.plots.scaling --q Nu --norm law
python -m rbc.plots.scaling --q Re
```

## Nu − 1 at fixed Ra or fixed Pr (figure D)

{mod}`rbc.plots.fixed` ({ref}`example <sphx_glr_gallery_auto_examples_plot_fixed.py>`)

Top, $Nu-1$ against Pr with one curve per Ra; bottom, against Ra with one curve per Pr,
each with a one-column legend beside it. `--slopes` puts each curve's single fitted
exponent in its legend (a summary only, since the curves bend); `--lit` overlays published
data from a CSV.

On the raw plot every curve rises by several decades, so differences of a few hundredths
in exponent cannot be seen. `--norm` divides them out against a reference law
$Nu - 1 \sim Pr^\alpha Ra^\beta$, by default the one power law fitted to all points
(`--ref sl` for Shishkina–Lohse, or `--alpha`, `--beta`), and `--fit linear|quadratic`
draws a fit through each curve instead of joining the points. Two of these are made by
`rbc.plots.all`:

`--norm comp --fit quadratic`
: Compensated by the plotted parameter: $(Nu-1)/Pr^\alpha$ against Pr, $(Nu-1)/Ra^\beta$
  against Ra. A curve is flat where its exponent equals the reference one, so curvature
  (an exponent that drifts) is easy to see, and the curves stay apart, one per level. The
  legend gives $m$, the departure from the reference exponent at the middle of each
  curve. Curves of fewer than four points get a straight fit.

`--norm collapse --fit linear`
: Divided by the fixed parameter instead: $(Nu-1)/Ra^\beta$ against Pr,
  $(Nu-1)/Pr^\alpha$ against Ra. If $Nu - 1$ is a product $f(Pr)\,g(Ra)$ the curves fall
  on top of one another, on the dashed line $C\,Pr^\alpha$ (or $C\,Ra^\beta$) of the
  reference law. A curve off the bundle is a level that does not follow the others; the
  slopes are the curves' own exponents, unchanged by the division.

The others: `--norm law` divides by the whole law, so every curve would be the dashed
horizontal line (the collapse with its trend removed, magnifying departures but crowding
the curves); `--norm common` divides each curve by its own geometric mean over the values
all curves share, with no reference law at all. Each choice writes its own file, named
after it (`fixed_Num1_pr_ra_comp_fit_quadratic.pdf`), with the slopes in a matching CSV.

```bash
python -m rbc.plots.fixed --slopes --norm comp --fit quadratic
python -m rbc.plots.fixed --slopes --norm collapse --fit linear --ref sl
```

## Local exponents

{mod}`rbc.plots.exponents` ({ref}`example <sphx_glr_gallery_auto_examples_plot_exponents.py>`)

The local exponents α and β of $Nu - 1 \sim Pr^\alpha Ra^\beta$, each against Ra and
against Pr. Two-point estimates from the same row (constant Pr, or constant Ra) are joined,
so each row's trend can be followed; the black line and the panel title give the trend
over all estimates. `--surface 2` prints the comparison of surface fits (`--draw-surface`
also draws one); `--min-q 1` leaves out estimates that use a run with $Nu - 1 < 1$, where
local exponents stop meaning much. Method and interpretation:
{doc}`../theory/exponents`.

## Boundary layers (figure F)

{mod}`rbc.plots.boundary_layers` ({ref}`example <sphx_glr_gallery_auto_examples_plot_boundary_layers.py>`)

`ratio`: $\delta_u/\delta_T$ against Pr with the Prandtl–Blasius asymptotes, and the
number of Chebyshev nodes inside each layer. `grid --points "5 0" "8 -3"`: near-wall
profiles with every node drawn (needs the runs' `profiles.npz`).

## Resolution (figure H)

{mod}`rbc.plots.resolution` ({ref}`example <sphx_glr_gallery_auto_examples_plot_resolution.py>`)

`study`: for every point with several resolutions, Nu from each identity relative to the
finest run, and the spectral decades, against $N_z$. `decades --point "6 -2"`: spectra of
the runs nearest 4, 5 and 6 decades, and the Nu error against decades with the
$10^{-D/2}$ rule. See {doc}`../theory/resolution`.

## Contours over the parameter plane (figure I)

{mod}`rbc.plots.contour` ({ref}`example <sphx_glr_gallery_auto_examples_plot_contour.py>`)

Filled contours of $\log_{10} Nu$ (or Re, or a compensated Nu with `--comp`) over
(Ra, Pr), interpolated on a triangulation of the runs; nothing is drawn outside their
convex hull. `--smooth` uses a smooth cubic interpolant, which can overshoot where runs
are sparse.

## Adding a figure

See {ref}`adding-a-figure`.
