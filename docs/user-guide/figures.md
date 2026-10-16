# Paper figures

Every figure is one module of {mod}`rbc.plots`, reading `data/results.csv` (and
`data/cases.csv`) and writing a PDF to `figures/` (and a PNG with `--png`).

```bash
python -m rbc.plots.all --png        # all of them
```

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

Left, $Nu-1$ against Pr with one curve per Ra; right, against Ra with one curve per Pr.
`--slopes` puts each curve's fitted exponent in its legend; `--lit` overlays published
data from a CSV.

## Local exponents

{mod}`rbc.plots.exponents` ({ref}`example <sphx_glr_gallery_auto_examples_plot_exponents.py>`)

The local exponents α and β of $Nu - 1 \sim Pr^\alpha Ra^\beta$, each against Ra and
against Pr, with a test of whether they are constant. `--surface 2` overlays a quadratic
surface fit and prints the model comparison. Method and interpretation:
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
