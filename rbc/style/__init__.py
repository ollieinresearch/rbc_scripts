"""Single source of truth for how every figure looks.

Two Matplotlib style sheets: ``paper`` (publication figures, made at printed size)
and ``diag`` (large per-run diagnostic figures, layered on ``paper``). Colours and
markers of the slices, of curves at one Ra or one Pr, run statuses and resolution stars
are defined here once, so every figure agrees: colours run from red (cheapest runs) to
purple (most expensive), see :data:`EXPENSE_LEVELS`. Usage::

    from rbc import style
    style.use("paper")                       # or "diag" for per-run cluster plots
    fig, ax = plt.subplots(figsize=style.figsize(0.5))
    ax.plot(ra, nu, **style.slice_kw(c))
    style.save(fig, outdir / "nu_slices")    # -> .pdf (+ .png with png=True)
"""
from __future__ import annotations

import shutil
from pathlib import Path

import matplotlib as mpl

_HERE = Path(__file__).resolve().parent

#: Text width of the paper in inches: \showthe\textwidth in the LaTeX source, in pt,
#: divided by 72.27. Placeholder until checked against the journal template.
TEXTWIDTH_IN = 5.3

#: When True, :func:`save` leaves figures open (the documentation gallery needs them).
KEEP_OPEN = False

#: Symbols used in every axis label and legend, so notation changes in one place.
SYM = {
    "Ra": r"$Ra$", "Pr": r"$Pr$", "Nu": r"$Nu$", "Re": r"$Re$",
    "Pr0": r"$Pr_0$", "Nz": r"$N_z$", "z": r"$z$", "t": r"$t$",
}

#: Paul Tol's discrete rainbow (Tol 2021, SRON/EPS/TN/09-002, fig. 22), by his numbering.
#: The colours sit at equal steps along a curve in CIELAB and alternate light and dark, so
#: neighbours differ in lightness as well as hue and stay apart for red-green colour
#: blindness. They are used as given, never interpolated.
TOL_RAINBOW = {9: "#882E72", 10: "#1965B0", 14: "#7BAFDE", 15: "#4EB265", 17: "#CAE0AB", 18: "#F7F056",
               21: "#F4A736", 23: "#EE8026", 24: "#E8601C", 26: "#DC050C", 28: "#72190E"}
#: Tol's choice of colours for n categories (n = 1 to 10), red first.
TOL_SETS = {1: [10], 2: [26, 10], 3: [26, 18, 10], 4: [26, 18, 15, 10], 5: [26, 18, 15, 14, 10],
            6: [26, 18, 17, 15, 14, 10], 7: [26, 18, 17, 15, 14, 10, 9], 8: [26, 23, 18, 17, 15, 14, 10, 9],
            9: [28, 26, 23, 18, 17, 15, 14, 10, 9], 10: [28, 26, 24, 21, 18, 17, 15, 14, 10, 9]}
_HALF = [round(0.5 * k, 1) for k in range(8, -9, -1)]          # 4, 3.5, ..., -4

#: Colours by expense. Every figure that colours runs, curves or slices by a parameter orders
#: them like the spectrum, red (cheapest) through orange, yellow, green and blue to purple
#: (most expensive). Higher Ra, lower Pr and lower slice Pr_0 are more expensive. Each kind
#: of parameter has fixed levels, cheapest first, each with its own colour from
#: :data:`TOL_RAINBOW`, so a value has the same colour in every figure. ``ra``: log10 Ra;
#: ``pr``: log10 Pr; ``slice-0.5`` and ``slice-1``: c = log10 Pr_0 of the slices
#: Pr ∝ Ra^-1/2 and Pr ∝ Ra^-1 (the two families share labels, not runs). A value between
#: levels takes the nearest one; beyond the ends, the end colour. Continuous fields
#: (contours) keep an ordinary colormap.
EXPENSE_LEVELS = {
    "ra": ([5.0, 6.0, 7.0, 8.0, 9.0], [26, 18, 15, 10, 9]),
    "pr": ([v for v in _HALF if -3.5 <= v <= 0.0], TOL_SETS[8]),
    "slice-0.5": ([v for v in _HALF if -3.0 <= v <= 1.0], TOL_SETS[9]),
    "slice-1": ([v for v in _HALF if -1.5 <= v <= 2.0], TOL_SETS[8]),
}
#: The values at the cheap and expensive ends of each kind.
EXPENSE_RANGE = {k: (v[0][0], v[0][-1]) for k, v in EXPENSE_LEVELS.items()}
#: Markers by position on the half-decade grid, so colour is never the only cue.
GRID_MARKERS = "osD^vph*PX<>d"
#: Colours whose contrast with white is below this get a dark outline on lines.
PALE_CONTRAST = 2.0


def _luminance(colour) -> float:
    """Relative luminance (WCAG 2) of a colour."""
    from matplotlib.colors import to_rgb
    c = [x / 12.92 if x <= 0.04045 else ((x + 0.055) / 1.055) ** 2.4 for x in to_rgb(colour)]
    return 0.2126 * c[0] + 0.7152 * c[1] + 0.0722 * c[2]


def is_pale(colour) -> bool:
    """Whether a colour is too light to read as a thin line on white.

    Parameters
    ----------
    colour : colour-like

    Returns
    -------
    bool
        True when its contrast ratio with white is below :data:`PALE_CONTRAST` (Tol's
        yellow and pale green).
    """
    return 1.05 / (_luminance(colour) + 0.05) < PALE_CONTRAST


def edge_colour(colour):
    """A darker shade of a colour, for marker edges and the outline of pale lines.

    Parameters
    ----------
    colour : colour-like

    Returns
    -------
    tuple
        RGB, the colour scaled towards black.
    """
    from matplotlib.colors import to_rgb
    return tuple(0.55 * x for x in to_rgb(colour))


def _outline(colour, extra=0.7):
    """Path effects that draw a line on top of a slightly wider one in its dark shade."""
    from matplotlib import patheffects as pe

    class Outline(pe.AbstractPathEffect):
        """Draws the path ``extra`` points wider, in the dark shade, under the original
        (lines and the rims of open markers); filled markers are skipped, as their edge is
        already dark."""

        def draw_path(self, renderer, gc, tpath, affine, rgbFace=None):
            if rgbFace is not None and (len(rgbFace) < 4 or rgbFace[3] > 0):
                return
            gc0 = renderer.new_gc()
            gc0.copy_properties(gc)
            gc0.set_linewidth(gc.get_linewidth() + extra)
            gc0.set_foreground(edge_colour(colour))
            renderer.draw_path(gc0, tpath, affine, None)
            gc0.restore()

    return [Outline(), pe.Normal()]


def line_kw(colour) -> dict:
    """Keywords for a line of a given colour: the colour, and a dark outline if it is pale.

    Parameters
    ----------
    colour : colour-like

    Returns
    -------
    dict
        ``color`` and, for a pale colour, ``path_effects``.
    """
    return dict(color=colour, path_effects=_outline(colour)) if is_pale(colour) else dict(color=colour)


def marker_kw(colour, filled: bool = True) -> dict:
    """Keywords for markers (and the line through them) of a given colour.

    Filled markers get an edge in the dark shade of their colour, so the pale ones show on
    white and neighbouring colours differ in two shades, not one. Open markers keep their
    colour, with the dark outline when it is pale.

    Parameters
    ----------
    colour : colour-like
    filled : bool, optional

    Returns
    -------
    dict
    """
    kw = line_kw(colour)
    if filled:
        kw.update(mfc=colour, mec=edge_colour(colour))
    else:
        kw.update(mfc="none", mec=colour)
    return kw


def expense_index(value: float, kind: str) -> int:
    """Index of the level nearest a parameter value (0 = cheapest).

    Parameters
    ----------
    value : float
        log10 Ra, log10 Pr or slice label c, according to ``kind``.
    kind : {"ra", "pr", "slice-0.5", "slice-1"}
        Which parameter.

    Returns
    -------
    int
    """
    levels = EXPENSE_LEVELS[kind][0]
    return min(range(len(levels)), key=lambda i: (abs(levels[i] - value), i))


def expense_colour(value: float, kind: str) -> str:
    """Colour of a parameter value (see :data:`EXPENSE_LEVELS`).

    Parameters
    ----------
    value : float
        log10 Ra, log10 Pr or slice label c.
    kind : {"ra", "pr", "slice-0.5", "slice-1"}
        Which parameter.

    Returns
    -------
    str
        Hex colour.
    """
    return TOL_RAINBOW[EXPENSE_LEVELS[kind][1][expense_index(value, kind)]]


def expense_norm(kind: str, values=None):
    """Discrete colormap and normalization for a colour bar of a parameter.

    One band per level, centred on it, so the bar shows exactly the colours of the points.

    Parameters
    ----------
    kind : {"ra", "pr", "slice-0.5", "slice-1"}
        Which parameter.
    values : array-like, optional
        Values in the figure; the bar then covers only the levels from the smallest to the
        largest.

    Returns
    -------
    tuple
        ``(cmap, norm)``: a ListedColormap and a BoundaryNorm over the parameter itself,
        increasing.
    """
    import numpy as np
    from matplotlib.colors import BoundaryNorm, ListedColormap
    levels, ids = EXPENSE_LEVELS[kind]
    pairs = sorted(zip(levels, ids))
    if values is not None and len(values):
        lo, hi = (expense_index(v, kind) for v in (np.min(values), np.max(values)))
        keep = {levels[i] for i in range(min(lo, hi), max(lo, hi) + 1)}
        pairs = [pr for pr in pairs if pr[0] in keep]
    lv = np.array([v for v, _ in pairs])
    step = np.diff(lv).min() if len(lv) > 1 else 1.0
    bounds = np.concatenate([[lv[0] - step / 2], (lv[1:] + lv[:-1]) / 2, [lv[-1] + step / 2]])
    return ListedColormap([TOL_RAINBOW[i] for _, i in pairs]), BoundaryNorm(bounds, len(pairs))


def ranked_colours(n: int) -> list:
    """``n`` colours from red to purple, for things ranked by expense within one figure
    (resolutions of one point, say): Tol's set for ``n``. Beyond 10 the set of 10 repeats.

    Parameters
    ----------
    n : int
        Number of colours.

    Returns
    -------
    list of str
    """
    ids = TOL_SETS[min(max(n, 1), 10)]
    return [TOL_RAINBOW[ids[i % len(ids)]] for i in range(n)]


def category_kw(value: float, kind: str, line: bool = False) -> dict:
    """Plot keywords for one curve labelled by a parameter value (one Ra, or one Pr).

    Parameters
    ----------
    value : float
        log10 Ra or log10 Pr.
    kind : {"ra", "pr"}
        Which parameter.
    line : bool, optional
        Line keywords only (:func:`line_kw`), instead of markers too.

    Returns
    -------
    dict
    """
    col = expense_colour(value, kind)
    if line:
        return line_kw(col)
    k = int(round(abs(value - EXPENSE_RANGE[kind][0]) / 0.5))
    return dict(marker_kw(col), marker=GRID_MARKERS[k % len(GRID_MARKERS)])


#: Markers of the slices Pr = Pr_0 (Ra/10^5)^a, keyed by c = log10(Pr_0); colours are
#: :func:`expense_colour` of c, so a lower Pr_0 (more expensive) is further towards purple.
SLICE_MARKERS = {2.0: "*", 1.5: "p", 1.0: "h", 0.5: "d", 0.0: "o", -0.5: "s", -1.0: "D", -1.5: "^",
                 -2.0: "v", -2.5: "P", -3.0: "X"}
SLICE_STYLE = {c: dict(color=expense_colour(c, "slice-0.5"), marker=m) for c, m in SLICE_MARKERS.items()}
#: Pr ∝ Ra^-1/2 slices are solid with filled markers; Pr ∝ Ra^-1 dashed with open markers.
FAMILY_STYLE = {-0.5: dict(ls="-", filled=True), -1.0: dict(ls="--", filled=False)}

#: Run status in parameter-space maps.
STATUS_STYLE = {
    "finished": dict(marker="o", mfc="k", mec="k"),
    "running":  dict(marker="o", mfc="none", mec="k"),
    "planned":  dict(marker="o", mfc="none", mec="0.6"),
}
#: Star (resolution study) size relative to the circles.
RES_CHECK_SCALE = 4 / 3


def _latex_ok() -> bool:
    """Whether Matplotlib's usetex will work: it needs latex and cm-super (type1ec.sty).

    Returns
    -------
    bool
    """
    import subprocess

    if not shutil.which("latex") or not shutil.which("kpsewhich"):
        return False
    r = subprocess.run(["kpsewhich", "type1ec.sty"], capture_output=True, text=True)
    return bool(r.stdout.strip())


def use(kind: str = "paper", usetex: bool | None = None) -> None:
    """Apply a named style.

    Parameters
    ----------
    kind : {"paper", "diag"}
        ``paper`` for publication figures; ``diag`` adds ``diag.mplstyle`` on top
        (bigger fonts and figures for per-run diagnostics).
    usetex : bool or None, optional
        LaTeX text rendering. None turns it on only where a working LaTeX exists
        (usually a laptop or login node, not a compute node).
    """
    import matplotlib.pyplot as plt

    sheets = [_HERE / "paper.mplstyle"]
    if kind != "paper":
        sheets.append(_HERE / f"{kind}.mplstyle")
    plt.style.use(sheets)
    if kind == "paper":
        if usetex is None:
            usetex = _latex_ok()
        mpl.rcParams["text.usetex"] = usetex
        if usetex:
            mpl.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"


def figsize(width: float = 1.0, aspect: float = 0.75) -> tuple[float, float]:
    """Figure size for a figure included at 100% in the paper.

    Parameters
    ----------
    width : float, optional
        Fraction of the text width (at most 1.5), or inches if larger than 1.5.
    aspect : float, optional
        Height / width.

    Returns
    -------
    tuple of float
        ``(width, height)`` in inches.
    """
    w = width * TEXTWIDTH_IN if width <= 1.5 else width
    return w, w * aspect


def save(fig, path, png: bool = False, close: bool = True) -> Path:
    """Save a figure as a vector PDF, and optionally a PNG preview.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The figure.
    path : str or pathlib.Path
        Output path; any suffix is replaced. Parent folders are created.
    png : bool, optional
        Also write a PNG.
    close : bool, optional
        Close the figure afterwards (ignored while :data:`KEEP_OPEN` is True).

    Returns
    -------
    pathlib.Path
        The PDF path.
    """
    import matplotlib.pyplot as plt

    path = Path(path).with_suffix("")
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path.with_suffix(".pdf"))
    if png:
        fig.savefig(path.with_suffix(".png"))
    if close and not KEEP_OPEN:
        plt.close(fig)
    return path.with_suffix(".pdf")


def slice_kw(c: float, a: float = -0.5, line: bool = False) -> dict:
    """Plot keyword arguments for one slice.

    Parameters
    ----------
    c : float
        Slice label log10 Pr_0.
    a : float, optional
        Slice slope (-0.5: filled markers; -1.0: open markers).
    line : bool, optional
        Return line keywords (colour, line style, outline if pale) instead of marker keywords.

    Returns
    -------
    dict
        Keywords for ``ax.plot`` or ``ax.errorbar``.
    """
    c = round(float(c) * 2) / 2
    kind = "slice-1" if a == -1.0 else "slice-0.5"
    col = expense_colour(c, kind)
    fam = FAMILY_STYLE.get(a, FAMILY_STYLE[-0.5])
    if line:
        return dict(line_kw(col), ls=fam["ls"])
    return dict(marker_kw(col, fam["filled"]), marker=SLICE_MARKERS.get(c, "x"), ls="none")


def pow10(e: float) -> str:
    """LaTeX for 10 to a power.

    Parameters
    ----------
    e : float
        Exponent (half-decades are written as ``10^{-2.5}``).

    Returns
    -------
    str
        LaTeX without dollar signs.
    """
    e = float(e)
    s = f"{e:g}"
    return rf"10^{{{s}}}"


def slice_label(c: float, a: float = -0.5) -> str:
    """Legend label of a slice.

    Parameters
    ----------
    c : float
        Slice label log10 Pr_0.
    a : float, optional
        Slice slope.

    Returns
    -------
    str
        Matplotlib mathtext/LaTeX label.
    """
    slope = {-0.5: r"Ra^{-1/2}", -1.0: r"Ra^{-1}"}.get(a, f"Ra^{{{a:g}}}")
    return rf"$Pr = {pow10(c)}\,(Ra/10^5)^{{{a:g}}}$" if a not in (-0.5, -1.0) else \
        rf"$Pr_0 = {pow10(c)}$, $Pr\propto {slope}$"
