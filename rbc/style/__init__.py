"""Single source of truth for how every figure looks.

Two Matplotlib style sheets: ``paper`` (publication figures, made at printed size)
and ``diag`` (large per-run diagnostic figures, layered on ``paper``). Colours and
markers of the slices, of curves at one Ra or one Pr, run statuses and resolution stars
are defined here once, so every figure agrees: colours run from red (cheapest runs) to
violet (most expensive), see :data:`EXPENSE_ANCHORS`. Usage::

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

#: Colours by expense. Every figure that colours runs, curves or slices by a parameter
#: orders the colours like the visible spectrum, red (cheapest) through orange, yellow,
#: green, blue and indigo to violet (most expensive), as photon energy rises. Higher Ra,
#: lower Pr and lower slice Pr_0 are more expensive. Yellow is darkened so it reads on
#: white. Continuous fields (contours) keep an ordinary colormap.
EXPENSE_ANCHORS = ["#d62728", "#f0811a", "#d9a90b", "#7aa821", "#239a49", "#14a0a6", "#2c6fcc",
                   "#4337b0", "#8a3ab9"]
#: The value each kind of parameter takes at the cheap (red) and expensive (violet) ends,
#: fixed so that a value has the same colour in every figure; values beyond are clipped.
#: ``ra``: log10 Ra; ``pr``: log10 Pr; ``slice-0.5`` and ``slice-1``: c = log10 Pr_0 of the
#: slices Pr ∝ Ra^-1/2 and Pr ∝ Ra^-1 (each family uses the whole spectrum; the two
#: families share their labels, not their runs).
EXPENSE_RANGE = {"ra": (5.0, 9.0), "pr": (0.0, -3.5), "slice-0.5": (1.0, -3.0), "slice-1": (2.0, -1.5)}
#: Markers by position on the half-decade grid, so colour is never the only cue.
GRID_MARKERS = "osD^vph*PX<>d"


def expense_cmap():
    """The red-to-violet colormap behind every expense colour.

    Returns
    -------
    matplotlib.colors.LinearSegmentedColormap
    """
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list("expense", EXPENSE_ANCHORS)


def expense_position(value: float, kind: str) -> float:
    """Position of a parameter value on the expense scale, 0 (cheapest) to 1.

    Parameters
    ----------
    value : float
        log10 Ra, log10 Pr or slice label c, according to ``kind``.
    kind : {"ra", "pr", "slice-0.5", "slice-1"}
        Which parameter.

    Returns
    -------
    float
    """
    cheap, dear = EXPENSE_RANGE[kind]
    return float(min(max((value - cheap) / (dear - cheap), 0.0), 1.0))


def expense_colour(value: float, kind: str):
    """Colour of a parameter value (see :data:`EXPENSE_ANCHORS`).

    Parameters
    ----------
    value : float
        log10 Ra, log10 Pr or slice label c.
    kind : {"ra", "pr", "slice-0.5", "slice-1"}
        Which parameter.

    Returns
    -------
    tuple
        RGBA.
    """
    return expense_cmap()(expense_position(value, kind))


def expense_norm(kind: str):
    """Colormap and normalization for a colour bar of a parameter, cheap end red.

    Parameters
    ----------
    kind : {"ra", "pr", "slice-0.5", "slice-1"}
        Which parameter.

    Returns
    -------
    tuple
        ``(cmap, norm)``; the norm runs over the parameter itself, increasing.
    """
    from matplotlib.colors import Normalize
    cheap, dear = EXPENSE_RANGE[kind]
    cmap = expense_cmap()
    return (cmap, Normalize(cheap, dear)) if dear > cheap else (cmap.reversed(), Normalize(dear, cheap))


def ranked_colours(n: int) -> list:
    """``n`` colours from red to violet, for things ranked by expense within one figure
    (resolutions of one point, say).

    Parameters
    ----------
    n : int
        Number of colours.

    Returns
    -------
    list of tuple
    """
    cmap = expense_cmap()
    return [cmap(i / max(n - 1, 1)) for i in range(n)]


def category_kw(value: float, kind: str, line: bool = False) -> dict:
    """Plot keywords for one curve labelled by a parameter value (one Ra, or one Pr).

    Parameters
    ----------
    value : float
        log10 Ra or log10 Pr.
    kind : {"ra", "pr"}
        Which parameter.
    line : bool, optional
        Colour only (for lines), instead of colour and marker.

    Returns
    -------
    dict
    """
    col = expense_colour(value, kind)
    if line:
        return dict(color=col)
    cheap, dear = EXPENSE_RANGE[kind]
    k = int(round(abs(value - cheap) / 0.5))
    return dict(color=col, marker=GRID_MARKERS[k % len(GRID_MARKERS)])


#: Markers of the slices Pr = Pr_0 (Ra/10^5)^a, keyed by c = log10(Pr_0); colours are
#: :func:`expense_colour` of c, so a lower Pr_0 (more expensive) is further towards violet.
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
        Return line keywords (colour, line style) instead of marker keywords.

    Returns
    -------
    dict
        Keywords for ``ax.plot`` or ``ax.errorbar``.
    """
    c = round(float(c) * 2) / 2
    kind = "slice-1" if a == -1.0 else "slice-0.5"
    base = dict(color=expense_colour(c, kind), marker=SLICE_MARKERS.get(c, "x"))
    fam = FAMILY_STYLE.get(a, FAMILY_STYLE[-0.5])
    if line:
        return dict(color=base["color"], ls=fam["ls"])
    return dict(color=base["color"], marker=base["marker"], ls="none",
                mfc=base["color"] if fam["filled"] else "none", mec=base["color"])


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
