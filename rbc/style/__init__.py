"""Single source of truth for how every figure looks.

Two Matplotlib style sheets: ``paper`` (publication figures, made at printed size)
and ``diag`` (large per-run diagnostic figures, layered on ``paper``). Colours and
markers of the slices, run statuses and resolution stars are defined here once, so
every figure agrees. Usage::

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

#: Slices Pr = Pr_0 (Ra/10^5)^{-1/2}, keyed by c = log10(Pr_0). Same colour and marker in
#: every figure (grid map, Nu/Re fits, contours, spectra comparisons). The palette was
#: checked for distinguishability under the common colour-vision deficiencies.
SLICE_STYLE = {
    2.0:  dict(color="#393b79", marker="*"),
    1.5:  dict(color="#8c564b", marker="p"),
    1.0:  dict(color="#17becf", marker="h"),
    0.5:  dict(color="#c49c2c", marker="d"),
    0.0:  dict(color="#1f77b4", marker="o"),
    -0.5: dict(color="#ff7f0e", marker="s"),
    -1.0: dict(color="#d62728", marker="D"),
    -1.5: dict(color="#9467bd", marker="^"),
    -2.0: dict(color="#2ca02c", marker="v"),
    -2.5: dict(color="#e377c2", marker="P"),
    -3.0: dict(color="#a5a61c", marker="X"),
}
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
    base = SLICE_STYLE.get(round(float(c) * 2) / 2, dict(color="0.3", marker="x"))
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
