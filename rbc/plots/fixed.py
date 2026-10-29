"""Nu-1 (or Nu, Re) against one parameter with the other held fixed (figure D).

Usage::

    python -m rbc.plots.fixed                          # Nu-1 vs Pr at fixed Ra | vs Ra at fixed Pr
    python -m rbc.plots.fixed --slopes                 # fitted exponent of each curve in the legend
    python -m rbc.plots.fixed --norm comp --fit linear # compensated, with a fit per curve
    python -m rbc.plots.fixed --norm law --ref sl      # against the Shishkina-Lohse law
    python -m rbc.plots.fixed --q Nu --along ra --lit data/literature.csv

Nu-1 is the convective part of the heat flux, the quantity the exact relations
<eps_u> = (Nu-1)/sqrt(Ra Pr) and the scaling theories refer to; at Nu ~ 3 the local
exponent of Nu is (Nu-1)/Nu = 2/3 of that of Nu-1, so the choice matters at low Ra.

Normalization (``--norm``), with reference exponents Q ~ Pr^alpha Ra^beta (``--ref``):

``none``
    Q itself.
``comp``
    compensated by the plotted parameter: Q / Pr^alpha against Pr, Q / Ra^beta against Ra.
    A curve is flat where its local exponent equals the reference one, and the fitted
    slope m is the departure from it.
``collapse``
    divided by the fixed parameter: Q / Ra^beta against Pr, Q / Pr^alpha against Ra. If
    Q separates as f(Pr) g(Ra) with that exponent, the curves fall on top of each other;
    the slopes are the curves' own exponents, unchanged. The dashed black line is the
    reference law with its best prefactor.
``law``
    Q / (Pr^alpha Ra^beta) in both panels: a single power law would make every curve the
    dashed horizontal line, so what remains is everything the law misses (``collapse``
    with the trend taken out, so the departures are magnified).
``common``
    each curve divided by its geometric mean over the values of the plotted parameter that
    every curve shares: no reference law at all, only the shapes are compared.

In the legend, alpha or beta is a curve's fitted exponent, and m (``comp`` and ``law``)
its departure from the reference exponent.

``--ref fit`` (default) takes alpha and beta from one power law fitted to all points;
``--ref sl`` uses the Shishkina-Lohse values 1/2, 1/2; ``--alpha`` and ``--beta`` override
either. ``--fit linear`` (or ``quadratic``) draws a log-log fit through each curve, as in
:mod:`rbc.plots.scaling`, and puts its slope in the legend.

--lit: optional CSV with columns source,Ra,Pr,Nu[,Re] (other studies) drawn in grey,
one marker per source, for checking agreement with published data.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from rbc import cases as rc
from rbc import style
from rbc.plots.common import decade_ticks, fit_loglog, plain_log_y

#: Markers cycled through the literature sources.
LIT_MARKERS = "osD^v<>PX*"
NORMS = ("none", "comp", "collapse", "law", "common")
#: Shishkina-Lohse reference exponents (alpha: Pr, beta: Ra) per quantity.
SL_REF = {"Nu": (0.5, 0.5), "Nu-1": (0.5, 0.5), "Re": (-0.5, 0.5)}


def qcol(df, q):
    """Values and errors of the plotted quantity.

    Parameters
    ----------
    df : pandas.DataFrame
        Results rows.
    q : {"Nu-1", "Nu", "Re"}
        Quantity.

    Returns
    -------
    tuple of pandas.Series
        Values and errors (errors None if the column is absent).
    """
    if q == "Nu-1":
        return df["Nu"] - 1, df.get("Nu_se")
    return df[q], df.get(f"{q}_se")


def _exp_tex(v: float) -> str:
    """An exponent for a label: a simple fraction if it is one, else 3 decimals."""
    f = Fraction(v).limit_denominator(12)
    if abs(float(f) - v) < 1e-9:
        return str(f.numerator) if f.denominator == 1 else f"{f.numerator}/{f.denominator}"
    return f"{v:.3f}"


def divisor(frame, along, norm, alpha, beta):
    """What each row's value is divided by for a normalization.

    Parameters
    ----------
    frame : pandas.DataFrame
        Rows with ``Pr`` and ``Ra``.
    along : {"pr", "ra"}
        Plotted parameter.
    norm : str
        One of :data:`NORMS`; ``common`` is per curve and handled in :func:`panel`.
    alpha, beta : float
        Reference exponents of Pr and Ra.

    Returns
    -------
    numpy.ndarray
    """
    pr, ra = frame["Pr"].to_numpy(float), frame["Ra"].to_numpy(float)
    if norm == "comp":
        return pr ** alpha if along == "pr" else ra ** beta
    if norm == "collapse":
        return ra ** beta if along == "pr" else pr ** alpha
    if norm == "law":
        return pr ** alpha * ra ** beta
    return np.ones(len(frame))


def ylabel(q, along, norm, alpha, beta):
    """Axis label for a quantity and normalization.

    Parameters
    ----------
    q : str
        Quantity.
    along : {"pr", "ra"}
        Plotted parameter.
    norm : str
        Normalization.
    alpha, beta : float
        Reference exponents.

    Returns
    -------
    str
    """
    qt = "Nu-1" if q == "Nu-1" else q
    pa, rb = rf"Pr^{{{_exp_tex(alpha)}}}", rf"Ra^{{{_exp_tex(beta)}}}"
    den = {"comp": pa if along == "pr" else rb, "collapse": rb if along == "pr" else pa,
           "law": pa + r"\," + rb}.get(norm)
    if norm == "common":
        return rf"$({qt})\,/$ curve mean"
    return rf"${qt}$" if den is None else rf"$({qt})\,/\,{den}$"


def prefactor(df, q, alpha, beta):
    """Best prefactor C of Q = C Pr^alpha Ra^beta with the exponents held fixed.

    Parameters
    ----------
    df : pandas.DataFrame
        Results rows.
    q : str
        Quantity.
    alpha, beta : float
        Exponents of Pr and Ra.

    Returns
    -------
    float
        Geometric mean of Q / (Pr^alpha Ra^beta) over all rows.
    """
    y, _ = qcol(df, q)
    return float(10 ** np.mean(np.log10(y.to_numpy() / divisor(df, "pr", "law", alpha, beta))))


def panel(ax, df, q, along, lit=None, slopes=False, leg_ax=None, norm="none", alpha=0.5, beta=0.5,
          fit=None, c=None):
    """One curve per value of the fixed parameter.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes for the curves.
    df : pandas.DataFrame
        Results rows.
    q : {"Nu-1", "Nu", "Re"}
        Quantity.
    along : {"pr", "ra"}
        Against Pr at fixed Ra, or against Ra at fixed Pr.
    lit : pandas.DataFrame, optional
        Published data to overlay.
    slopes : bool, optional
        Add each curve's fitted slope (of the plotted, normalized values) to its legend entry.
    leg_ax : matplotlib.axes.Axes, optional
        Empty axes beside ``ax`` that holds the legend as one column (default: inside
        ``ax``).
    norm : str, optional
        One of :data:`NORMS` (see the module docstring).
    alpha, beta : float, optional
        Reference exponents of Pr and Ra for the normalization.
    fit : {"linear", "quadratic"}, optional
        Draw a log-log fit through each curve (points are then not joined). A quadratic
        needs four points; shorter curves get a linear fit. The slope reported is the one
        at the middle of the curve (geometric mean of its x values).
    c : float, optional
        Prefactor of the reference law (:func:`prefactor`), drawn as a dashed black line
        for ``collapse`` and ``law``.

    Returns
    -------
    list of dict
        One row per curve: fixed value, fitted slope and its error.
    """
    other, x = ("log_ra", "Pr") if along == "pr" else ("log_pr", "Ra")
    vals = sorted(df[other].unique())
    pkind = "ra" if other == "log_ra" else "pr"     # colours by expense of the fixed parameter
    sym = r"\alpha" if along == "pr" else r"\beta"
    ref = alpha if along == "pr" else beta
    departure = norm in ("comp", "law")      # slopes are departures from the reference
    if departure:
        sym = "m"
    shared = None
    if norm == "common":
        sets = [set(df.loc[df[other] == v, x]) for v in vals if (df[other] == v).sum() >= 2]
        shared = set.intersection(*sets) if sets else set()
        if not shared:
            print(f"no {x} value shared by every curve; normalizing each by its own mean instead")
    rows = []
    for v in vals:
        g = df[df[other] == v].sort_values(x)
        if len(g) < 2:
            continue
        ckw = style.category_kw(v, pkind)
        col = ckw["color"]
        y, err = qcol(g, q)
        div = divisor(g, along, norm, alpha, beta)
        if norm == "common":
            sel = g[x].isin(shared).to_numpy() if shared else np.ones(len(g), bool)
            div = np.full(len(g), 10 ** np.mean(np.log10(y.to_numpy()[sel])))
        y = y / div
        err = None if err is None else err / div
        lab = rf"${style.pow10(v)}$"
        f = None
        if slopes or fit:
            kind = fit if fit == "quadratic" and len(g) >= 4 else "linear"
            # fitted about the curve's middle, so a quadratic's slope is the one there
            xc = 10 ** np.mean(np.log10(g[x].to_numpy(float)))
            f = fit_loglog(g[x].to_numpy() / xc, y.to_numpy(), None if err is None else err.to_numpy(), kind)
            rows.append(dict(fixed=other, value=v, slope=f["coef"][-2], slope_err=f["err"][-2], fit=kind,
                             n=len(g)))
        if slopes:      # value first, so the numbers line up in the legend column
            lab = rf"${sym} = {f['coef'][-2]:+.3f}$ at {lab}" if departure else \
                rf"${sym} = {f['coef'][-2]:.3f}$ at {lab}"
        ax.errorbar(g[x], y, yerr=err, color=col, marker=ckw["marker"], ms=3, lw=0 if fit else 0.9,
                    capsize=1.5, label=lab)
        if fit:
            xx = np.logspace(np.log10(g[x].min()), np.log10(g[x].max()), 60)
            ax.plot(xx, f["predict"](xx / xc), color=col, lw=0.9)
    if lit is not None:
        for i, (src, g) in enumerate(lit.groupby("source")):
            if q == "Nu-1" and "Nu" not in g or q != "Nu-1" and q not in g:
                continue
            y, _ = qcol(g, q)
            if norm not in ("none", "common"):
                y = y / divisor(g, along, norm, alpha, beta)
            ax.plot(g[x], y, ls="none", marker=LIT_MARKERS[i % len(LIT_MARKERS)], mfc="none",
                    mec="0.45", ms=4, label=src)
    if norm == "common":
        ax.axhline(1, color="0.5", lw=0.6, ls="--")
    elif c is not None and norm in ("collapse", "law"):
        xx = np.logspace(np.log10(df[x].min()), np.log10(df[x].max()), 60)
        yy = np.full_like(xx, c) if norm == "law" else c * xx ** ref
        var = "Pr" if along == "pr" else "Ra"
        ax.plot(xx, yy, color="k", lw=0.8, ls="--", zorder=1,
                label=rf"$C\,{var}^{{{_exp_tex(ref)}}}$" if norm == "collapse" else "reference law")
    ax.set_xscale("log")
    ax.set_yscale("log")
    decade_ticks(ax, "x")
    plain_log_y(ax)
    ax.set_xlabel(style.SYM["Pr" if along == "pr" else "Ra"])
    ax.set_ylabel(ylabel(q, along, norm, alpha, beta))
    fixed = "Ra" if along == "pr" else "Pr"
    title = rf"fixed ${fixed}$"
    if slopes:
        own = r"\alpha" if along == "pr" else r"\beta"
        title += rf": $m = {own} - {ref:.3f}$" if departure else " (fitted exponent)"
        if fit == "quadratic":
            title += "\nslopes at mid-curve"
    h, l = ax.get_legend_handles_labels()
    if leg_ax is None:
        ax.legend(h, l, title=title, fontsize=6.5, title_fontsize=7, loc="best")
    else:
        leg_ax.axis("off")
        leg_ax.legend(h, l, title=title, loc="center left", borderaxespad=0, fontsize=7,
                      title_fontsize=7, handlelength=1.4, handletextpad=0.5, labelspacing=0.35,
                      alignment="left")
    return rows


def reference_exponents(df, q, ref):
    """Reference exponents (alpha, beta) for the normalization.

    Parameters
    ----------
    df : pandas.DataFrame
        Results rows.
    q : str
        Quantity.
    ref : {"fit", "sl"}
        One power law fitted to all points, or the Shishkina-Lohse values.

    Returns
    -------
    tuple of float
    """
    if ref == "sl":
        return SL_REF[q]
    from rbc.analysis import exponents as ex
    s = ex.Surface(ex.prepare(df, q), 1)
    return float(s.exponent(0, 0, "alpha")[0]), float(s.exponent(0, 0, "beta")[0])


def get_parser() -> argparse.ArgumentParser:
    """Command-line options (see the module docstring).

    Returns
    -------
    argparse.ArgumentParser
    """
    p = argparse.ArgumentParser(prog="python -m rbc.plots.fixed", description=__doc__.split("\n")[0])
    p.add_argument("--q", choices=("Nu-1", "Nu", "Re"), default="Nu-1", help="quantity on the y axis")
    p.add_argument("--along", choices=("pr", "ra"), nargs="+", default=["pr", "ra"], help="pr: against Pr at fixed Ra; ra: against Ra at fixed Pr")
    p.add_argument("--norm", choices=NORMS, default="none", help="normalization (see the module docstring)")
    p.add_argument("--ref", choices=("fit", "sl"), default="fit",
                   help="reference exponents: one power law fitted to all points, or Shishkina-Lohse")
    p.add_argument("--alpha", type=float, help="reference exponent of Pr (overrides --ref)")
    p.add_argument("--beta", type=float, help="reference exponent of Ra (overrides --ref)")
    p.add_argument("--fit", choices=("linear", "quadratic"), help="draw a log-log fit through each curve")
    p.add_argument("--slopes", action="store_true", help="put each curve's fitted slope in the legend")
    p.add_argument("--lit", type=Path, help="CSV of published data (source,Ra,Pr,Nu[,Re]) drawn in grey")
    p.add_argument("--results", type=Path, default=rc.RESULTS_CSV, help="results table (default: data/results.csv)")
    p.add_argument("--out", type=Path, default=Path("figures"), help="output folder")
    p.add_argument("--png", action="store_true", help="also write PNG previews")
    p.add_argument("--no-tex", action="store_true", help="never use LaTeX for text (default: use it where it is installed)")
    return p


def main(argv=None) -> None:
    """Command-line entry point.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments (default: ``sys.argv[1:]``).
    """
    a = get_parser().parse_args(argv)

    style.use("paper", usetex=False if a.no_tex else None)
    df = rc.load_results(a.results)
    rc.require_columns(df, ["Nu" if a.q == "Nu-1" else a.q], "the fixed-parameter figure")
    df = df[(df["status"] == "finished") & np.isfinite(df["Nu" if a.q == "Nu-1" else a.q])]
    alpha, beta = reference_exponents(df, a.q, a.ref) if a.norm not in ("none", "common") else (0.0, 0.0)
    alpha = a.alpha if a.alpha is not None else alpha
    beta = a.beta if a.beta is not None else beta
    if a.norm not in ("none", "common"):
        print(f"normalizing with alpha = {alpha:.4f} (Pr), beta = {beta:.4f} (Ra)")
    lit = pd.read_csv(a.lit) if a.lit else None
    n = len(a.along)
    w, _ = style.figsize(1.0)
    fig = plt.figure(figsize=(w, 0.37 * w * n + 0.1))
    gs = fig.add_gridspec(n, 2, width_ratios=[2.5, 1])
    c = prefactor(df, a.q, alpha, beta) if a.norm in ("collapse", "law") else None
    rows = []
    for k, along in enumerate(a.along):
        rows += panel(fig.add_subplot(gs[k, 0]), df, a.q, along, lit, a.slopes, fig.add_subplot(gs[k, 1]),
                      a.norm, alpha, beta, a.fit, c)
    tag = f"fixed_{a.q.replace('-', 'm')}_{'_'.join(a.along)}" + ("" if a.norm == "none" else f"_{a.norm}") \
        + (f"_{a.ref}" if a.norm in ("comp", "collapse", "law") else "") + (f"_{a.fit}" if a.fit else "")
    print(style.save(fig, a.out / tag, png=a.png))
    if rows:
        a.out.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(rows).to_csv(a.out / f"{tag}_slopes.csv", index=False, float_format="%.6g")


if __name__ == "__main__":
    main()
