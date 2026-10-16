"""Compensated Nu, Nu-1 or Re along slices Pr = Pr_0 (Ra/10^5)^a (figures A and B).

Reference law Q ~ Pr^alpha Ra^beta. Along a slice Pr ~ Ra^a, Q ~ Ra^(beta + alpha a), so
Q / Ra^(beta + alpha a) is constant if the law holds; the fitted slope m of the compensated
data is the departure from it. The defaults are the low-Pr ultimate sub-regime of
Shishkina & Lohse (2024), up to logarithms::

    Nu ~ Pr^1/2 Ra^1/2   ->  Ra^1/4 for a = -1/2,  Ra^0 for a = -1
    Re ~ Pr^-1/2 Ra^1/2  ->  Ra^3/4 for a = -1/2,  Ra^1 for a = -1

Vertical normalization (``--norm``). Dividing a slice by a constant shifts it on a log
axis and never changes its slope; the choice only decides how the slices sit relative to
each other:

``none``
    compensated values as they are (slice offsets carry the Pr dependence).
``first``
    each slice divided by its first point.
``ref``
    each slice divided by its value at Ra = 10^(``--ref-ra``): the data point if present,
    else the slice's own fit there.
``common``
    each slice divided by its geometric mean over the Ra values all slices share: the
    least-squares vertical alignment, no single point privileged.
``law``
    Q / (Pr^alpha Ra^beta) for every point, no per-slice constant: if the reference law
    holds, all slices collapse onto one horizontal line, so the residual offsets test the
    Pr exponent and the slopes test the slice exponent.

Slices are every line in results.csv with at least ``--min-points`` finished runs;
running points are drawn open and excluded from fits. Usage::

    python -m rbc.plots.scaling --q Nu-1 --norm common
    python -m rbc.plots.scaling --q Re --family -0.5 --fit linear quadratic --norm law
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from rbc import cases as rc
from rbc import style
from rbc.plots.common import decade_ticks, fit_loglog, proxy, slope_text, tidy_log_y

#: Reference exponents (alpha: Pr, beta: Ra) of the low-Pr ultimate sub-regime.
REF = {"Nu": (0.5, 0.5), "Nu-1": (0.5, 0.5), "Re": (-0.5, 0.5)}
FAMILY_TITLE = {-0.5: r"$Pr\propto Ra^{-1/2}$", -1.0: r"$Pr\propto Ra^{-1}$"}
NORMS = ("none", "first", "ref", "common", "law")


def exponent_tex(g: float) -> str:
    """An exponent as a simple fraction for axis labels (0.25 -> ``1/4``).

    Parameters
    ----------
    g : float
        Exponent.

    Returns
    -------
    str
    """
    from fractions import Fraction
    f = Fraction(g).limit_denominator(12)
    return str(f.numerator) if f.denominator == 1 else rf"{f.numerator}/{f.denominator}"


def qvals(df, q):
    """Values and standard errors of the plotted quantity.

    Parameters
    ----------
    df : pandas.DataFrame
        Results rows.
    q : {"Nu", "Nu-1", "Re"}
        Quantity.

    Returns
    -------
    tuple of pandas.Series
        Values, and errors (None unless every row has a finite error).
    """
    v, e = (df["Nu"] - 1, df.get("Nu_se")) if q == "Nu-1" else (df[q], df.get(f"{q}_se"))
    if e is not None and not np.all(np.isfinite(e)):
        e = None
    return v, e


def panel(ax, df, q, a, alpha, beta, fit_kinds, norm, ref_ra, min_points, include_running):
    """One family of compensated slices on one axes.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes.
    df : pandas.DataFrame
        Results rows.
    q : {"Nu", "Nu-1", "Re"}
        Quantity.
    a : float
        Slice slope (Pr ∝ Ra^a).
    alpha, beta : float
        Reference exponents of Pr and Ra; the data are divided by Ra^(beta + alpha a)
        (or by Pr^alpha Ra^beta with ``norm="law"``).
    fit_kinds : list of str
        ``"linear"`` and/or ``"quadratic"`` fits of each slice.
    norm : {"none", "first", "ref", "common", "law"}
        Vertical normalization (see the module docstring).
    ref_ra : float
        log10 Ra of the reference point for ``norm="ref"``.
    min_points : int
        Fewest finished points for a slice.
    include_running : bool
        Fit running points too.

    Returns
    -------
    list of dict
        One row per (slice, fit) for the fits CSV.
    """
    g = beta + alpha * a
    use = df if include_running else df[df["status"] == "finished"]
    slices = rc.slices(use, a, min_points)

    def comp(frame):
        """Compensated values (and errors) of some rows."""
        v, e = qvals(frame, q)
        div = frame["Pr"] ** alpha * frame["Ra"] ** beta if norm == "law" else frame["Ra"] ** g
        return v / div, (None if e is None else e / div)

    common = None
    if norm == "common":
        sets = [set(s["log_ra"]) for s in slices.values()]
        common = sorted(set.intersection(*sets)) if sets else []
        if not common:
            print(f"no Ra shared by all slices of family {a}; using 'first'")
            norm = "first"

    rows, handles, labels = [], [], []
    for c, s in slices.items():
        allpts = df[np.isclose(rc.slice_c(df["log_ra"], df["log_pr"], a), c)].sort_values("log_ra")
        y, sy = comp(s)
        y_all, _ = comp(allpts)
        fits = [fit_loglog(s["Ra"].values, y.values, None if sy is None else sy.values, k) for k in fit_kinds]
        if norm == "first":
            scale = y.iloc[0]
        elif norm == "ref":
            hit = np.isclose(s["log_ra"], ref_ra)
            scale = y[hit].iloc[0] if hit.any() else float(fits[0]["predict"](10.0 ** ref_ra))
        elif norm == "common":
            scale = 10 ** np.mean(np.log10(y[s["log_ra"].isin(common)]))
        else:
            scale = 1.0
        kw = style.slice_kw(c, a)
        done = allpts["status"] == "finished" if not include_running else np.ones(len(allpts), bool)
        ax.errorbar(s["Ra"], y / scale, yerr=None if sy is None else sy / scale, capsize=2, **kw)
        if (~done).any():
            ax.plot(allpts["Ra"][~done], (y_all / scale)[~done], **dict(kw, mfc="none"), alpha=0.6)
        xx = np.logspace(s["log_ra"].min(), s["log_ra"].max(), 100)
        for i, (kind, f) in enumerate(zip(fit_kinds, fits)):
            ax.plot(xx, f["predict"](xx) / scale, color=kw["color"], ls=["-", ":"][i], lw=1)
            rows.append(dict(q=q, a=a, c=c, kind=kind, slice_exponent=g, n=len(s), slope=f["coef"][-2],
                             slope_err=f["err"][-2], chi2_dof=f["chi2_dof"], aic=f["aic"], scale=scale,
                             **({"q2": f["coef"][0]} if kind == "quadratic" else {})))
        handles.append(proxy(kw))
        labels.append(f"${style.pow10(c)}$: {slope_text(fits[0])}")
        if norm == "none":
            ax.axhline(y.iloc[0], color=kw["color"], lw=0.5, ls=(0, (1, 2)))
    if norm in ("first", "ref", "common"):
        ax.axhline(1, color="0.4", lw=0.6, ls="--")
    if norm == "ref":
        ax.axvline(10 ** ref_ra, color="0.7", lw=0.6, ls=":")
    if norm == "common" and common:
        ax.axvspan(10 ** min(common), 10 ** max(common), color="0.92", zorder=0, lw=0)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(style.SYM["Ra"])
    qtex = "Nu-1" if q == "Nu-1" else q
    if norm == "law":
        ylab = rf"${qtex}\,/\,(Pr^{{{exponent_tex(alpha)}}}Ra^{{{exponent_tex(beta)}}})$"
    else:
        ylab = rf"${qtex}\,/\,Ra^{{{exponent_tex(g)}}}$" if g else rf"${qtex}$"
        ylab += {"first": " / first point", "ref": rf" / value at $Ra=10^{{{ref_ra:g}}}$",
                 "common": " / shared-Ra mean"}.get(norm, "")
    ax.set_ylabel(ylab)
    ax.set_title(FAMILY_TITLE.get(a, rf"$Pr\propto Ra^{{{a:g}}}$"))
    tidy_log_y(ax)
    decade_ticks(ax, "x")
    ax.legend(handles, labels, title=r"$Pr_0$: fitted slope $m$ of compensated data", loc="upper left",
              bbox_to_anchor=(0.0, -0.2), ncol=1, fontsize=7, title_fontsize=7, alignment="left")
    return rows


def get_parser() -> argparse.ArgumentParser:
    """Command-line options (see the module docstring).

    Returns
    -------
    argparse.ArgumentParser
    """
    p = argparse.ArgumentParser(prog="python -m rbc.plots.scaling", description=__doc__.split("\n")[0])
    p.add_argument("--q", choices=("Nu", "Nu-1", "Re"), default="Nu", help="quantity to compensate")
    p.add_argument("--family", type=float, nargs="+", default=[-0.5], help="slice slopes a to draw (one panel each)")
    p.add_argument("--alpha", type=float, help="reference exponent of Pr")
    p.add_argument("--beta", type=float, help="reference exponent of Ra")
    p.add_argument("--fit", nargs="+", choices=("linear", "quadratic"), default=["linear"], help="fits of log Q against log Ra along each slice")
    p.add_argument("--norm", choices=NORMS, default="none", help="vertical normalization of each slice (see the module docstring)")
    p.add_argument("--ref-ra", type=float, default=7.0, help="log10 Ra for --norm ref")
    p.add_argument("--min-points", type=int, default=3, help="fewest points for a slice to be drawn")
    p.add_argument("--include-running", action="store_true", help="also use points whose cases.csv status is not 'finished'")
    p.add_argument("--results", type=Path, default=rc.RESULTS_CSV, help="results table (default: data/results.csv)")
    p.add_argument("--out", type=Path, default=Path("figures"), help="output folder")
    p.add_argument("--width", type=float, default=1.0, help="figure width as a fraction of the text width")
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
    alpha, beta = REF[a.q]
    alpha = a.alpha if a.alpha is not None else alpha
    beta = a.beta if a.beta is not None else beta
    df = rc.load_results(a.results)
    df = df[np.isfinite(df["Nu" if a.q == "Nu-1" else a.q])]
    n = len(a.family)
    fig, axes = plt.subplots(1, n, figsize=style.figsize(a.width, 0.55 if n > 1 else 0.95), squeeze=False)
    rows = []
    for ax, fam in zip(axes[0], a.family):
        rows += panel(ax, df, a.q, fam, alpha, beta, a.fit, a.norm, a.ref_ra, a.min_points, a.include_running)
    tag = f"{a.q.replace('-', 'm')}_{'_'.join(a.fit)}_{a.norm}"
    print(style.save(fig, a.out / f"scaling_{tag}", png=a.png))
    fits = pd.DataFrame(rows)
    fits.to_csv(a.out / f"scaling_{tag}_fits.csv", index=False, float_format="%.6g")
    with pd.option_context("display.width", 140, "display.precision", 4):
        print(fits)


if __name__ == "__main__":
    main()
