"""How Nu - 1 decays with Pr: (Nu-1) / (Nu-1)|_{Pr=1} against Pr, one curve per Ra.

Dividing each Ra by its own value at Pr = 1 (``--ref-pr``) puts every curve through
(1, 1), so only the shape of the Pr dependence is compared, not the Ra dependence. The
lower panel is the local exponent of that decay,

    alpha_local = d log(Nu-1) / d log Pr,

from two-point secants between neighbouring Pr at each Ra (drawn at their geometric
midpoint). Horizontal guides mark the Pr exponents of the low-Pr regimes of
Grossmann & Lohse (2001, table I): I_l (Nu ~ Ra^1/4 Pr^1/8), II_l (Ra^1/5 Pr^1/5) and
IV_l (Ra^1/2 Pr^1/2, the same Pr exponent as the low-Pr ultimate sub-regime of Shishkina
& Lohse). A curve whose local exponent climbs past every guide as Pr falls is leaving
these turbulent regimes altogether: Nu - 1 is heading to 0, as convection weakens.

A Ra without a run at the reference Pr is left out (and listed); nothing is extrapolated.
Usage::

    python -m rbc.plots.pr_decay
    python -m rbc.plots.pr_decay --ref-pr -0.5 --guides 0.2 0.5
    python -m rbc.plots.pr_decay --results nusselt.csv          # a simple Ra,Pr,Nusselt file
    python -m rbc.plots.pr_decay --errorbars none

Error bars: ±1 standard error, and with ``--errorbars all`` a whisker adding the
worst-case systematic bound from ``Nu_sys`` (:mod:`rbc.plots.common`).
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from rbc import cases as rc
from rbc import style
from rbc.plots.common import add_errorbar_arg, decade_ticks, errorbar, plain_log_y, sys_values

LN10 = np.log(10.0)
#: Pr exponents of the low-Pr regimes of Grossmann & Lohse (2001, table I).
GL_GUIDES = {0.125: r"I$_l$", 0.2: r"II$_l$", 0.5: r"IV$_l$"}


def qvals(df, q):
    """Values and standard errors (None where unknown) of the plotted quantity.

    Parameters
    ----------
    df : pandas.DataFrame
        Results rows.
    q : {"Nu-1", "Nu"}
        Quantity.

    Returns
    -------
    tuple of numpy.ndarray
        Values, errors (None unless every row has a finite error) and systematic bounds
        (NaN where unknown).
    """
    v = (df["Nu"] - 1 if q == "Nu-1" else df["Nu"]).to_numpy(float)
    e = df["Nu_se"].to_numpy(float) if "Nu_se" in df else None
    if e is not None and not np.all(np.isfinite(e)):
        e = None
    return v, e, sys_values(df, q)


def curves(df, q, ref_pr, min_points=2):
    """Normalized curve and local exponents for every Ra that has the reference Pr.

    Parameters
    ----------
    df : pandas.DataFrame
        Results rows (``log_ra``, ``log_pr``, ``Pr``, ``Nu``, optionally ``Nu_se``).
    q : {"Nu-1", "Nu"}
        Quantity.
    ref_pr : float
        log10 of the reference Pr.
    min_points : int, optional
        Fewest Pr values for a curve.

    Returns
    -------
    rows : list of dict
        Per Ra: ``log_ra``, ``Pr``, ratio ``r``, its error ``sr`` and systematic bound
        ``br``, and the secants ``Pr_mid``, ``alpha``, ``alpha_err``, ``alpha_sys``.
        Systematic bounds add linearly (worst case, :mod:`rbc.analysis.exponents`).
    skipped : list of float
        log10 Ra values without a run at the reference Pr.
    """
    rows, skipped = [], []
    for lra, g in df.groupby("log_ra"):
        g = g.sort_values("log_pr")
        if g["log_pr"].nunique() < min_points:
            continue
        ref = np.isclose(g["log_pr"], ref_pr)
        if not ref.any():
            skipped.append(float(lra))
            continue
        v, e, b = qvals(g, q)
        v0 = v[ref][0]
        brel = b / v
        rel = None if e is None else e / v                     # relative errors
        r = v / v0
        # independent runs: relative errors add in quadrature; the reference point is exactly 1
        sr = None if rel is None else np.where(ref, 0.0, r * np.sqrt(rel ** 2 + rel[ref][0] ** 2))
        L = np.log10(v)
        x = np.log10(g["Pr"].to_numpy(float))
        dx = np.diff(x)
        alpha = np.diff(L) / dx
        aerr = None if rel is None else np.sqrt((rel[1:] / LN10) ** 2 + (rel[:-1] / LN10) ** 2) / dx
        br = np.where(ref, 0.0, r * (brel + brel[ref][0]))
        rows.append(dict(log_ra=float(lra), Pr=g["Pr"].to_numpy(float), r=r, sr=sr, br=br,
                         Pr_mid=10 ** (0.5 * (x[1:] + x[:-1])), alpha=alpha, alpha_err=aerr,
                         alpha_sys=(brel[1:] + brel[:-1]) / (LN10 * dx)))
    return rows, skipped


def plot(ax, axl, rows, q, ref_pr, guides, bars="all"):
    """Draw the normalized curves (``ax``) and their local exponents (``axl``).

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes for (Q / Q at the reference Pr) against Pr.
    axl : matplotlib.axes.Axes or None
        Axes for the local exponent (None: not drawn).
    rows : list of dict
        From :func:`curves`.
    q : str
        Quantity.
    ref_pr : float
        log10 of the reference Pr.
    guides : list of float
        Exponents drawn as horizontal guides on ``axl``.
    bars : {"all", "stat", "none"}, optional
        Error bars drawn.
    """
    for r in rows:
        kw = style.category_kw(r["log_ra"], "ra")
        lab = rf"$Ra={style.pow10(r['log_ra'])}$"
        errorbar(ax, r["Pr"], r["r"], r["sr"], r["br"], bars, ms=3, lw=0.9, capsize=1.5, label=lab, **kw)
        if axl is not None:
            errorbar(axl, r["Pr_mid"], r["alpha"], r["alpha_err"], r["alpha_sys"], bars, ms=3, lw=0.9,
                     capsize=1.5, **kw)
    qt = "Nu-1" if q == "Nu-1" else "Nu"
    ref_tex = style.pow10(ref_pr)
    ax.axhline(1, color="0.6", lw=0.6, ls=":")
    ax.axvline(10 ** ref_pr, color="0.6", lw=0.6, ls=":")
    ax.set_xscale("log")
    ax.set_yscale("log")
    decade_ticks(ax, "x")
    plain_log_y(ax)
    ax.set_ylabel(rf"$({qt})\,/\,({qt})_{{Pr={ref_tex}}}$")
    ax.legend(fontsize=7, loc="lower right")
    if axl is None:
        ax.set_xlabel(style.SYM["Pr"])
        return
    for gexp in guides:
        name = GL_GUIDES.get(round(gexp, 3))
        axl.axhline(gexp, color="0.5", lw=0.7, ls="--", zorder=0)
        axl.annotate(f"{gexp:g}" + (f" ({name})" if name else ""), (1.01, gexp), xycoords=("axes fraction", "data"),
                     fontsize=6.5, va="center", color="0.35", annotation_clip=False)
    axl.set_xscale("log")
    decade_ticks(axl, "x")
    axl.set_xlabel(style.SYM["Pr"])
    axl.set_ylabel(r"local $\alpha = d\log(" + qt + r")/d\log Pr$", fontsize=8)


def get_parser() -> argparse.ArgumentParser:
    """Command-line options (see the module docstring).

    Returns
    -------
    argparse.ArgumentParser
    """
    p = argparse.ArgumentParser(prog="python -m rbc.plots.pr_decay", description=__doc__.split("\n")[0])
    p.add_argument("--q", choices=("Nu-1", "Nu"), default="Nu-1", help="quantity")
    p.add_argument("--ref-pr", type=float, default=0.0, help="log10 of the Pr each curve is divided by its value at")
    p.add_argument("--guides", type=float, nargs="*", default=sorted(GL_GUIDES),
                   help="exponents drawn as guides on the local-exponent panel (default: GL I_l, II_l, IV_l)")
    p.add_argument("--no-local", action="store_true", help="leave out the local-exponent panel")
    p.add_argument("--min-points", type=int, default=2, help="fewest Pr values for a curve")
    p.add_argument("--include-running", action="store_true", help="also use points whose status is not 'finished'")
    p.add_argument("--results", type=Path, default=rc.RESULTS_CSV, help="results table (default: data/results.csv)")
    p.add_argument("--out", type=Path, default=Path("figures"), help="output folder")
    p.add_argument("--png", action="store_true", help="also write PNG previews")
    p.add_argument("--no-tex", action="store_true", help="never use LaTeX for text (default: use it where it is installed)")
    add_errorbar_arg(p)
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
    rc.require_columns(df, ["Nu"], "the Pr-decay figure")
    df = df[np.isfinite(df["Nu"])]
    if not a.include_running:
        df = df[df["status"] == "finished"]
    rows, skipped = curves(df, a.q, a.ref_pr, a.min_points)
    if skipped:
        print("no run at the reference Pr for Ra = " + ", ".join(f"1e{s:g}" for s in skipped) + "; left out")
    if not rows:
        raise SystemExit(f"no Ra has a run at Pr = 10^{a.ref_pr:g}")
    w, _ = style.figsize(0.8)
    if a.no_local:
        fig, ax = plt.subplots(figsize=(w, 0.8 * w), layout="constrained")
        axl = None
    else:
        fig, (ax, axl) = plt.subplots(2, 1, figsize=(w, 1.05 * w), sharex=True, layout="constrained",
                                      height_ratios=[1.3, 1])
    plot(ax, axl, rows, a.q, a.ref_pr, a.guides, a.errorbars)
    tag = f"pr_decay_{a.q.replace('-', 'm')}" + ("" if a.ref_pr == 0 else f"_ref{a.ref_pr:g}")
    print(style.save(fig, a.out / tag, png=a.png))
    table = pd.DataFrame([dict(log_ra=r["log_ra"], Pr_mid=pm, alpha=al,
                               alpha_err=np.nan if r["alpha_err"] is None else ae, alpha_sys=asy)
                          for r in rows for pm, al, ae, asy in zip(r["Pr_mid"], r["alpha"],
                                                                   r["alpha_err"] if r["alpha_err"] is not None
                                                                   else r["alpha"], r["alpha_sys"])])
    a.out.mkdir(parents=True, exist_ok=True)
    table.to_csv(a.out / f"{tag}_local.csv", index=False, float_format="%.6g")


if __name__ == "__main__":
    main()
