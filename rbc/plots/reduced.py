"""Nu - 1 along slices against one combined variable x = Ra Pr^p, and its local exponent.

For slices Pr ∝ Ra^a (``--family``), every point is placed at x = Ra Pr^p, by default
with p = a: x = Ra Pr^-1/2 for the slices Pr ∝ Ra^-1/2, and x = Ra Pr^-1 = Ra/Pr = Re_ff^2
(the free-fall Reynolds number squared) for Pr ∝ Ra^-1. ``--pr-power`` sets p for every
family. Top row: Nu - 1 against x, one curve per slice. Bottom row: the local exponent of
Nu - 1 ~ c x^gamma,

    gamma = d log(Nu-1) / d log x,

from two-point secants between neighbouring points of each slice (drawn at the geometric
mean of their x), and for slices of four points or more the derivative of a quadratic fit
of log(Nu-1) against log x.

What gamma is. Along a slice log x = (1 + p a) log Ra + const, so gamma is the slice
exponent m = d log(Nu-1) / d log Ra (the slope in :mod:`rbc.plots.scaling`) divided by
(1 + p a): 5/4 for p = a = -1/2, and 2 for p = a = -1. In terms of the local exponents of
Nu - 1 ~ Pr^alpha Ra^beta, gamma = (beta + a alpha) / (1 + p a). What the combined variable
adds is a common horizontal axis for every slice: if the gamma(x) of all slices fall on
one curve, x is what organizes the approach to the asymptotic exponent; if they stay
ordered by slice, it is not. Nu - 1 itself does not collapse unless alpha = p beta.

Reference lines: the exponent the low-Pr ultimate law of Shishkina & Lohse (Nu ~ Pr^1/2
Ra^1/2, also Grossmann & Lohse's IV_l) gives, (1 + a) / (2 (1 + p a)), which is 0 for
a = -1 because that law depends on Ra Pr only, which is fixed along those slices; and
the one from a single power law fitted to all points (``--no-fit-line`` leaves it out).

Usage::

    python -m rbc.plots.reduced                                  # both families
    python -m rbc.plots.reduced --family -1 --min-points 3
    python -m rbc.plots.reduced --pr-power -1                    # x = Ra/Pr for every family
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from rbc import cases as rc
from rbc import style
from rbc.plots.common import fit_loglog, plain_log_y, proxy

LN10 = np.log(10.0)
FAMILY_TITLE = {-0.5: r"$Pr\propto Ra^{-1/2}$", -1.0: r"$Pr\propto Ra^{-1}$"}


def _frac(v: float) -> str:
    """An exponent for a label: a simple fraction if it is one."""
    from fractions import Fraction
    f = Fraction(v).limit_denominator(12)
    return str(f.numerator) if f.denominator == 1 else f"{f.numerator}/{f.denominator}"


def sparse_decades(ax, lo, hi):
    """Labelled ticks at every decade of the x axis, or every second or third when wide.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes (log x).
    lo, hi : float
        Data range of x.
    """
    from matplotlib import ticker
    k0, k1 = int(np.floor(np.log10(lo))), int(np.ceil(np.log10(hi)))
    step = 1 if k1 - k0 <= 4 else 2 if k1 - k0 <= 8 else 3
    ax.xaxis.set_major_locator(ticker.FixedLocator([10.0 ** k for k in range(k0, k1 + 1) if k % step == 0]))
    ax.xaxis.set_major_formatter(ticker.LogFormatterMathtext())
    ax.xaxis.set_minor_locator(ticker.LogLocator(base=10, subs=np.arange(2, 10) * 0.1, numticks=60))
    ax.xaxis.set_minor_formatter(ticker.NullFormatter())


def x_label(p: float) -> str:
    """Axis label of x = Ra Pr^p.

    Parameters
    ----------
    p : float
        Power of Pr.

    Returns
    -------
    str
    """
    if p == -1:
        return r"$Ra\,Pr^{-1} = Re_{\rm ff}^2$"
    return rf"$Ra\,Pr^{{{_frac(p)}}}$"


def slice_exponents(s, p, min_q=0.0):
    """Secant and quadratic-fit local exponents of Nu - 1 against x along one slice.

    Parameters
    ----------
    s : pandas.DataFrame
        The slice's rows, sorted by Ra.
    p : float
        Power of Pr in x = Ra Pr^p.
    min_q : float, optional
        Leave out points with Nu - 1 below this (weak convection).

    Returns
    -------
    dict
        ``x``, ``y`` (Nu - 1), ``sy`` (errors or None); secants ``x_mid``, ``gamma``,
        ``gamma_err`` (or None); and ``fit`` (:func:`rbc.plots.common.fit_loglog`
        quadratic, or None for fewer than four points).
    """
    s = s[(s["Nu"] - 1) >= min_q] if min_q > 0 else s
    x = (s["Ra"] * s["Pr"] ** p).to_numpy(float)
    y = (s["Nu"] - 1).to_numpy(float)
    sy = s["Nu_se"].to_numpy(float) if "Nu_se" in s else None
    if sy is not None and not np.all(np.isfinite(sy)):
        sy = None
    lx, ly = np.log10(x), np.log10(y)
    d = np.diff(lx)
    gamma = np.diff(ly) / d
    gerr = None
    if sy is not None:
        sl = sy / (y * LN10)                       # error of log10(Nu-1)
        gerr = np.sqrt(sl[1:] ** 2 + sl[:-1] ** 2) / d
    fit = fit_loglog(x, y, sy, "quadratic") if len(x) >= 4 else None
    return dict(x=x, y=y, sy=sy, x_mid=10 ** (0.5 * (lx[1:] + lx[:-1])), gamma=gamma, gamma_err=gerr, fit=fit)


def global_exponents(df):
    """alpha, beta of one power law Nu - 1 = C Pr^alpha Ra^beta fitted to all points.

    Parameters
    ----------
    df : pandas.DataFrame
        Results rows.

    Returns
    -------
    tuple of float
        ``(alpha, beta)``.
    """
    from rbc.analysis import exponents as ex
    s = ex.Surface(ex.prepare(df, "Nu-1"), 1)
    return float(s.exponent(0, 0, "alpha")[0]), float(s.exponent(0, 0, "beta")[0])


def family(ax, axg, df, a, p, min_points, min_q, ref_lines):
    """One family of slices: Nu - 1 against x (``ax``) and gamma against x (``axg``).

    Parameters
    ----------
    ax, axg : matplotlib.axes.Axes
        Axes for Nu - 1 and for the local exponent.
    df : pandas.DataFrame
        Results rows (finished).
    a : float
        Slice slope (Pr ∝ Ra^a).
    p : float
        Power of Pr in x = Ra Pr^p.
    min_points : int
        Fewest points for a slice.
    min_q : float
        Leave out points with Nu - 1 below this.
    ref_lines : dict
        Label -> gamma, drawn as horizontal lines on ``axg``.

    Returns
    -------
    list of dict
        One row per secant, for the CSV.
    """
    rows, handles, labels, xs = [], [], [], []
    for c, s in rc.slices(df, a, min_points).items():
        r = slice_exponents(s, p, min_q)
        if len(r["x"]) < 2:
            continue
        xs += list(r["x"])
        kw = style.slice_kw(c, a)
        line = style.slice_kw(c, a, line=True)
        ax.errorbar(r["x"], r["y"], yerr=r["sy"], capsize=1.5, **kw)
        ax.plot(r["x"], r["y"], lw=0.8, **line)
        axg.errorbar(r["x_mid"], r["gamma"], yerr=r["gamma_err"], capsize=1.5, **kw)
        axg.plot(r["x_mid"], r["gamma"], lw=0.8, **line)
        if r["fit"] is not None:
            q2, m, _ = r["fit"]["coef"]
            xx = np.logspace(np.log10(r["x"][0]), np.log10(r["x"][-1]), 60)
            axg.plot(xx, m + 2 * q2 * np.log10(xx), color=line["color"], lw=0.6, ls=":")
        handles.append(proxy(kw, line["ls"]))
        labels.append(f"${style.pow10(c)}$ ({len(r['x'])})")
        rows += [dict(a=a, p=p, c=c, x_mid=xm, gamma=g, gamma_err=np.nan if r["gamma_err"] is None else ge)
                 for xm, g, ge in zip(r["x_mid"], r["gamma"],
                                      r["gamma"] * np.nan if r["gamma_err"] is None else r["gamma_err"])]
    for (lab, gv), ls in zip(ref_lines.items(), ("--", "-.", ":")):
        axg.axhline(gv, color="0.45", lw=0.8, ls=ls, zorder=0, label=lab)
    for axis in (ax, axg):
        axis.set_xscale("log")
    if xs:
        sparse_decades(axg, min(xs), max(xs))
    ax.set_yscale("log")
    plain_log_y(ax)
    ax.set_ylabel(r"$Nu-1$")
    ax.set_title(FAMILY_TITLE.get(a, rf"$Pr\propto Ra^{{{a:g}}}$"), fontsize=9)
    axg.set_xlabel(x_label(p))
    axg.set_ylabel(r"local $\gamma = d\log(Nu-1)/d\log x$", fontsize=8)
    ref_h, ref_l = axg.get_legend_handles_labels()
    axg.legend(handles + ref_h, labels + ref_l,
               title=r"slice $Pr_0$ (points); reference lines", fontsize=6, title_fontsize=6, loc="upper center",
               bbox_to_anchor=(0.5, -0.3), ncol=2, handlelength=2.2, alignment="left")
    return rows


def get_parser() -> argparse.ArgumentParser:
    """Command-line options (see the module docstring).

    Returns
    -------
    argparse.ArgumentParser
    """
    p = argparse.ArgumentParser(prog="python -m rbc.plots.reduced", description=__doc__.split("\n")[0])
    p.add_argument("--family", type=float, nargs="+", default=[-0.5, -1.0], help="slice slopes a (one column each)")
    p.add_argument("--pr-power", type=float, help="p in x = Ra Pr^p for every family (default: p = a)")
    p.add_argument("--min-points", type=int, default=3, help="fewest points for a slice")
    p.add_argument("--min-q", type=float, default=0.0, help="leave out points with Nu-1 below this (weak convection)")
    p.add_argument("--no-fit-line", action="store_true", help="leave out the line from one power law fitted to all points")
    p.add_argument("--include-running", action="store_true", help="also use points whose status is not 'finished'")
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
    rc.require_columns(df, ["Nu"], "the reduced-variable figure")
    df = df[np.isfinite(df["Nu"])]
    if not a.include_running:
        df = df[df["status"] == "finished"]
    alpha, beta = (None, None) if a.no_fit_line else global_exponents(df)
    if alpha is not None:
        print(f"one power law for all points: Nu-1 ~ Pr^{alpha:.4f} Ra^{beta:.4f}")
    n = len(a.family)
    w, _ = style.figsize(1.0)
    fig, axes = plt.subplots(2, n, figsize=(w if n > 1 else 0.7 * w, 1.2 * w if n > 1 else 1.4 * w),
                             squeeze=False, sharex="col", layout="constrained", height_ratios=[1.2, 1])
    rows = []
    for k, fam in enumerate(a.family):
        p = fam if a.pr_power is None else a.pr_power
        span = 1 + p * fam                          # d log x / d log Ra along a slice
        refs = {rf"SL, GL IV$_l$: {(1 + fam) / (2 * span):.3g}": (1 + fam) / (2 * span)}
        if alpha is not None:
            g = (beta + fam * alpha) / span
            refs[rf"one law for all points: {g:.3g}"] = g
        rows += family(axes[0, k], axes[1, k], df, fam, p, a.min_points, a.min_q, refs)
    if not rows:
        raise SystemExit(f"no slice with {a.min_points} points or more in families {a.family}")
    def num(v):                                     # file-name form: -0.5 -> m0p5
        return f"{v:g}".replace("-", "m").replace(".", "p")

    tag = "reduced_" + "_".join(num(f) for f in a.family) + ("" if a.pr_power is None else f"_p{num(a.pr_power)}")
    print(style.save(fig, a.out / tag, png=a.png))
    pd.DataFrame(rows).to_csv(a.out / f"{tag}_gamma.csv", index=False, float_format="%.6g")


if __name__ == "__main__":
    main()
