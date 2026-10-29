"""The scaling plot of the slices against x = Ra Pr^p instead of Ra, with its local exponent.

For slices Pr ∝ Ra^a (``--family``), every point is placed at x = Ra Pr^p, by default
with p = a: x = Ra Pr^-1/2 for the slices Pr ∝ Ra^-1/2, and x = Ra Pr^-1 = Ra/Pr = Re_ff^2
(the free-fall Reynolds number squared) for Pr ∝ Ra^-1. ``--pr-power`` sets p for every
family.

Top row: exactly :mod:`rbc.plots.scaling` against x. Nu - 1 is divided by the power of x
the reference law Nu - 1 ∝ Pr^alpha Ra^beta predicts along the slice (default the low-Pr
ultimate law of Shishkina & Lohse, alpha = beta = 1/2; ``--alpha``, ``--beta``), so the
prediction is a horizontal line, with the same straight-line (and ``--fit quadratic``)
fits and the same vertical normalizations (``--norm``). Along a slice
log x = (1 + p a) log Ra + const, so that power is (beta + a alpha) / (1 + p a): x^(1/5)
for the Pr ∝ Ra^-1/2 slices, and x^0 (no compensation) for Pr ∝ Ra^-1, where the
reference law predicts a constant because it depends on Ra Pr only.

Bottom row: the local exponent of Nu - 1 ~ c x^gamma, gamma = d log(Nu-1) / d log x, from
two-point secants between neighbouring points of each slice (drawn at the geometric mean of
their x), the derivative of a quadratic fit for slices of four points or more (dotted), and
each straight-line fit of the top row as a band (its slope plus the compensating power).

What gamma is: along a slice it is the slice exponent d log(Nu-1)/d log Ra divided by
(1 + p a), i.e. (beta + a alpha) / (1 + p a) in terms of the local exponents of
Nu - 1 ~ Pr^alpha Ra^beta. What the combined variable adds is a common horizontal axis:
if the gamma(x) of all slices fall on one curve, x organizes the approach to the
asymptotic exponent. Nu - 1 itself collapses only if alpha = p beta.

Usage::

    python -m rbc.plots.reduced                                  # both families
    python -m rbc.plots.reduced --family -1 --fit linear quadratic
    python -m rbc.plots.reduced --pr-power -1                    # x = Ra/Pr for every family
    python -m rbc.plots.reduced --norm common
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from rbc import cases as rc
from rbc import style
from rbc.plots.common import fit_loglog

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


def family(ax, axg, df, a, p, ref, fit_kinds, norm, ref_ra, min_points, include_running, ref_lines):
    """One family of slices: the compensated scaling plot against x (``ax``) and the local
    exponent gamma against x (``axg``).

    Parameters
    ----------
    ax, axg : matplotlib.axes.Axes
        Axes for compensated Nu - 1 and for the local exponent.
    df : pandas.DataFrame
        Results rows.
    a : float
        Slice slope (Pr ∝ Ra^a).
    p : float
        Power of Pr in x = Ra Pr^p.
    ref : tuple of float
        (alpha, beta) of the reference law Nu - 1 ∝ Pr^alpha Ra^beta used to compensate.
    fit_kinds : list of str
        ``"linear"`` and/or ``"quadratic"`` fits of each slice, as in :mod:`rbc.plots.scaling`.
    norm : str
        Vertical normalization of the top panel (:data:`rbc.plots.scaling.NORMS`).
    ref_ra : float
        log10 Ra for ``norm="ref"``.
    min_points : int
        Fewest points for a slice.
    include_running : bool
        Also use points whose status is not 'finished'.
    ref_lines : dict
        Label -> gamma, drawn as horizontal lines on ``axg``.

    Returns
    -------
    fits : list of dict
        One row per (slice, fit) of the top panel, slopes against x.
    secants : list of dict
        One row per secant of the lower panel.
    """
    from rbc.plots import scaling

    items = []
    fits = scaling.panel(ax, df, "Nu-1", a, ref[0], ref[1], fit_kinds, norm, ref_ra, min_points, include_running,
                         p=p, legend=False, legend_items=items)
    use = df if include_running else df[df["status"] == "finished"]
    secants, xs = [], []
    for c, s in rc.slices(use, a, min_points).items():
        r = slice_exponents(s, p)
        xs += list(r["x"])
        kw = style.slice_kw(c, a)
        line = style.slice_kw(c, a, line=True)
        axg.errorbar(r["x_mid"], r["gamma"], yerr=r["gamma_err"], capsize=1.5, **kw)
        axg.plot(r["x_mid"], r["gamma"], lw=0.8, **line)
        if r["fit"] is not None:
            q2, m, _ = r["fit"]["coef"]
            xx = np.logspace(np.log10(r["x"][0]), np.log10(r["x"][-1]), 60)
            axg.plot(xx, m + 2 * q2 * np.log10(xx), lw=0.6, ls=":", **style.line_kw(line["color"]))
        for f in fits:     # the straight-line fit's exponent in x, over the slice's range
            if f["c"] == c and f["kind"] == "linear":
                axg.plot([r["x"][0], r["x"][-1]], [f["slope"] + f["slice_exponent"]] * 2, color=line["color"],
                         lw=2.2, alpha=0.3, solid_capstyle="butt", zorder=1)
        secants += [dict(a=a, p=p, c=c, x_mid=xm, gamma=g, gamma_err=np.nan if r["gamma_err"] is None else ge)
                    for xm, g, ge in zip(r["x_mid"], r["gamma"],
                                         r["gamma"] * np.nan if r["gamma_err"] is None else r["gamma_err"])]
    for (lab, gv), ls in zip(ref_lines.items(), ("--", "-.", ":")):
        axg.axhline(gv, color="0.45", lw=0.8, ls=ls, zorder=0, label=lab)
    axg.set_xscale("log")
    if xs:
        sparse_decades(axg, min(xs), max(xs))
    axg.set_xlabel(x_label(p))
    axg.set_ylabel(r"local $\gamma = d\log(Nu-1)/d\log x$", fontsize=8)
    ref_h, ref_l = axg.get_legend_handles_labels()
    from matplotlib.lines import Line2D
    band = Line2D([], [], color="0.5", lw=2.2, alpha=0.3)
    handles = [h for h, _ in items] + [band] + ref_h
    labels = [lab for _, lab in items] + ["straight-line fit, as an exponent in $x$"] + ref_l
    axg.legend(handles, labels, title=r"slice $Pr_0$: fitted slope $m$ of the compensated data against $x$",
               fontsize=6, title_fontsize=6, loc="upper center", bbox_to_anchor=(0.5, -0.3), ncol=1,
               handlelength=2.2, alignment="left")
    return fits, secants


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
    p.add_argument("--fit", nargs="+", choices=("linear", "quadratic"), default=["linear"],
                   help="fits of log(Nu-1) against log x along each slice, as in rbc.plots.scaling")
    p.add_argument("--norm", choices=("none", "first", "ref", "common", "law"), default="none",
                   help="vertical normalization of the top row (see rbc.plots.scaling)")
    p.add_argument("--ref-ra", type=float, default=7.0, help="log10 Ra for --norm ref")
    p.add_argument("--alpha", type=float, default=0.5, help="Pr exponent of the reference law used to compensate")
    p.add_argument("--beta", type=float, default=0.5, help="Ra exponent of the reference law used to compensate")
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
    if a.min_q > 0:
        df = df[df["Nu"] - 1 >= a.min_q]
    fitted = None if a.no_fit_line else global_exponents(
        df if a.include_running else df[df["status"] == "finished"])
    if fitted is not None:
        print(f"one power law for all points: Nu-1 ~ Pr^{fitted[0]:.4f} Ra^{fitted[1]:.4f}")
    n = len(a.family)
    w, _ = style.figsize(1.0)
    fig, axes = plt.subplots(2, n, figsize=(w if n > 1 else 0.7 * w, 1.3 * w if n > 1 else 1.5 * w),
                             squeeze=False, sharex="col", layout="constrained", height_ratios=[1.2, 1])
    fits, secants = [], []
    for k, fam in enumerate(a.family):
        p = fam if a.pr_power is None else a.pr_power
        span = 1 + p * fam                          # d log x / d log Ra along a slice
        g_ref = (a.beta + fam * a.alpha) / span
        refs = {rf"reference law ($\alpha$={a.alpha:g}, $\beta$={a.beta:g}; SL, GL IV$_l$): {g_ref:.3g}": g_ref}
        if fitted is not None:
            refs[rf"one law for all points: {(fitted[1] + fam * fitted[0]) / span:.3g}"] = \
                (fitted[1] + fam * fitted[0]) / span
        f, s_ = family(axes[0, k], axes[1, k], df, fam, p, (a.alpha, a.beta), a.fit, a.norm, a.ref_ra,
                       a.min_points, a.include_running, refs)
        axes[0, k].set_xlabel("")
        fits += f
        secants += s_
    if not secants:
        raise SystemExit(f"no slice with {a.min_points} points or more in families {a.family}")

    def num(v):                                     # file-name form: -0.5 -> m0p5
        return f"{v:g}".replace("-", "m").replace(".", "p")

    tag = "reduced_" + "_".join(num(f) for f in a.family) + ("" if a.pr_power is None else f"_p{num(a.pr_power)}")
    print(style.save(fig, a.out / tag, png=a.png))
    pd.DataFrame(secants).to_csv(a.out / f"{tag}_gamma.csv", index=False, float_format="%.6g")
    fits = pd.DataFrame(fits)
    fits.to_csv(a.out / f"{tag}_fits.csv", index=False, float_format="%.6g")
    with pd.option_context("display.width", 140, "display.precision", 4):
        print(fits)


if __name__ == "__main__":
    main()
