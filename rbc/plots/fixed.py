"""Nu-1 (or Nu, Re) against one parameter with the other held fixed (figure D).

Usage::

    python -m rbc.plots.fixed                    # Nu-1 vs Pr at fixed Ra | Nu-1 vs Ra at fixed Pr
    python -m rbc.plots.fixed --q Nu --along ra --lit data/literature.csv
    python -m rbc.plots.fixed --slopes           # fitted exponent of each curve in the legend

Nu-1 is the convective part of the heat flux, the quantity the exact relations
<eps_u> = (Nu-1)/sqrt(Ra Pr) and the scaling theories refer to; at Nu ~ 3 the local
exponent of Nu is (Nu-1)/Nu = 2/3 of that of Nu-1, so the choice matters at low Ra.

--lit: optional CSV with columns source,Ra,Pr,Nu[,Re] (other studies) drawn in grey,
one marker per source, for checking agreement with published data.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import colormaps

from rbc import cases as rc
from rbc import style
from rbc.plots.common import decade_ticks, fit_loglog, tidy_log_y

#: Markers cycled through the literature sources.
LIT_MARKERS = "osD^v<>PX*"


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


def panel(ax, df, q, along, lit=None, slopes=False, leg_ax=None):
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
        Add each curve's fitted exponent (a single power law through all its points) to
        its legend entry.
    leg_ax : matplotlib.axes.Axes, optional
        Empty axes beside ``ax`` that holds the legend as one column (default: inside
        ``ax``).
    """
    other, x = ("log_ra", "Pr") if along == "pr" else ("log_pr", "Ra")
    vals = sorted(df[other].unique())
    cmap = colormaps["viridis"]
    lo, hi = min(vals), max(vals)
    sym = r"\alpha" if along == "pr" else r"\beta"
    for v in vals:
        g = df[df[other] == v].sort_values(x)
        if len(g) < 2:
            continue
        col = cmap(0.9 * (v - lo) / max(hi - lo, 1e-9))
        lab = rf"${style.pow10(v)}$"
        y, err = qcol(g, q)
        if slopes:      # value first, so the numbers line up in the legend column
            f = fit_loglog(g[x].values, y.values, None if err is None else err.values)
            lab = rf"${sym} = {f['coef'][-2]:.3f}$ at {lab}"
        ax.errorbar(g[x], y, yerr=err, color=col, marker="o", ms=3, lw=0.9, capsize=1.5, label=lab)
    if lit is not None:
        for i, (src, g) in enumerate(lit.groupby("source")):
            y, _ = qcol(g, q) if (q != "Nu-1" or "Nu" in g) else (None, None)
            if y is not None:
                ax.plot(g[x], y, ls="none", marker=LIT_MARKERS[i % len(LIT_MARKERS)], mfc="none",
                        mec="0.45", ms=4, label=src)
    ax.set_xscale("log")
    ax.set_yscale("log")
    decade_ticks(ax, "x")
    tidy_log_y(ax)
    ax.set_xlabel(style.SYM["Pr" if along == "pr" else "Ra"])
    ax.set_ylabel(r"$Nu-1$" if q == "Nu-1" else style.SYM[q])
    fixed = "Ra" if along == "pr" else "Pr"
    title = rf"fixed ${fixed}$" + (r" (fitted slope)" if slopes else "")
    h, l = ax.get_legend_handles_labels()
    if leg_ax is None:
        ax.legend(h, l, title=title, fontsize=6.5, title_fontsize=7, loc="best")
    else:
        leg_ax.axis("off")
        leg_ax.legend(h, l, title=title, loc="center left", borderaxespad=0, fontsize=7,
                      title_fontsize=7, handlelength=1.4, handletextpad=0.5, labelspacing=0.35,
                      alignment="left")


def get_parser() -> argparse.ArgumentParser:
    """Command-line options (see the module docstring).

    Returns
    -------
    argparse.ArgumentParser
    """
    p = argparse.ArgumentParser(prog="python -m rbc.plots.fixed", description=__doc__.split("\n")[0])
    p.add_argument("--q", choices=("Nu-1", "Nu", "Re"), default="Nu-1", help="quantity on the y axis")
    p.add_argument("--along", choices=("pr", "ra"), nargs="+", default=["pr", "ra"], help="pr: against Pr at fixed Ra; ra: against Ra at fixed Pr")
    p.add_argument("--slopes", action="store_true", help="put each curve's fitted exponent in the legend")
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
    lit = pd.read_csv(a.lit) if a.lit else None
    n = len(a.along)
    w, _ = style.figsize(1.0)
    fig = plt.figure(figsize=(w, 0.37 * w * n + 0.1))
    gs = fig.add_gridspec(n, 2, width_ratios=[2.5, 1])
    for k, along in enumerate(a.along):
        panel(fig.add_subplot(gs[k, 0]), df, a.q, along, lit, a.slopes, leg_ax=fig.add_subplot(gs[k, 1]))
    print(style.save(fig, a.out / f"fixed_{a.q.replace('-', 'm')}_{'_'.join(a.along)}", png=a.png))


if __name__ == "__main__":
    main()