"""Filled contours of Nu (or Re, or compensated Nu) over the (Ra, Pr) plane (figure I).

Interpolates log10 Q over a Delaunay triangulation in (log Ra, log Pr), so it works
on the irregular set of points actually simulated and never extrapolates past
their convex hull (no NaN -> 1 filling). Refined with a C1 cubic interpolant for
smooth contours with --smooth (cubic; can overshoot where points are sparse).
Flat, sliver triangles on the hull edge are masked. Usage::

    python -m rbc.plots.contour
    python -m rbc.plots.contour --q Nu --comp 0.25 --levels 10
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
from matplotlib import ticker

from rbc import cases as rc
from rbc import style


def get_parser() -> argparse.ArgumentParser:
    """Command-line options (see the module docstring).

    Returns
    -------
    argparse.ArgumentParser
    """
    p = argparse.ArgumentParser(prog="python -m rbc.plots.contour", description=__doc__.split("\n")[0])
    p.add_argument("--q", choices=("Nu", "Re"), default="Nu", help="quantity to contour")
    p.add_argument("--comp", type=float, default=0.0, help="plot Q / Ra^comp")
    p.add_argument("--levels", type=int, default=12, help="number of filled levels")
    p.add_argument("--cmap", default="viridis", help="Matplotlib colormap")
    p.add_argument("--subdiv", type=int, default=4, help="refinement levels of the triangulation")
    p.add_argument("--smooth", action="store_true", help="cubic instead of linear interpolation")
    p.add_argument("--results", type=Path, default=rc.RESULTS_CSV, help="results table (default: data/results.csv)")
    p.add_argument("--out", type=Path, default=Path("figures"), help="output folder")
    p.add_argument("--width", type=float, default=0.62, help="figure width as a fraction of the text width")
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
    rc.require_columns(df, [a.q], "the contour figure")
    df = df[(df["status"] == "finished") & np.isfinite(df[a.q])]
    x, y = df["log_ra"].values.astype(float), df["log_pr"].values.astype(float)
    v = np.log10(df[a.q].values / df["Ra"].values ** a.comp)

    tri = mtri.Triangulation(x, y)
    tri.set_mask(mtri.TriAnalyzer(tri).get_flat_tri_mask(min_circle_ratio=0.05))
    interp = (mtri.CubicTriInterpolator(tri, v, kind="min_E") if a.smooth
              else mtri.LinearTriInterpolator(tri, v))
    fine, vf = mtri.UniformTriRefiner(tri).refine_field(v, triinterpolator=interp, subdiv=a.subdiv)
    vf = np.ma.masked_invalid(vf)
    bad = np.any(np.ma.getmaskarray(vf)[fine.triangles], axis=1)
    if fine.mask is not None:
        bad |= fine.mask
    fine_plot = mtri.Triangulation(10 ** fine.x, 10 ** fine.y, fine.triangles, mask=bad)
    levels = np.linspace(v.min(), v.max(), a.levels + 1)   # linear interpolation stays in range
    vf = vf.filled(np.nan)

    fig, ax = plt.subplots(figsize=style.figsize(a.width, 0.8))
    cf = ax.tricontourf(fine_plot, vf, levels=levels, cmap=a.cmap)
    cl = ax.tricontour(fine_plot, vf, levels=levels[::2], colors="white", linewidths=0.4, alpha=0.7)
    ax.clabel(cl, fmt=lambda L: f"{10 ** L:.3g}", fontsize=6, inline_spacing=3)

    for fam in rc.FAMILIES:
        for c, s in rc.slices(df, fam).items():
            ax.plot(s["Ra"], s["Pr"], color="white", lw=0.8, ls=style.FAMILY_STYLE[fam]["ls"], alpha=0.9)
    ax.plot(df["Ra"], df["Pr"], "o", ms=2.5, mfc="white", mec="k", mew=0.4)

    pad = 0.08
    ax.set(xscale="log", yscale="log", xlabel=style.SYM["Ra"], ylabel=style.SYM["Pr"],
           xlim=(10 ** (x.min() - pad), 10 ** (x.max() + pad)), ylim=(10 ** (y.min() - pad), 10 ** (y.max() + pad)))
    lab = a.q if not a.comp else rf"{a.q}\,/\,Ra^{{{a.comp:g}}}"
    cb = fig.colorbar(cf, ax=ax, pad=0.02, aspect=25)
    cb.set_label(rf"$\log_{{10}}\,{lab}$" if a.comp else rf"$\log_{{10}} {a.q}$")
    cb.ax.tick_params(direction="in", which="both")
    cb.ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.2f"))
    tag = a.q + (f"_comp{a.comp:g}" if a.comp else "")
    print(style.save(fig, a.out / f"contour_{tag}", png=a.png))


if __name__ == "__main__":
    main()