"""Map of the (Ra, Pr) points simulated: filled = finished, open = running,
star = resolution study, coloured lines = slices Pr ∝ Ra^{-1/2} (dashed Pr ∝ Ra^{-1}
with --ra-inv).

Edit ``data/cases.csv`` to move or add points; nothing here needs to change. Usage::

    python -m rbc.plots.param_grid                 # 2d and 3d
    python -m rbc.plots.param_grid --dim 3 --legend --png
    python -m rbc.plots.param_grid --ra-inv        # also draw Pr ∝ Ra^-1 slices
    python -m rbc.plots.param_grid --span finished # lines stop at the last finished point
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from rbc import cases as rc
from rbc import style

#: Marker size of the circles.
MS = 6.5


def _limits(log_vals, pad):
    """Axis limits padded by ``pad`` decades around log10 values."""
    return 10 ** (np.min(log_vals) - pad), 10 ** (np.max(log_vals) + pad)


def draw(ax, df, families=(-0.5,), span="all", min_points=3, ms=MS, count="all"):
    """Draw the parameter map on an axes.

    Slices are every line in the table with at least ``min_points`` points.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes.
    df : pandas.DataFrame
        Cases of one dimension (:func:`rbc.cases.load_cases`).
    families : tuple of float, optional
        Slice slopes to draw (-0.5 solid, -1.0 dashed).
    span : {"all", "finished"}, optional
        A slice line runs through every point on it, or stops at the last finished one.
    min_points : int, optional
        Fewest points for a slice.
    ms : float, optional
        Marker size.
    count : {"all", "finished"}, optional
        Which points count towards ``min_points``.

    Returns
    -------
    dict
        Legend label -> handle, one per status (and the resolution stars).
    """
    base = df if count == "all" else df[df["status"] == "finished"]
    for a in families:
        for c, on in rc.slices(base, a, min_points).items():
            if span == "finished":
                on = on[on["status"] == "finished"]
            if len(on) < 2:
                continue
            lr = np.array([on["log_ra"].min(), on["log_ra"].max()])
            ax.plot(10 ** lr, 10 ** (c + a * (lr - rc.RA_REF)), lw=1.0, zorder=1,
                    **style.slice_kw(c, a, line=True))

    handles = {}
    for status, kw in style.STATUS_STYLE.items():
        for star in (False, True):
            sub = df[(df["status"] == status) & (df["res_check"] == star)]
            if sub.empty:
                continue
            k = dict(kw, ls="none", ms=ms, mew=0.8, zorder=3)
            if star:
                k.update(marker="*", ms=ms * style.RES_CHECK_SCALE * 1.15)
            label = "Resolution study" if star else status.capitalize()
            h, = ax.plot(sub["Ra"], sub["Pr"], **k)
            handles.setdefault(label, h)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(style.SYM["Ra"])
    ax.set_ylabel(style.SYM["Pr"], rotation=0, labelpad=8)
    return handles


def get_parser() -> argparse.ArgumentParser:
    """Command-line options (see the module docstring).

    Returns
    -------
    argparse.ArgumentParser
    """
    p = argparse.ArgumentParser(prog="python -m rbc.plots.param_grid", description=__doc__.split("\n")[0])
    p.add_argument("--dim", type=int, choices=(2, 3), nargs="*", default=[2, 3], help="which dimensions to draw (one figure each)")
    p.add_argument("--csv", type=Path, default=rc.CASES_CSV, help="the cases registry")
    p.add_argument("--out", type=Path, default=Path("figures"), help="output folder")
    p.add_argument("--span", choices=("all", "finished"), default="all", help="slice lines run through all points, or stop at the last finished one")
    p.add_argument("--ra-inv", action="store_true", help="also draw slices Pr ∝ Ra^-1")
    p.add_argument("--min-points", type=int, default=3, help="fewest points for a slice to be drawn")
    p.add_argument("--count", choices=("all", "finished"), default="all",
                   help="which points count towards a slice's min-points")
    p.add_argument("--width", type=float, default=0.62, help="fraction of text width")
    p.add_argument("--pad", type=float, default=0.25, help="axis padding in decades")
    p.add_argument("--legend", action="store_true", help="add a legend above the axes")
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
    allcases = rc.load_cases(a.csv)
    ylim = _limits(allcases["log_pr"], a.pad)   # shared so 2d/3d panels line up

    for dim in a.dim:
        df = allcases[allcases["dim"] == dim]
        fig, ax = plt.subplots(figsize=style.figsize(a.width, 0.75))
        families = (-0.5, -1.0) if a.ra_inv else (-0.5,)
        h = draw(ax, df, families, a.span, a.min_points, count=a.count)
        ax.set_xlim(_limits(df["log_ra"], a.pad))
        ax.set_ylim(ylim)
        if a.legend:
            ax.legend(h.values(), h.keys(), loc="lower center", bbox_to_anchor=(0.5, 1.0),
                      ncol=len(h), numpoints=1, markerscale=0.8, columnspacing=1.0,
                      handletextpad=0.3)
        print(style.save(fig, a.out / f"param_grid_{dim}d", png=a.png))


if __name__ == "__main__":
    main()
