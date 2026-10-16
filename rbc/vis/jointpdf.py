"""How alike are w and theta?  (V5)

(a) joint PDF of (w/w_rms, theta/theta_rms) on one plane, accumulated over snapshots.
    Mass along the diagonal (warm fluid rising, cool sinking) is what carries heat;
    a round blob means w and theta are unrelated.
(b) correlation coefficient C(z) = <w theta>_h / sqrt(<w^2>_h <theta^2>_h), theta = T - <T>_h,
    averaged over snapshots. C = 1: theta is a scaled copy of w. Since
    <w theta>_h = (Nu - 1)/sqrt(Ra Pr) away from the walls, a falling C at low Pr with
    fixed flux means w has variance that theta does not share: the extra small scales.

Usage::

    python -m rbc.vis.jointpdf RUN_PR1 RUN_PR001 --height mid --last 20 --jobs 8
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np

from rbc import style
from rbc.vis import fields as fl


def one_write(run, k, height, bins, lim, zchunk=16):
    """Histogram and correlation profile of one snapshot.

    Parameters
    ----------
    run : pathlib.Path
        Run directory.
    k : int
        Snapshot index.
    height : str or float
        Plane for the histogram.
    bins : int
        Bins per axis.
    lim : float
        Histogram range in rms units.
    zchunk : int, optional
        Levels read at a time for C(z) (memory).

    Returns
    -------
    H : numpy.ndarray
        (bins, bins) counts of (w, theta).
    z : numpy.ndarray
        Grid heights.
    C : numpy.ndarray
        Correlation coefficient at each height.
    """
    info = fl.RunInfo.load(run)
    wr = fl.writes(info.path)[k]
    with h5py.File(wr.file, "r") as h:
        z0 = info.height(height)
        th = fl.normalized(info, "T", fl.plane(h, "T", wr.index, z0), z0)
        w = fl.normalized(info, "w", fl.plane(h, "w", wr.index, z0), z0)
        H, _, _ = np.histogram2d(w.ravel(), th.ravel(), bins=bins, range=[[-lim, lim], [-lim, lim]])
        dT, cT = fl._task(h, "T")
        dw, cw = fl._task(h, "w")
        z = fl.grids(h)[2]
        C = np.empty(z.size)
        for a in range(0, z.size, zchunk):
            s = slice(a, min(a + zchunk, z.size))
            T = dT[(wr.index,) + ((cT,) if cT is not None else ()) + (slice(None), slice(None), s)]
            W = dw[(wr.index,) + ((cw,) if cw is not None else ()) + (slice(None), slice(None), s)]
            th3 = T - T.mean(axis=(0, 1))
            C[s] = (W * th3).mean(axis=(0, 1)) / np.sqrt((W ** 2).mean(axis=(0, 1)) * (th3 ** 2).mean(axis=(0, 1)))
    return H, z, C


def accumulate(run, height, last, bins, lim, jobs):
    """Joint PDF and mean C(z) over many snapshots, in parallel.

    Parameters
    ----------
    run : pathlib.Path
        Run directory.
    height : str or float
        Plane for the PDF.
    last : int
        Use the last N snapshots (0: all).
    bins : int
        Bins per axis.
    lim : float
        Range in rms units.
    jobs : int
        Parallel processes.

    Returns
    -------
    tuple
        ``(P, z, C, n)``: normalized histogram, heights, mean C(z), snapshots used.
    """
    n = len(fl.writes(Path(run).resolve()))
    ks = range(max(0, n - last), n) if last else range(n)
    with ProcessPoolExecutor(jobs) as ex:
        res = list(ex.map(one_write, *zip(*[(run, k, height, bins, lim) for k in ks])))
    H = sum(r[0] for r in res)
    return H / H.sum(), res[0][1], np.mean([r[2] for r in res], axis=0), len(res)


def get_parser() -> argparse.ArgumentParser:
    """Command-line options (see the module docstring).

    Returns
    -------
    argparse.ArgumentParser
    """
    p = argparse.ArgumentParser(prog="python -m rbc.vis.jointpdf", description=__doc__.split("\n")[0])
    p.add_argument("runs", type=Path, nargs="+", help="run directories (one panel each)")
    p.add_argument("--height", default="mid", help="mid, bl (1/(4 Nu)), bl-edge (1/(2 Nu)) or a number")
    p.add_argument("--last", type=int, default=0, help="use the last N snapshots (0 = all)")
    p.add_argument("--bins", type=int, default=121, help="histogram bins per axis")
    p.add_argument("--lim", type=float, default=4.0, help="axis limits in rms units")
    p.add_argument("--jobs", type=int, default=4, help="parallel processes")
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

    n = len(a.runs)
    fig = plt.figure(figsize=style.figsize(1.0, 0.42 if n <= 2 else 0.34))
    gs = fig.add_gridspec(1, n + 1, width_ratios=[1] * n + [1.1])
    axc = fig.add_subplot(gs[0, n])
    from matplotlib import colormaps
    cm = colormaps["viridis"]
    edges = np.linspace(-a.lim, a.lim, a.bins + 1)
    mid = 0.5 * (edges[1:] + edges[:-1])
    levels = np.arange(-5, 0.01, 0.5)
    for i, run in enumerate(a.runs):
        info = fl.RunInfo.load(run)
        P, z, C, nw = accumulate(run, a.height, a.last, a.bins, a.lim, a.jobs)
        dens = P / np.diff(edges)[0] ** 2
        ax = fig.add_subplot(gs[0, i])
        with np.errstate(divide="ignore"):
            cf = ax.contourf(mid, mid, np.log10(dens).T, levels=levels, cmap="Greys", extend="min")
        ax.plot([-a.lim, a.lim], [-a.lim, a.lim], color="0.5", lw=0.5, ls="--")
        z0 = info.height(a.height)
        Cz = np.interp(z0, z, C)
        t = [-int(a.lim) + 1, 0, int(a.lim) - 1]
        ax.set(xlabel=r"$w/w_{\rm rms}$", aspect="equal", xticks=t, yticks=t,
               title=f"{info.label()}\n$C={Cz:.2f}$ ({nw} snapshots)")
        ax.title.set_fontsize(7)
        if i == 0:
            ax.set_ylabel(r"$\theta/\theta_{\rm rms}$")
        col = cm(0.85 * i / max(n - 1, 1))
        axc.plot(C, z, color=col, label=info.label())
        axc.plot([Cz], [z0], "o", color=col, ms=3)
    pdf_axes = [ax for ax in fig.axes if ax is not axc]
    cb = fig.colorbar(cf, ax=pdf_axes, location="right", shrink=0.8, pad=0.02, aspect=20)
    cb.set_label(r"$\log_{10}$ PDF", fontsize=7)
    cb.ax.tick_params(labelsize=6)
    axc.set(xlabel=r"$C(z)$", ylabel=r"$z$", xlim=(0, 1), ylim=(0, 1))
    axc.legend(fontsize=6)
    print(style.save(fig, a.out / f"jointpdf_{a.height}", png=a.png))


if __name__ == "__main__":
    main()
