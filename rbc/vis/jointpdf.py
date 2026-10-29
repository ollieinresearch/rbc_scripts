"""How alike are w and theta?  (V5)

(a) joint PDF of (w/w_rms, theta/theta_rms) on one plane, accumulated over snapshots.
    Mass along the diagonal (warm fluid rising, cool sinking) is what carries heat;
    a round blob means w and theta are unrelated.
(b) correlation coefficient C(z) = <w theta>_h / sqrt(<w^2>_h <theta^2>_h), theta = T - <T>_h,
    averaged over snapshots. C = 1: theta is a scaled copy of w. Since
    <w theta>_h = (Nu - 1)/sqrt(Ra Pr) away from the walls, a falling C at low Pr with
    fixed flux means w has variance that theta does not share: the extra small scales.

Usage::

    python -m rbc.vis.jointpdf RUN_A RUN_B RUN_C --height mid --last 20

Any number of runs, one panel each. The snapshots of all runs are processed by one pool of
workers, as many as the cores and memory allow (``--jobs 0``, see :mod:`rbc.vis.parallel`).
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
from rbc.vis import parallel

#: Peak memory of one snapshot: a fixed part, and bytes per horizontal grid point (the
#: 8-level plane reads and 16-level chunks of T and w for C(z), in float64, with copies).
TASK_BASE = 0.4 * parallel.GB
BYTES_PER_HPOINT = 1000


def one_write(info, wr, height, bins, lim, zchunk=16):
    """Histogram and correlation profile of one snapshot.

    Parameters
    ----------
    info : rbc.vis.fields.RunInfo
        The run.
    wr : rbc.vis.fields.Write
        The snapshot write.
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


def _one(args):
    """Worker: :func:`one_write` on a tuple of arguments."""
    return one_write(*args)


def accumulate(infos, height, last, bins, lim, jobs):
    """Joint PDF and mean C(z) of every run, over many snapshots, in one pool.

    Parameters
    ----------
    infos : list of rbc.vis.fields.RunInfo
        The runs.
    height : str or float
        Plane for the PDF.
    last : int
        Use the last N snapshots of each run (0: all).
    bins : int
        Bins per axis.
    lim : float
        Range in rms units.
    jobs : int
        Parallel processes (0: as many as cores and memory allow).

    Returns
    -------
    list of tuple
        Per run ``(P, z, C, n)``: normalized histogram, heights, mean C(z), snapshots used.
    """
    tasks, owner = [], []
    for i, info in enumerate(infos):
        ws = fl.writes(info.path)
        if not ws:
            raise SystemExit(f"{info.path}: no snapshots")
        for wr in (ws[-last:] if last else ws):
            tasks.append((info, wr, height, bins, lim))
            owner.append(i)
    hpoints = max(np.prod(fl.grid_shape(t[1])[:2]) for t in tasks)
    jobs = parallel.auto_jobs(len(tasks), TASK_BASE + BYTES_PER_HPOINT * hpoints, jobs, what="snapshots")
    with ProcessPoolExecutor(jobs) as ex:
        res = list(ex.map(_one, tasks, chunksize=max(1, len(tasks) // (4 * jobs))))
    out = []
    for i in range(len(infos)):
        r = [x for x, o in zip(res, owner) if o == i]
        H = sum(x[0] for x in r)
        out.append((H / H.sum(), r[0][1], np.mean([x[2] for x in r], axis=0), len(r)))
    return out


def get_parser() -> argparse.ArgumentParser:
    """Command-line options (see the module docstring).

    Returns
    -------
    argparse.ArgumentParser
    """
    p = argparse.ArgumentParser(prog="python -m rbc.vis.jointpdf", description=__doc__.split("\n")[0])
    p.add_argument("runs", type=Path, nargs="+", help="run directories, any number (one panel each)")
    p.add_argument("--height", default="mid", help="mid, bl (1/(4 Nu)), bl-edge (1/(2 Nu)) or a number")
    p.add_argument("--last", type=int, default=0, help="use the last N snapshots (0 = all)")
    p.add_argument("--bins", type=int, default=121, help="histogram bins per axis")
    p.add_argument("--lim", type=float, default=4.0, help="axis limits in rms units")
    p.add_argument("--labels", nargs="+", help="panel labels, one per run (default: Ra, Pr and, if needed, N_z)")
    p.add_argument("--columns", type=int, default=3, help="PDF panels per row (more runs wrap onto further rows)")
    p.add_argument("--jobs", type=int, default=0, help="parallel processes (0 = as many as cores and memory allow)")
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
    from rbc.vis.slices import row_labels

    a = get_parser().parse_args(argv)
    style.use("paper", usetex=False if a.no_tex else None)
    infos = [fl.RunInfo.load(r) for r in a.runs]
    labels = row_labels(infos, a.labels)
    results = accumulate(infos, a.height, a.last, a.bins, a.lim, a.jobs)

    n = len(infos)
    ncol = min(n, a.columns)
    nrow = -(-n // ncol)
    w_in, _ = style.figsize(1.0)
    panel = w_in / (ncol + 1.3)
    fig = plt.figure(figsize=(w_in, nrow * (panel * 1.25 + 0.15) + 0.55), layout="constrained")
    gs = fig.add_gridspec(nrow, ncol + 1, width_ratios=[1] * ncol + [1.2])
    axc = fig.add_subplot(gs[:, ncol])
    # runs ranked by expense: free-fall Reynolds number squared Ra/Pr, then N_z
    order = sorted(range(len(infos)), key=lambda i: (infos[i].Ra / infos[i].Pr,
                                                      (fl.parse_run_path(infos[i].path) or {}).get("nz") or 0))
    ranked = style.ranked_colours(len(infos))
    colour = {i: ranked[k] for k, i in enumerate(order)}
    edges = np.linspace(-a.lim, a.lim, a.bins + 1)
    mid = 0.5 * (edges[1:] + edges[:-1])
    levels = np.arange(-5, 0.01, 0.5)
    t = [-int(a.lim) + 1, 0, int(a.lim) - 1]
    pdf_axes = []
    for i, (info, label, (P, z, C, nw)) in enumerate(zip(infos, labels, results)):
        dens = P / np.diff(edges)[0] ** 2
        ax = fig.add_subplot(gs[i // ncol, i % ncol])
        pdf_axes.append(ax)
        with np.errstate(divide="ignore"):
            cf = ax.contourf(mid, mid, np.log10(dens).T, levels=levels, cmap="Greys", extend="min")
        ax.plot([-a.lim, a.lim], [-a.lim, a.lim], color="0.5", lw=0.5, ls="--")
        z0 = info.height(a.height)
        Cz = np.interp(z0, z, C)
        col = colour[i]
        ax.set(aspect="equal", xticks=t, yticks=t)
        ax.set_title(f"{label}\n$C={Cz:.2f}$, {nw} snapshots", fontsize=6, color=col if n > 1 else "k")
        ax.tick_params(labelsize=6)
        if i % ncol == 0:
            ax.set_ylabel(r"$\theta/\theta_{\rm rms}$", fontsize=7)
        if i >= n - ncol:
            ax.set_xlabel(r"$w/w_{\rm rms}$", fontsize=7)
        axc.plot(C, z, color=col, label=label.replace("\n", ", "))
        axc.plot([Cz], [z0], "o", color=col, ms=3)
    cb = fig.colorbar(cf, ax=pdf_axes, location="bottom", shrink=0.6, aspect=30)
    cb.set_label(r"$\log_{10}$ PDF", fontsize=7)
    cb.ax.tick_params(labelsize=6)
    axc.set(xlim=(0, 1), ylim=(0, 1))
    axc.set_xlabel(r"$C(z)$", fontsize=7)
    axc.set_ylabel(r"$z$", fontsize=7)
    axc.tick_params(labelsize=6)
    axc.legend(fontsize=5, loc="upper center", bbox_to_anchor=(0.5, -0.3), frameon=False)
    print(style.save(fig, a.out / f"jointpdf_{a.height}", png=a.png))


if __name__ == "__main__":
    main()
