"""Horizontal slices of several runs side by side (V1, V2).

Rows are runs, columns are (field, height) pairs, so midplane and boundary-layer
slices can share one figure. Colour limits are fixed in normalized units (see
:mod:`rbc.vis.fields`), so rows at different Pr are directly comparable and a movie
never rescales. Usage::

    python -m rbc.vis.slices RUN_PR1 RUN_PR001 --fields T w --heights mid bl
    python -m rbc.vis.slices RUN_A RUN_B --movie --jobs 16

Heights: mid (z = 1/2), bl (z = delta_T/2 = 1/(4 Nu)), bl-edge (z = delta_T), or a number.
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np

from rbc import style
from rbc.vis import fields as fl

CMAP = {"T": "RdBu_r", "w": "PuOr_r", "u": "PuOr_r", "v": "PuOr_r", "vort": "magma"}
LABEL = {"T": r"$\theta/\theta_{\rm rms}$", "w": r"$w/w_{\rm rms}$", "u": r"$u/u_{\rm rms}$",
         "v": r"$v/v_{\rm rms}$", "vort": r"$\log_{10}|\omega|/\omega_{\rm rms}$"}
HEIGHT_LABEL = {"mid": r"$z=1/2$", "bl": r"$z=\delta_T/2$", "bl-edge": r"$z=\delta_T$"}


def frame(runs, fields_, heights, k, lim, out_png=None, width=1.0, keep_negative=False):
    """One figure: rows are runs, columns are (field, height) pairs.

    Parameters
    ----------
    runs : list of pathlib.Path
        Run directories.
    fields_ : list of str
        Fields (keys of ``CMAP``).
    heights : list of str or float
        Heights (see :meth:`rbc.vis.fields.RunInfo.height`).
    k : int
        Snapshot index in each run (-1: last).
    lim : float
        Colour limit in rms units (vorticity uses a fixed log range).
    out_png : pathlib.Path, optional
        Save here and close the figure.
    width : float, optional
        Fraction of the text width.
    keep_negative : bool, optional
        Keep negative ``|ω|²`` points.

    Returns
    -------
    matplotlib.figure.Figure
    """
    infos = [fl.RunInfo.load(r) for r in runs]
    cols = [(f, h) for h in heights for f in fields_]
    nr, nc = len(runs), len(cols)
    w_in, _ = style.figsize(width)
    fig, axes = plt.subplots(nr, nc, figsize=(w_in, w_in / (nc + 0.6) * nr * 1.08 + 0.25), squeeze=False,
                             layout="constrained")
    ims = {}
    times = []
    for i, info in enumerate(infos):
        wr = fl.writes(info.path)[k]
        times.append(wr.time)
        with h5py.File(wr.file, "r") as h:
            x, y, _ = fl.grids(h)
            for j, (f, hgt) in enumerate(cols):
                z0 = info.height(hgt)
                d = fl.normalized(info, f, fl.plane(h, f, wr.index, z0), z0, keep_negative)
                if f == "vort":
                    d = np.log10(np.clip(d, 1e-3, None))
                    vmin, vmax = -1.0, 1.5
                else:
                    vmin, vmax = -lim, lim
                ax = axes[i, j]
                ims[f] = ax.imshow(d.T, origin="lower", cmap=CMAP[f], vmin=vmin, vmax=vmax,
                                   extent=(0, x.max() + x[1], 0, y.max() + y[1]), interpolation="nearest")
                ax.set_xticks([])
                ax.set_yticks([])
                if i == 0:
                    ax.set_title(f"{LABEL[f]}, {HEIGHT_LABEL.get(hgt, f'$z={z0:.3g}$')}", fontsize=7)
                if hgt != "mid":
                    ax.text(0.03, 0.03, rf"$z={z0:.3g}$", transform=ax.transAxes, fontsize=6,
                            color="k", bbox=dict(fc="white", ec="none", alpha=0.7, pad=1))
            axes[i, 0].set_ylabel(info.label(), fontsize=7)
    for f in reversed(fields_):
        cb = fig.colorbar(ims[f], ax=axes, location="right", fraction=0.025, pad=0.01, aspect=25 * nr)
        cb.set_label(LABEL[f], fontsize=7)
        cb.ax.tick_params(labelsize=6)
    fig.suptitle(", ".join(f"$t={t:.2f}$" for t in times), fontsize=7)
    if out_png:
        fig.savefig(out_png, dpi=200)
        plt.close(fig)
    return fig


def _frame_job(args):
    """Worker for the movie: one frame."""
    frame(*args)
    return args[5]


def get_parser() -> argparse.ArgumentParser:
    """Command-line options (see the module docstring).

    Returns
    -------
    argparse.ArgumentParser
    """
    p = argparse.ArgumentParser(prog="python -m rbc.vis.slices", description=__doc__.split("\n")[0])
    p.add_argument("runs", type=Path, nargs="+", help="run directories (one row each)")
    p.add_argument("--fields", nargs="+", default=["T", "w"], choices=list(CMAP), help="fields to show (columns)")
    p.add_argument("--heights", nargs="+", default=["mid", "bl"], help="mid, bl (1/(4 Nu)), bl-edge (1/(2 Nu)) or numbers")
    p.add_argument("--write", type=int, default=-1, help="snapshot index (per run; -1 = last)")
    p.add_argument("--lim", type=float, default=3.0, help="colour limit in rms units")
    p.add_argument("--keep-negative", action="store_true", help="keep points where the squared vorticity is negative (aliasing) instead of clipping them")
    p.add_argument("--movie", action="store_true", help="render every snapshot and assemble an MP4 (needs ffmpeg)")
    p.add_argument("--jobs", type=int, default=4, help="parallel processes")
    p.add_argument("--fps", type=int, default=30, help="movie frame rate")
    p.add_argument("--out", type=Path, default=Path("figures"), help="output folder")
    p.add_argument("--name", default="slices", help="output file name")
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
    style.use("paper", usetex=False if (a.no_tex or a.movie) else None)
    a.out.mkdir(parents=True, exist_ok=True)

    if not a.movie:
        fig = frame(a.runs, a.fields, a.heights, a.write, a.lim, keep_negative=a.keep_negative)
        print(style.save(fig, a.out / a.name, png=True))
        return
    n = min(len(fl.writes(Path(r).resolve())) for r in a.runs)
    fdir = a.out / f"{a.name}_frames"
    fdir.mkdir(exist_ok=True)
    jobs = [(a.runs, a.fields, a.heights, k, a.lim, fdir / f"frame_{k:05d}.png", 1.0, a.keep_negative)
            for k in range(n)]
    with ProcessPoolExecutor(a.jobs) as ex:
        for done in ex.map(_frame_job, jobs):
            pass
    if shutil.which("ffmpeg"):
        subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-r", str(a.fps), "-i", str(fdir / "frame_%05d.png"),
                        "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2", "-pix_fmt", "yuv420p",
                        str(a.out / f"{a.name}.mp4")], check=True)
        print(a.out / f"{a.name}.mp4")


if __name__ == "__main__":
    main()
