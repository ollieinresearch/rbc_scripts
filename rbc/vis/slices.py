"""Horizontal slices of any number of runs side by side (V1, V2).

Rows are runs, in the order given; columns are (field, height) pairs, so midplane and
boundary-layer slices can share one figure. Colour limits are fixed in normalized units
(see :mod:`rbc.vis.fields`), so rows at different Pr are directly comparable and a movie
never rescales. Usage::

    python -m rbc.vis.slices RUN_A RUN_B RUN_C --fields T w --heights mid bl
    python -m rbc.vis.slices RUN_A RUN_B RUN_C RUN_D --movie
    python -m rbc.vis.slices RUN_A RUN_B --write -1 12          # a chosen snapshot per run
    python -m rbc.vis.slices RUN_A RUN_B --labels "Pr = 1" "Pr = 10^{-3.5}"

Heights: mid (z = 1/2), bl (z = delta_T/2 = 1/(4 Nu)), bl-edge (z = delta_T), or a number.

Each row is labelled with its run (Ra, Pr, and N_z when two rows would otherwise read
the same, or ``--labels``) and the time of its snapshot. Runs rarely have the same
number of snapshots, so a movie pairs them by position: ``--align end`` (default) uses
the last n snapshots of every run, n being the shortest run's count, and ``--align start``
the first n.

Each plane is Fourier-interpolated to at least ``--points`` points per side (exact for these
periodic fields: the resolved field drawn smoothly rather than as blocks), and drawn with
antialiasing, so a plane with more points than the panel has pixels is averaged rather
than subsampled. Movie frames are written at ``--dpi``. They render in parallel, as many
at once as the cores and memory allow (:mod:`rbc.vis.parallel`, ``--jobs 0``); frames
already written are kept, so an interrupted movie resumes (``--redo`` starts over).
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np

from rbc import style
from rbc.vis import fields as fl
from rbc.vis import parallel

CMAP = {"T": "RdBu_r", "w": "PuOr_r", "u": "PuOr_r", "v": "PuOr_r", "vort": "magma"}
LABEL = {"T": r"$\theta/\theta_{\rm rms}$", "w": r"$w/w_{\rm rms}$", "u": r"$u/u_{\rm rms}$",
         "v": r"$v/v_{\rm rms}$", "vort": r"$\log_{10}|\omega|/\omega_{\rm rms}$"}
HEIGHT_LABEL = {"mid": r"$z=1/2$", "bl": r"$z=\delta_T/2$", "bl-edge": r"$z=\delta_T$"}
#: Fixed memory of one frame worker (Python, NumPy, Matplotlib), bytes; measured ~0.25 GB.
BASE_MEMORY = 0.4 * parallel.GB


def row_labels(infos, labels=None) -> list[str]:
    """Row labels: ``labels`` if given, else Ra and Pr, with N_z added where needed.

    Parameters
    ----------
    infos : list of rbc.vis.fields.RunInfo
        The runs.
    labels : list of str, optional
        One label per run (TeX allowed); overrides the default.

    Returns
    -------
    list of str

    Raises
    ------
    SystemExit
        If the number of labels does not match the number of runs.
    """
    if labels:
        if len(labels) != len(infos):
            raise SystemExit(f"--labels: {len(labels)} labels for {len(infos)} runs")
        return list(labels)
    base = [i.label() for i in infos]
    out = list(base)
    for k, info in enumerate(infos):
        if base.count(base[k]) > 1:                  # same (Ra, Pr): tell them apart by N_z
            nz = (fl.parse_run_path(info.path) or {}).get("nz")
            out[k] = base[k] + "\n" + (rf"$N_z={nz}$" if nz else info.path.name)
    return out


def frame(rows, fields_, heights, lim, out_png=None, width=1.0, keep_negative=False, points=512, dpi=300):
    """One figure: rows are runs, columns are (field, height) pairs.

    Parameters
    ----------
    rows : list of tuple
        ``(info, write, label)`` per row: the run (:class:`rbc.vis.fields.RunInfo`), the
        snapshot write to show (:class:`rbc.vis.fields.Write`) and the row label.
    fields_ : list of str
        Fields (keys of ``CMAP``).
    heights : list of str or float
        Heights (see :meth:`rbc.vis.fields.RunInfo.height`).
    lim : float
        Colour limit in rms units (vorticity uses a fixed log range).
    out_png : pathlib.Path, optional
        Save here (under a temporary name, renamed when complete) and close the figure.
    width : float, optional
        Fraction of the text width.
    keep_negative : bool, optional
        Keep negative ``|ω|²`` points.
    points : int, optional
        Fourier-interpolate each plane to at least this many points per side (0: off).
    dpi : int, optional
        Resolution of ``out_png``.

    Returns
    -------
    matplotlib.figure.Figure
    """
    cols = [(f, h) for h in heights for f in fields_]
    nr, nc = len(rows), len(cols)
    w_in, _ = style.figsize(width)
    fig, axes = plt.subplots(nr, nc, figsize=(w_in, w_in / (nc + 0.6) * nr * 1.08 + 0.25), squeeze=False,
                             layout="constrained")
    ims = {}
    for i, (info, wr, label) in enumerate(rows):
        with h5py.File(wr.file, "r") as h:
            x, y, _ = fl.grids(h)
            for j, (f, hgt) in enumerate(cols):
                z0 = info.height(hgt)
                raw = fl.fourier_resample(fl.plane(h, f, wr.index, z0), points)
                d = fl.normalized(info, f, raw, z0, keep_negative)
                if f == "vort":
                    d = np.log10(np.clip(d, 1e-3, None))
                    vmin, vmax = -1.0, 1.5
                else:
                    vmin, vmax = -lim, lim
                ax = axes[i, j]
                ims[f] = ax.imshow(d.T, origin="lower", cmap=CMAP[f], vmin=vmin, vmax=vmax,
                                   extent=(0, x.max() + x[1], 0, y.max() + y[1]), interpolation="antialiased")
                ax.set_xticks([])
                ax.set_yticks([])
                if i == 0:
                    ax.set_title(f"{LABEL[f]}, {HEIGHT_LABEL.get(hgt, f'$z={z0:.3g}$')}", fontsize=7)
                if hgt != "mid":
                    ax.text(0.03, 0.03, rf"$z={z0:.3g}$", transform=ax.transAxes, fontsize=6,
                            color="k", bbox=dict(fc="white", ec="none", alpha=0.7, pad=1))
        axes[i, 0].set_ylabel(f"{label}\n$t={wr.time:.2f}$", fontsize=7)
    for f in reversed(fields_):
        cb = fig.colorbar(ims[f], ax=axes, location="right", fraction=0.025, pad=0.01, aspect=25 * nr)
        cb.set_label(LABEL[f], fontsize=7)
        cb.ax.tick_params(labelsize=6)
    if out_png:
        out_png = Path(out_png)
        tmp = out_png.with_name(out_png.stem + ".part.png")
        fig.savefig(tmp, dpi=dpi)
        plt.close(fig)
        os.replace(tmp, out_png)
    return fig


def _frame_job(args):
    """Worker for the movie: one frame (``args`` as for :func:`frame`)."""
    frame(*args)


def select(writes, align, n=None) -> list[list]:
    """The writes of each run to pair up, frame by frame.

    Parameters
    ----------
    writes : list of list of rbc.vis.fields.Write
        Every run's writes, in time order.
    align : {"end", "start"}
        Pair the last (or first) n writes of every run.
    n : int, optional
        Number of frames (default: the shortest run's count).

    Returns
    -------
    list of list
        ``frames[k][i]``: the write of run ``i`` in frame ``k``.
    """
    n = min(len(w) for w in writes) if n is None else n
    pick = (lambda w: w[len(w) - n:]) if align == "end" else (lambda w: w[:n])
    return [list(r) for r in zip(*[pick(w) for w in writes])]


def frame_memory(writes, fields_, heights, points) -> float:
    """Estimated peak memory of one frame, in bytes.

    Per plane: the 8 levels read for the interpolation in z (float64) and the
    Fourier-interpolated image with Matplotlib's working copies.

    Parameters
    ----------
    writes : list of rbc.vis.fields.Write
        One write per run (for the grid sizes).
    fields_, heights : list
        Columns of the figure.
    points : int
        Minimum points per side after interpolation.

    Returns
    -------
    float
    """
    total = BASE_MEMORY
    for wr in writes:
        nx, ny, _ = fl.grid_shape(wr)
        m = max(nx, ny, points)
        total += len(fields_) * len(heights) * (8 * 8 * nx * ny + 8 * 8 * m * m)
    return total


def get_parser() -> argparse.ArgumentParser:
    """Command-line options (see the module docstring).

    Returns
    -------
    argparse.ArgumentParser
    """
    p = argparse.ArgumentParser(prog="python -m rbc.vis.slices", description=__doc__.split("\n")[0])
    p.add_argument("runs", type=Path, nargs="+", help="run directories, any number (one row each, in this order)")
    p.add_argument("--fields", nargs="+", default=["T", "w"], choices=list(CMAP), help="fields to show (columns)")
    p.add_argument("--heights", nargs="+", default=["mid", "bl"], help="mid, bl (1/(4 Nu)), bl-edge (1/(2 Nu)) or numbers")
    p.add_argument("--write", type=int, nargs="+", default=[-1],
                   help="snapshot index: one for every run, or one per run (-1 = last)")
    p.add_argument("--labels", nargs="+", help="row labels, one per run (default: Ra, Pr and, if needed, N_z)")
    p.add_argument("--align", choices=("end", "start"), default="end",
                   help="movie: pair the last (end) or first (start) snapshots of the runs")
    p.add_argument("--lim", type=float, default=3.0, help="colour limit in rms units")
    p.add_argument("--keep-negative", action="store_true", help="keep points where the squared vorticity is negative (aliasing) instead of clipping them")
    p.add_argument("--movie", action="store_true", help="render every snapshot and assemble an MP4 (needs ffmpeg)")
    p.add_argument("--jobs", type=int, default=0, help="parallel processes (0 = as many as cores and memory allow)")
    p.add_argument("--redo", action="store_true", help="movie: re-render frames that already exist")
    p.add_argument("--fps", type=int, default=30, help="movie frame rate")
    p.add_argument("--points", type=int, default=512,
                   help="Fourier-interpolate each plane to at least this many points per side (0: off)")
    p.add_argument("--dpi", type=int, default=300, help="resolution of movie frames")
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
    infos = [fl.RunInfo.load(r) for r in a.runs]
    labels = row_labels(infos, a.labels)
    writes = [fl.writes(i.path) for i in infos]          # listed once, not in every frame
    for info, w in zip(infos, writes):
        if not w:
            raise SystemExit(f"{info.path}: no snapshots")

    if not a.movie:
        if len(a.write) not in (1, len(infos)):
            raise SystemExit(f"--write: give one index, or one per run ({len(infos)})")
        ks = a.write * len(infos) if len(a.write) == 1 else a.write
        rows = [(info, w[k], lab) for info, w, k, lab in zip(infos, writes, ks, labels)]
        fig = frame(rows, a.fields, a.heights, a.lim, keep_negative=a.keep_negative, points=a.points)
        print(style.save(fig, a.out / a.name, png=True))
        return

    frames = select(writes, a.align)
    counts = ", ".join(str(len(w)) for w in writes)
    print(f"{len(frames)} frames ({a.align}-aligned; snapshots per run: {counts})")
    fdir = a.out / f"{a.name}_frames"
    settings = dict(runs=[str(i.path) for i in infos], labels=labels, fields=a.fields, heights=a.heights,
                    align=a.align, lim=a.lim, keep_negative=a.keep_negative, points=a.points, dpi=a.dpi)
    parallel.prepare_frames(fdir, settings, len(frames), a.redo)
    items = [(fdir / f"frame_{k:05d}.png",
              ([(info, wr, lab) for info, wr, lab in zip(infos, ws, labels)], a.fields, a.heights, a.lim,
               fdir / f"frame_{k:05d}.png", 1.0, a.keep_negative, a.points, a.dpi))
             for k, ws in enumerate(frames)]
    jobs = parallel.auto_jobs(len(items), frame_memory(frames[0], a.fields, a.heights, a.points), a.jobs)
    lost = parallel.render_all(_frame_job, items, jobs)
    if lost:
        raise SystemExit(f"{len(lost)} of {len(items)} frames failed (first: {lost[0].name}); no movie assembled. "
                         "Rerun the same command to render only those.")
    movie = parallel.make_movie(fdir, "frame_%05d.png", a.out / f"{a.name}.mp4", a.fps)
    if movie:
        print(movie)


if __name__ == "__main__":
    main()
