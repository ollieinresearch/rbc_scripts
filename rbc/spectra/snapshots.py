"""Power spectra from saved grid snapshots, and per-write plots of any spectra file.

New runs compute spectra inside the simulation (``spectra/spectra_s<N>.h5``), so this
is mainly for

1. runs made before that existed: it converts ``snapshots_s<N>.h5`` into
   ``snapshots_s<N>_spectra.h5`` with exactly the same definitions, and
2. per-write PNGs (from either kind of spectra file), in ``res_check_3d/`` and
   ``res_check_temp/``.

The snapshot grid values at scales = 1 define the Dedalus coefficients exactly, so the
same forward transforms Dedalus uses are applied (real FFT in x and y, DCT-II in z; see
:mod:`rbc.spectra.core`): no interpolation, no artificial periodic extension in z, no
Nyquist artefact, and the result equals, to round-off, what the simulation would have
written. (An earlier script assumed Gauss–Lobatto nodes while Dedalus uses Gauss
nodes, so every value was interpolated from the wrong height.)

Files are split across MPI ranks with ``[rank::size]``, so::

    srun -n 32 python -m rbc.spectra.snapshots snapshots/snapshots_s*.h5 --plot-every 0

processes 32 files at once. Missing y-limits are chosen automatically.
"""

import argparse
import re
from pathlib import Path

import h5py as h5
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from rbc import style
from rbc.spectra import core as cs

#: Domain size; must match the simulation.
LX, LY, LZ = 2.0, 2.0, 1.0

style.use("diag")


# ===========================================================================
#  Snapshot file -> spectra file
# ===========================================================================

def _load_velocity(f, ti):
    """Velocity components of write ``ti`` (saved a few different ways over the project).

    Returns
    -------
    tuple
        List of the three component arrays, and the dataset they came from.
    """
    tasks = f["tasks"]
    if "u" in tasks and tasks["u"].ndim == 5:          # (t, 3, nx, ny, nz)
        d = tasks["u"]
        return [np.asarray(d[ti, i]) for i in range(3)], d
    if "x" in tasks:
        return [np.asarray(tasks[k][ti]) for k in ("x", "y", "w")], tasks["x"]
    return [np.asarray(tasks[k][ti]) for k in ("u", "v", "w")], tasks["u"]


def _task_scales(dset):
    """Grid scales a snapshot task was written at (1 if not recorded)."""
    s = dset.attrs.get("scales", None)
    if s is None:
        return (1, 1, 1)
    s = tuple(float(v) for v in np.atleast_1d(s))
    return s if len(s) == 3 else (s[0],) * 3


def _check_z_grid(dset, nz):
    """Check the file's z scale (if attached) against the Gauss–Chebyshev grid.

    Raises
    ------
    ValueError
        If they differ (wrong ``LZ`` or bounds).
    """
    try:
        z_file = np.asarray(dset.dims[dset.ndim - 1][0])
    except Exception:
        return
    if z_file.size != nz:
        return
    z_gauss = cs.chebyshev_grid(nz, (0.0, LZ))
    if not np.allclose(z_file, z_gauss, atol=1e-10 * LZ):
        raise ValueError("Snapshot z grid is not the Gauss-Chebyshev grid on "
                         f"(0, {LZ}); check LZ / bounds before trusting spectra.")


def snapshot_to_spectra(fp, out, overwrite=False):
    """Compute spectra for every write in one snapshot file.

    The output is written to ``<out>.part`` and renamed when complete.

    Parameters
    ----------
    fp : pathlib.Path
        ``snapshots_s<N>.h5``.
    out : pathlib.Path
        Spectra file to write.
    overwrite : bool, optional
        Recompute even if ``out`` exists.
    """
    if out.exists() and not overwrite:
        print(f"Skipping {fp.name}: {out.name} already exists.")
        return
    tmp = out.with_suffix(".h5.part")                  # never leave a half file behind
    with h5.File(fp, "r") as f:
        times = f["scales/sim_time"][:]
        iters = f["scales/iteration"][:] if "iteration" in f["scales"] else np.zeros(len(times), int)
        T_d = f["tasks"]["temperature"]
        scales = _task_scales(T_d)
        nz_grid = T_d.shape[-1]
        _check_z_grid(T_d, nz_grid)
        z = cs.chebyshev_grid(nz_grid, (0.0, LZ))
        b0 = (LZ - z)[None, None, :]                   # conductive profile, as in rbc_script

        red = None
        for ti, t in enumerate(times):
            if np.isclose(t, 0.0):
                continue                               # initial condition: nothing to see
            comps, _ = _load_velocity(f, ti)
            uc = [cs.grid_to_coeff(c, scales) for c in comps]
            gms_vel = float(sum(np.mean(c**2) for c in comps))
            b = np.asarray(T_d[ti]) - b0               # b = T - b0, same as in-sim
            bc = cs.grid_to_coeff(b, scales)
            gms_temp = float(np.mean(b**2))
            del comps, b

            if red is None:
                Nx, Ny, Nz = bc.shape
                red = cs.SpectralReducer(np.arange(Nx), np.arange(Ny), np.arange(Nz), Nx, Ny, Nz)
                cs.create_spectra_file(tmp, red, dict(source="snapshots", snapshot_file=fp.name,
                                                      Nx=Nx, Ny=Ny, Nz=Nz, Lx=LX, Ly=LY, Lz=LZ))
            vel = red.reduce(red.power(uc))
            temp = red.reduce(red.power([bc]))
            del uc, bc
            cs.append_spectra(tmp, t, int(iters[ti]), vel, temp, gms_vel, gms_temp)
    if tmp.exists():
        tmp.rename(out)
        print(f"Wrote {out}")


# ===========================================================================
#  Per-write plots
# ===========================================================================

def plot_write(t, modes, spec, attrs, savepath, title, ylims=None, gms=np.nan):
    """Four-panel figure (3D shells, z, horizontal, horizontal shells) of one write.

    Parameters
    ----------
    t : float
        Simulation time.
    modes, spec, attrs : dict
        Mode tables, one write's spectra and file attributes.
    savepath : pathlib.Path
        PNG to write.
    title : str
        Figure title.
    ylims : list of tuple, optional
        Log10 y limits per panel (default: automatic).
    gms : float, optional
        Grid mean square, for the Parseval check in the title.
    """
    series = cs.panel_series(modes, spec, attrs)
    fig, axes = plt.subplots(2, 2, figsize=(16, 14), layout="constrained")
    for i, key in enumerate(("3d", "nz", "h", "perp")):
        (r, c), ptitle, xlabel = cs.PANEL_META[key]
        ax = axes[r, c]
        n, E = series[key]
        m = (n >= 1) & (E > 0)
        ax.loglog(n[m], E[m], ".-", linewidth=1.2, markersize=4, label="fluctuations")
        extra = [E[m]]
        if key == "nz":
            Em = spec["mean_E_nz"]
            mm = (modes["nz"] >= 1) & (Em > 0)
            ax.loglog(modes["nz"][mm], Em[mm], "--", linewidth=1.2, label="horizontal mean")
            ax.legend(fontsize=14, loc="lower left")
            extra.append(Em[mm])
        d = cs.decades(n, E)
        ax.text(0.97, 0.95, f"D = {d['D']:.2f}", transform=ax.transAxes, ha="right",
                va="top", fontsize=16, bbox=dict(boxstyle="round", fc="white", alpha=0.75))
        ax.set_title(ptitle)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("E")
        lo, hi = ylims[i] if ylims is not None else cs.auto_ylim(*extra)
        cs.setup_yaxis(ax, lo, hi)
    total = spec["E_kx"].sum() + spec["mean_E_nz"].sum()
    check = "" if not np.isfinite(gms) else \
        f"\nParseval: sum E / grid mean square - 1 = {total / gms - 1:.1e}"
    fig.suptitle(f"{title}  |  t = {t:.4f}{check}")
    fig.savefig(savepath, dpi=150)
    plt.close(fig)


def plot_file(spec_fp, basepath, every, vel_lims, temp_lims):
    """Per-write figures of one spectra file.

    Parameters
    ----------
    spec_fp : pathlib.Path
        Spectra file.
    basepath : pathlib.Path
        Run directory (figures go to ``res_check_3d/`` and ``res_check_temp/``).
    every : int
        Plot every n-th write (0: none).
    vel_lims, temp_lims : list of tuple or None
        Log10 y limits per panel.
    """
    if every <= 0:
        return
    with h5.File(spec_fp, "r") as f:
        modes = {k: f["modes"][k][:] for k in f["modes"]}
        attrs = dict(f.attrs)
        times = f["scales/sim_time"][:]
        iters = f["scales/iteration"][:]
        gms = {fld: f["scales"][f"grid_mean_square_{fld}"][:] for fld in cs.FIELDS}
        for fld, folder, lims, name in (("vel", "res_check_3d", vel_lims, "Velocity"),
                                        ("temp", "res_check_temp", temp_lims, "Temperature fluctuation")):
            outdir = basepath / folder
            outdir.mkdir(exist_ok=True)
            for ti in range(0, len(times), every):
                spec = {k: f[fld][k][ti] for k in cs.SPEC_KEYS}
                plot_write(times[ti], modes, spec, attrs,
                           outdir / f"iter_{int(iters[ti]):09d}.png",
                           f"{name} spectra", lims, gms[fld][ti])


def _lims(mins, maxs):
    """Parse comma-separated limits into four (min, max) pairs."""
    if mins is None or maxs is None:
        return None
    lo = [float(v) for v in mins.split(",")]
    hi = [float(v) for v in maxs.split(",")]
    return list(zip(lo, hi))


def process(path, every, vel_lims, temp_lims, overwrite):
    """Convert (if a snapshot file) and plot one file.

    Parameters
    ----------
    path : str or pathlib.Path
        ``snapshots_s<N>.h5`` or ``spectra_s<N>.h5``.
    every : int
        Plot every n-th write.
    vel_lims, temp_lims : list of tuple or None
        Log10 y limits.
    overwrite : bool
        Recompute existing snapshot spectra.
    """
    fp = Path(path)
    if re.search(r"snapshots_s\d+\.h5$", fp.name):
        out = fp.parent / f"{fp.stem}_spectra.h5"
        snapshot_to_spectra(fp, out, overwrite)
        basepath = fp.parents[1]
        if out.exists():
            plot_file(out, basepath, every, vel_lims, temp_lims)
    elif re.search(r"spectra_s\d+\.h5$", fp.name):
        plot_file(fp, fp.parents[1], every, vel_lims, temp_lims)
    else:
        print(f"Ignoring {fp} (not a snapshots_sN.h5 or spectra_sN.h5 file)")


def get_parser() -> argparse.ArgumentParser:
    """Command-line options.

    Returns
    -------
    argparse.ArgumentParser
    """
    p = argparse.ArgumentParser(prog="python -m rbc.spectra.snapshots",
                                description="Spectra from snapshots, and per-write spectra plots.")
    p.add_argument("files", nargs="+", type=Path,
                   help="snapshots_s*.h5 (computed, then plotted) or spectra_s*.h5 (plotted)")
    p.add_argument("--plot-every", type=int, default=1, help="plot every n-th write (0 = no plots)")
    p.add_argument("--vmins", help="velocity log10 y minima, 4 comma-separated values (3D, z, horizontal, perp)")
    p.add_argument("--vmaxs", help="velocity log10 y maxima (4 values)")
    p.add_argument("--tmins", help="temperature log10 y minima (4 values)")
    p.add_argument("--tmaxs", help="temperature log10 y maxima (4 values)")
    p.add_argument("--overwrite", action="store_true", help="recompute existing snapshot spectra")
    return p


def main(argv=None) -> None:
    """Process the given files, split across MPI ranks if run under MPI.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments (default: ``sys.argv[1:]``).
    """
    a = get_parser().parse_args(argv)
    vel_lims = _lims(a.vmins, a.vmaxs)
    temp_lims = _lims(a.tmins, a.tmaxs)
    try:
        from mpi4py import MPI
        rank, size = MPI.COMM_WORLD.rank, MPI.COMM_WORLD.size
    except Exception:          # no MPI available: run serially
        rank, size = 0, 1
    all_files = sorted(a.files)
    if rank == 0:
        print(f"{len(all_files)} file(s) across {size} rank(s).")
    for fp in all_files[rank::size]:
        process(fp, a.plot_every, vel_lims, temp_lims, a.overwrite)


if __name__ == "__main__":
    main()
