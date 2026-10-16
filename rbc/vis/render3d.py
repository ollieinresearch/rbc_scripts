"""3D views of one run: temperature and vorticity, as isosurfaces or volume renders (V6),
optionally with the in-simulation spectra at the same time underneath (movie frames).

``iso``
    theta/theta_rms(z) = +-ISO_T (warm red and cool blue sheets: plumes), and the vorticity
    magnitude at ISO_OMEGA times omega_rms, coloured by w/w_rms (vortex tubes, and whether
    they sit in rising or sinking fluid). Thresholds are in rms units, so the same numbers
    mean the same thing at every Pr and in every frame.
``volume``
    theta/theta_rms transparent within +-1 and opaque by +-3; log10 of the vorticity
    magnitude over omega_rms transparent below 0 and opaque by 1.5. Values past the end of
    a transfer function are drawn at its end colour, so a fixed range never rescales.

Usage::

    python -m rbc.vis.render3d RUN                       # last snapshot, iso
    python -m rbc.vis.render3d RUN --mode volume --spectra --movie --jobs 16
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

from rbc import style
from rbc.vis import fields as fl

#: Camera position, focal point and up vector (view direction only; the box is refitted).
CAMERA = [(4.0, 4.0, 0.5), (1.0, 1.0, 0.35), (0.0, 0.0, 1.0)]


def _pv():
    """Import pyvista, starting a virtual display on headless nodes if possible."""
    import pyvista as pv
    if not os.environ.get("DISPLAY") and shutil.which("Xvfb"):
        try:
            pv.start_xvfb()
        except Exception:
            pass
    pv.global_theme.allow_empty_mesh = True
    return pv


def _uniform_z(a, z, nzu):
    """Linear resample of the last axis from z (ascending) to nzu uniform points."""
    zu = np.linspace(z[0], z[-1], nzu)
    j = np.clip(np.searchsorted(z, zu) - 1, 0, z.size - 2)
    w = ((zu - z[j]) / (z[j + 1] - z[j])).astype(a.dtype)
    return a[..., j] * (1 - w) + a[..., j + 1] * w, zu


def _set_opacity(vol, values, alpha):
    """Opacity transfer function in data units, set on VTK directly (pyvista's
    array convention has changed between versions: 0-1 vs 0-255)."""
    from vtkmodules.vtkCommonDataModel import vtkPiecewiseFunction
    f = vtkPiecewiseFunction()
    for v, a in zip(values, alpha):
        f.AddPoint(float(v), float(a))
    vol.GetProperty().SetScalarOpacity(f)


def render(run, k, o) -> tuple[np.ndarray, float]:
    """Render one snapshot (temperature left, vorticity right) off screen.

    Parameters
    ----------
    run : pathlib.Path
        Run directory.
    k : int
        Snapshot index.
    o : argparse.Namespace
        Options (mode, thresholds, stride, window, zoom, keep_negative, nz_uniform).

    Returns
    -------
    image : numpy.ndarray
        RGB screenshot.
    time : float
        Simulation time of the snapshot.
    """
    pv = _pv()
    info = fl.RunInfo.load(run)
    wr = fl.writes(info.path)[k]
    s = o.stride
    with h5py.File(wr.file, "r") as h:
        x, y, z = fl.grids(h)
        x, y = x[::s], y[::s]
        T = fl.volume(h, "T", wr.index)[::s, ::s]
        W = fl.volume(h, "w", wr.index)[::s, ::s]
        V = fl.volume(h, "vort", wr.index)[::s, ::s]
    th = fl.normalized(info, "T", T, z).astype(np.float32)
    wn = fl.normalized(info, "w", W, z).astype(np.float32)
    om = fl.normalized(info, "vort", V, z, keep_negative=True).astype(np.float32)
    neg = om < 0
    om = np.maximum(om, 0)

    p = pv.Plotter(shape=(1, 2), off_screen=True, border=False, window_size=o.window)
    p.set_background("white")
    if o.mode == "iso":
        g = pv.RectilinearGrid(x, y, z)
        g.point_data["theta"] = th.ravel(order="F")
        g.point_data["omega"] = om.ravel(order="F")
        g.point_data["w"] = wn.ravel(order="F")
        p.subplot(0, 0)
        p.add_mesh(g.outline(), color="gray", line_width=1)
        surf = g.contour([-o.iso_T, o.iso_T], scalars="theta")
        p.add_mesh(surf, scalars="theta", cmap="RdBu_r", clim=[-3, 3], show_scalar_bar=False, smooth_shading=True)
        p.subplot(0, 1)
        p.add_mesh(g.outline(), color="gray", line_width=1)
        tubes = g.contour([o.iso_omega], scalars="omega")
        if tubes.n_points:
            tubes = tubes.sample(g) if "w" not in tubes.point_data else tubes
        p.add_mesh(tubes, scalars="w", cmap="PuOr_r", clim=[-3, 3], show_scalar_bar=False, smooth_shading=True)
    else:
        nzu = o.nz_uniform or 2 * z.size
        thu, zu = _uniform_z(th, z, nzu)
        omu, _ = _uniform_z(np.log10(np.maximum(om, 1e-3)), z, nzu)
        sp = (x[1] - x[0], y[1] - y[0], zu[1] - zu[0])
        lin = np.linspace(-3, 3, 256)
        lin_o = np.linspace(-1, 1.5, 256)
        for col, (data, cmap, clim, op) in enumerate((
                (thu, "RdBu_r", [-3, 3], np.clip((np.abs(lin) - 1) / 2, 0, 1) * 0.8),
                (omu, "magma", [-1, 1.5], np.clip(lin_o / 1.5, 0, 1) * 0.8))):
            img = pv.ImageData(dimensions=data.shape, spacing=sp, origin=(x[0], y[0], zu[0]))
            img.point_data["f"] = data.ravel(order="F")
            p.subplot(0, col)
            p.add_mesh(img.outline(), color="gray", line_width=1)
            vol = p.add_volume(img, scalars="f", cmap=cmap, clim=clim, opacity="linear", show_scalar_bar=False)
            _set_opacity(vol, np.linspace(*clim, op.size), op)
    if o.keep_negative and neg.any():
        idx = np.argwhere(neg)
        pts = np.c_[x[idx[:, 0]], y[idx[:, 1]], z[idx[:, 2]]]
        p.subplot(0, 1)
        p.add_points(pts, color="magenta", point_size=3, render_points_as_spheres=True)
    for c in (0, 1):
        p.subplot(0, c)
        p.camera_position = CAMERA      # sets the view direction only;
        p.reset_camera()                # fit the box and reset the clipping range
        p.camera.zoom(o.zoom)
    img = p.screenshot(return_img=True)
    p.close()
    return img, wr.time


def spectra_at(run, t):
    """In-simulation spectra nearest a time, with their time means.

    Parameters
    ----------
    run : pathlib.Path
        Run directory.
    t : float
        Simulation time.

    Returns
    -------
    tuple or None
        ``(spectra, time)``: ``spectra[(field, panel)] = (n, E(t), mean E)`` for the
        ``nz`` and ``h`` panels; None if the run has no spectra.
    """
    from rbc.spectra import core as cs
    files = cs.spectra_files(run)
    if not files:
        return None
    time, modes, data, attrs = cs.load_spectra(files)
    i = int(np.argmin(np.abs(time - t)))
    out = {}
    for fld in cs.FIELDS:
        series = cs.panel_series(modes, data[fld], attrs)
        for key in ("nz", "h"):
            n, E = series[key]
            out[(fld, key)] = (n, E[i], cs.time_average(time, E))
    return out, time[i]


def compose(run, k, o, path):
    """Render, add colour bars (and spectra), and save one frame.

    Parameters
    ----------
    run : pathlib.Path
        Run directory.
    k : int
        Snapshot index.
    o : argparse.Namespace
        Options.
    path : pathlib.Path
        PNG to write.

    Returns
    -------
    pathlib.Path
    """
    info = fl.RunInfo.load(run)
    img, t = render(run, k, o)
    rows = 2 if o.spectra else 1
    w_in = 7.0
    fig = plt.figure(figsize=(w_in, w_in * img.shape[0] / img.shape[1] + (2.0 if o.spectra else 0) + 0.5),
                     layout="constrained")
    gs = fig.add_gridspec(rows, 4, height_ratios=[img.shape[0] / img.shape[1] * 4, 1.3][:rows])
    ax = fig.add_subplot(gs[0, :])
    ax.imshow(img)
    ax.axis("off")
    ax.set_title(f"{info.label()},  $t = {t:.3f}$", fontsize=9)
    if o.mode == "iso":
        bars = [(r"$\theta/\theta_{\rm rms}$ (surfaces at $\pm" + f"{o.iso_T:g}" + "$)", "RdBu_r", (-3, 3)),
                (r"$w/w_{\rm rms}$ on $|\omega|/\omega_{\rm rms}=" + f"{o.iso_omega:g}" + "$", "PuOr_r", (-3, 3))]
    else:
        bars = [(r"$\theta/\theta_{\rm rms}$", "RdBu_r", (-3, 3)),
                (r"$\log_{10}|\omega|/\omega_{\rm rms}$", "magma", (-1, 1.5))]
    for side, (lab, cmap, lim) in zip(("left", "right"), bars):
        cax = ax.inset_axes([0.08 if side == "left" else 0.58, 0.07, 0.34, 0.025])
        cb = fig.colorbar(ScalarMappable(Normalize(*lim), cmap), cax=cax, orientation="horizontal")
        cb.set_label(lab, fontsize=7)
        cb.ax.tick_params(labelsize=6)
    if o.spectra:
        sp = spectra_at(info.path, t)
        if sp is not None:
            spec, ts = sp
            for j, (fld, key) in enumerate([("vel", "nz"), ("vel", "h"), ("temp", "nz"), ("temp", "h")]):
                a = fig.add_subplot(gs[1, j])
                n, E, Em = spec[(fld, key)]
                m = (n >= 1) & (E > 0)
                a.loglog(n[m], Em[m], color="0.7", lw=1, label="time mean")
                a.loglog(n[m], E[m], color="C0" if fld == "vel" else "C3", lw=0.9, label=f"$t={ts:.2f}$")
                a.set_xlabel(r"$n_z$" if key == "nz" else r"$k_h$", fontsize=7)
                a.set_title(("velocity" if fld == "vel" else "temperature"), fontsize=7)
                a.tick_params(labelsize=6)
                if j == 0:
                    a.legend(fontsize=6)
    fig.savefig(path, dpi=o.dpi)
    plt.close(fig)
    return path


def _job(args):
    """Worker for the movie: one frame."""
    return compose(*args)


def get_parser() -> argparse.ArgumentParser:
    """Command-line options (see the module docstring).

    Returns
    -------
    argparse.ArgumentParser
    """
    p = argparse.ArgumentParser(prog="python -m rbc.vis.render3d", description=__doc__.split("\n")[0])
    p.add_argument("run", type=Path, help="run directory")
    p.add_argument("--mode", choices=("iso", "volume"), default="iso", help="isosurfaces or volume rendering")
    p.add_argument("--iso-T", type=float, default=1.5, help="temperature isosurfaces at +- this many theta_rms")
    p.add_argument("--iso-omega", type=float, default=3.0, help="vorticity isosurface level in omega_rms")
    p.add_argument("--keep-negative", action="store_true", help="mark points where the squared vorticity is negative (aliasing)")
    p.add_argument("--spectra", action="store_true", help="add the spectra at the same time underneath")
    p.add_argument("--write", type=int, default=-1, help="snapshot index (-1 = last)")
    p.add_argument("--movie", action="store_true", help="render every snapshot and assemble an MP4 (needs ffmpeg)")
    p.add_argument("--jobs", type=int, default=4, help="parallel processes")
    p.add_argument("--stride", type=int, default=1, help="horizontal subsampling (memory/speed)")
    p.add_argument("--nz-uniform", type=int, default=0, help="volume mode z resample (default 2 nz)")
    p.add_argument("--window", type=int, nargs=2, default=[1600, 700], help="render size in pixels (width height)")
    p.add_argument("--zoom", type=float, default=1.0, help="camera zoom")
    p.add_argument("--dpi", type=int, default=200, help="output resolution")
    p.add_argument("--fps", type=int, default=30, help="movie frame rate")
    p.add_argument("--out", type=Path, help="output folder")
    return p


def main(argv=None) -> None:
    """Command-line entry point.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments (default: ``sys.argv[1:]``).
    """
    a = get_parser().parse_args(argv)
    style.use("diag")
    run = a.run.resolve()
    out = a.out or run / "visualization"
    out.mkdir(parents=True, exist_ok=True)
    tag = f"{a.mode}{'_spectra' if a.spectra else ''}"
    if not a.movie:
        print(compose(run, a.write, a, out / f"render_{tag}.png"))
        return
    fdir = out / f"frames_{tag}"
    fdir.mkdir(exist_ok=True)
    n = len(fl.writes(run))
    with ProcessPoolExecutor(a.jobs) as ex:
        list(ex.map(_job, [(run, k, a, fdir / f"frame_{k:05d}.png") for k in range(n)]))
    if shutil.which("ffmpeg"):
        subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-r", str(a.fps), "-i", str(fdir / "frame_%05d.png"),
                        "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2", "-pix_fmt", "yuv420p",
                        str(out / f"movie_{tag}.mp4")], check=True)
        print(out / f"movie_{tag}.mp4")


if __name__ == "__main__":
    main()
