"""3D views of one run: temperature and vorticity, as isosurfaces or volume renders (V6),
optionally with the in-simulation spectra at the same time underneath (movie frames).

``iso``
    theta/theta_rms(z) = +-ISO_T (warm red and cool blue sheets: plumes), and the vorticity
    magnitude at ISO_OMEGA times omega_rms, coloured by w/w_rms (vortex tubes, and whether
    they sit in rising or sinking fluid). Thresholds are in rms units, so the same numbers
    mean the same thing at every Pr and in every frame.
``volume``
    the transfer functions of the original ``plotting_v3.py``. Temperature T on [0, 1],
    opacity 1 - exp(-1e8 (s - 1/2)^12) with s the position on that range: transparent
    while T is within about 0.2 of 1/2 and opaque beyond about 0.25, so only the hot and
    cold extremes (plumes and boundary layers) show. Vorticity: the saved ``|ω|²`` on an
    asinh scale (linear width 1.5), mapped over its range in that frame, opacity
    1 - exp(-3e11 s^22) with s from -1/2 to 1/2 and zero from 0 up to the middle: only
    the top fifth of the range, the strongest filaments, is opaque. Colour map ``jet``
    for both by default (``--cmap-T``, ``--cmap-omega``). The Chebyshev z points are
    resampled to a uniform grid (``--stretch-z`` places them evenly instead, as the
    original did, which thickens the boundary layers on screen).

Usage::

    python -m rbc.vis.render3d RUN                       # last snapshot, iso
    python -m rbc.vis.render3d RUN --mode volume --spectra --movie
    python -m rbc.vis.render3d RUN --movie --redo        # re-render every frame

Headless nodes: unless ``DISPLAY`` is set, or ``VTK_DEFAULT_OPENGL_WINDOW`` asks for an
EGL or OSMesa window (VTK 9.4 or later), one virtual X server (Xvfb) is started for the
whole job on a free display, before any worker, and stopped at the end; its log is
``xvfb.log`` in the output folder.

Movies: each frame is written to a temporary name and renamed when complete, and frames
already present are kept, so an interrupted movie (or the next post-processing job, with
more snapshots) only renders what is missing. Frames render in parallel, as many at once
as the cores and memory allow (``--jobs 0``, see :mod:`rbc.vis.parallel`): a frame needs
about 0.6 GB plus 100 (volume) to 150 (iso) bytes per grid point, so memory, not cores,
limits large runs; the cores left over go to each worker's software OpenGL threads
(``LP_NUM_THREADS``). The snapshot list and the spectra are read once, in the parent, and
each worker gets only what its frame needs. Workers are replaced every
``--frames-per-process`` frames. If a worker dies (for instance when it loses its X
connection) the pool is rebuilt, Xvfb restarted if it died, and the missing frames tried
again, up to ``--retries`` times; the MP4 is only assembled when every frame exists.
"""
from __future__ import annotations

import argparse
import os
import select
import shutil
import subprocess
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import AsinhNorm, Normalize
from matplotlib.ticker import FuncFormatter, NullLocator

from rbc import style
from rbc.vis import fields as fl
from rbc.vis import parallel

#: Peak memory of one frame: fixed part (Python, VTK, Matplotlib) and bytes per grid point
#: by mode, measured on 72^2 x 36 and 288^2 x 36 grids (0.5 GB + 95 and 58 B/point) and
#: rounded up, since isosurfaces of a turbulent field have many more triangles.
FRAME_BASE = 0.6 * parallel.GB
BYTES_PER_POINT = {"iso": 150, "volume": 100}

#: Camera position, focal point and up vector (view direction only; the box is refitted).
CAMERA = [(4.0, 4.0, 0.5), (1.0, 1.0, 0.35), (0.0, 0.0, 1.0)]


#: Render windows that need no X server (VTK 9.4+ selects them at run time).
HEADLESS_WINDOWS = ("vtkEGLRenderWindow", "vtkOSOpenGLRenderWindow")


class VirtualDisplay:
    """One Xvfb for this process and every worker it starts.

    ``pyvista.start_xvfb`` launches ``Xvfb :99`` from whichever process calls it and
    leaves it running unowned; called from every worker of a pool, that is one race per
    worker for the same display, and a server no process watches. Here the parent starts
    a single server, lets Xvfb choose a free display (``-displayfd``), exports ``DISPLAY``
    before any worker exists, and stops the server at the end.

    Parameters
    ----------
    size : tuple of int
        Screen size in pixels (width, height); renders are off screen, so this only needs
        to be about the window size.
    log : pathlib.Path, optional
        File for the server's own messages.
    """

    def __init__(self, size=(1600, 900), log=None):
        self.size, self.log, self.proc = size, log, None

    @staticmethod
    def needed() -> bool:
        """True if there is no display and VTK was not asked for a headless window."""
        return not os.environ.get("DISPLAY") and \
            os.environ.get("VTK_DEFAULT_OPENGL_WINDOW", "") not in HEADLESS_WINDOWS

    def alive(self) -> bool:
        """True if the server this object started is still running."""
        return self.proc is not None and self.proc.poll() is None

    def start(self, timeout=20.0) -> None:
        """Start Xvfb and set ``DISPLAY``; does nothing if one is not needed.

        Raises
        ------
        RuntimeError
            If Xvfb is not installed or does not report a display in ``timeout`` s.
        """
        if self.alive() or not (self.needed() or self.proc is not None):
            return
        if not shutil.which("Xvfb"):
            raise RuntimeError("no DISPLAY and no Xvfb: run on a node with Xvfb, or set "
                               "VTK_DEFAULT_OPENGL_WINDOW=vtkOSOpenGLRenderWindow with VTK 9.4+")
        r, w = os.pipe()
        err = open(self.log, "a") if self.log else subprocess.DEVNULL
        self.proc = subprocess.Popen(
            ["Xvfb", "-displayfd", str(w), "-screen", "0", f"{self.size[0]}x{self.size[1]}x24",
             "-nolisten", "tcp", "+extension", "GLX"],
            pass_fds=(w,), stdout=err, stderr=err, start_new_session=True)
        os.close(w)
        num = b""
        try:
            while not num.endswith(b"\n"):
                if not select.select([r], [], [], timeout)[0]:
                    raise RuntimeError(f"Xvfb did not start within {timeout:g} s (see {self.log})")
                chunk = os.read(r, 16)
                if not chunk:
                    raise RuntimeError(f"Xvfb exited at start (see {self.log})")
                num += chunk
        finally:
            os.close(r)
        os.environ["DISPLAY"] = f":{int(num)}"
        print(f"Xvfb started on DISPLAY={os.environ['DISPLAY']} (pid {self.proc.pid})")

    def ensure(self) -> None:
        """Restart the server if this object started one and it has died."""
        if self.proc is not None and not self.alive():
            print(f"Xvfb (pid {self.proc.pid}) exited with code {self.proc.returncode}; restarting")
            os.environ.pop("DISPLAY", None)
            self.start()

    def stop(self) -> None:
        """Stop the server, if this object started one."""
        if self.alive():
            self.proc.terminate()
            try:
                self.proc.wait(5)
            except subprocess.TimeoutExpired:
                self.proc.kill()

    def __enter__(self):
        self.start()
        return self

    def __exit__(self, *exc):
        self.stop()


def _pv():
    """Import pyvista (the display is set up by :class:`VirtualDisplay` in the parent)."""
    import pyvista as pv
    pv.global_theme.allow_empty_mesh = True
    return pv


def probe() -> None:
    """Render one tiny scene, so a broken OpenGL set-up fails here with a clear message.

    Raises
    ------
    SystemExit
        If rendering fails.
    """
    try:
        pv = _pv()
        p = pv.Plotter(off_screen=True, window_size=(64, 64))
        p.add_mesh(pv.Sphere())
        p.screenshot(return_img=True)
        p.close()
    except Exception as e:      # pragma: no cover - depends on the node
        raise SystemExit(f"OpenGL rendering does not work here ({e!r}); DISPLAY="
                         f"{os.environ.get('DISPLAY')!r}. See Troubleshooting in the docs.")


def _uniform_z(a, z, nzu):
    """Linear resample of the last axis from z (ascending) to nzu uniform points."""
    zu = np.linspace(z[0], z[-1], nzu)
    j = np.clip(np.searchsorted(z, zu) - 1, 0, z.size - 2)
    w = ((zu - z[j]) / (z[j + 1] - z[j])).astype(a.dtype)
    return a[..., j] * (1 - w) + a[..., j + 1] * w, zu


#: Linear width of the asinh scale applied to |ω|² in volume mode (as in plotting_v3.py).
ASINH_WIDTH = 1.5


def asinh_scale(a, width=ASINH_WIDTH):
    """``width * asinh(a / width)``: linear below ``width``, logarithmic above (as
    :class:`matplotlib.scale.AsinhTransform`)."""
    return width * np.arcsinh(np.asarray(a, dtype=np.float64) / width)


def opacity_T(n=255):
    """Temperature opacity of the original volume renders, sampled across T in [0, 1].

    1 - exp(-1e8 (s - 1/2)^12): 0.5 where T is 0.21 from 1/2, 0.99 at 0.245.

    Parameters
    ----------
    n : int, optional
        Number of samples (pyvista interpolates them over its lookup table).

    Returns
    -------
    numpy.ndarray
    """
    s = np.linspace(0, 1, n)
    return 1.0 - np.exp(-1e8 * (s - 0.5) ** 12)


def opacity_omega(lo, hi, n=255):
    """Vorticity opacity of the original volume renders, sampled across [lo, hi].

    1 - exp(-3e11 s^22) for s from -1/2 to 1/2 (0.5 at s = ±0.30, 0.99 at ±0.32), and
    zero from the value 0 up to the middle, so only the top fifth of the range is opaque
    (and the bottom fifth too, if negative values of ``|ω|²`` are kept).

    Parameters
    ----------
    lo, hi : float
        Range of the (asinh-scaled) data in this frame.
    n : int, optional
        Number of samples.

    Returns
    -------
    numpy.ndarray
    """
    s = np.linspace(-0.5, 0.5, n)
    op = 1.0 - np.exp(-3e11 * s ** 22)
    z0 = int(np.searchsorted(np.linspace(lo, hi, 256), 0))
    if z0 < 125:
        op[max(z0 - 1, 0):125] = 0.0
    return op


def render(info, wr, o) -> tuple[np.ndarray, float, dict]:
    """Render one snapshot (temperature left, vorticity right) off screen.

    Parameters
    ----------
    info : rbc.vis.fields.RunInfo
        The run.
    wr : rbc.vis.fields.Write
        The snapshot write.
    o : argparse.Namespace
        Options (mode, thresholds, stride, window, zoom, keep_negative, nz_uniform).

    Returns
    -------
    image : numpy.ndarray
        RGB screenshot.
    time : float
        Simulation time of the snapshot.
    meta : dict
        For the colour bars: ``omega_range``, the range of ``|ω|²`` in volume mode.
    """
    pv = _pv()
    s = o.stride
    with h5py.File(wr.file, "r") as h:
        x, y, z = fl.grids(h)
        x, y = x[::s], y[::s]
        T = fl.volume(h, "T", wr.index)[::s, ::s]
        W = fl.volume(h, "w", wr.index)[::s, ::s]
        V = fl.volume(h, "vort", wr.index)[::s, ::s]
    meta = {}
    neg = V < 0

    p = pv.Plotter(shape=(1, 2), off_screen=True, border=False, window_size=o.window)
    p.set_background("white")
    if o.mode == "iso":
        th = fl.normalized(info, "T", T, z).astype(np.float32)
        wn = fl.normalized(info, "w", W, z).astype(np.float32)
        om = fl.normalized(info, "vort", V, z).astype(np.float32)
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
        # Transfer functions go through pyvista's own ``opacity=`` (as in plotting_v3.py):
        # pyvista 0.44 rescales the scalars to 0-255 before VTK sees them, so an opacity
        # function set on VTK in data units would not line up with them.
        Vr = V if o.keep_negative else np.maximum(V, 0)
        if o.stretch_z:                     # Chebyshev points placed evenly, as originally
            zu = np.linspace(0.0, 1.0, z.size)
            Tu, Vu = T, Vr
        else:
            nzu = o.nz_uniform or 2 * z.size
            Tu, zu = _uniform_z(T, z, nzu)
            Vu, _ = _uniform_z(Vr, z, nzu)
        A = asinh_scale(np.nan_to_num(Vu)).astype(np.float32)
        lo, hi = float(A.min()), float(A.max())
        hi = hi if hi > lo else lo + 1.0
        meta["omega_range"] = (float(np.sinh(lo / ASINH_WIDTH) * ASINH_WIDTH),
                               float(np.sinh(hi / ASINH_WIDTH) * ASINH_WIDTH))
        sp = (x[1] - x[0], y[1] - y[0], zu[1] - zu[0])
        # pyvista's default opacity unit distance for the original full grid, so the same
        # transfer function gives the same opacity whatever the stride or z resampling
        box = (sp[0] * x.size, sp[1] * y.size, 1.0)
        unit = float(np.linalg.norm(box) / (np.mean([x.size * s, y.size * s, z.size]) - 1))
        for col, (data, cmap, clim, op) in enumerate((
                (np.clip(np.nan_to_num(Tu, nan=0.5), 0.0, 1.0).astype(np.float32), o.cmap_T, [0.0, 1.0],
                 opacity_T()),
                (A, o.cmap_omega, [lo, hi], opacity_omega(lo, hi)))):
            img = pv.ImageData(dimensions=data.shape, spacing=sp, origin=(x[0], y[0], zu[0]))
            img.point_data["f"] = data.ravel(order="F")
            p.subplot(0, col)
            p.add_mesh(img.outline(), color="gray", line_width=1)
            vol = p.add_volume(img, scalars="f", cmap=cmap, clim=clim, opacity=op, opacity_unit_distance=unit,
                               show_scalar_bar=False)
            vol.mapper.interpolate_before_map = False
    if o.mode == "iso" and o.keep_negative and neg.any():
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
    return img, wr.time, meta


class Spectra:
    """The in-simulation spectra of a run, read once, for the panels under each frame.

    Parameters
    ----------
    run : pathlib.Path
        Run directory.
    """

    def __init__(self, run):
        from rbc.spectra import core as cs
        self.time, self.series, self.mean = None, {}, {}
        files = cs.spectra_files(run)
        if not files:
            return
        time, modes, data, attrs = cs.load_spectra(files)
        self.time = time
        for fld in cs.FIELDS:
            series = cs.panel_series(modes, data[fld], attrs)
            for key in ("nz", "h"):
                n, E = series[key]
                self.series[(fld, key)] = (n, E)
                self.mean[(fld, key)] = cs.time_average(time, E)

    def at(self, t):
        """Spectra nearest time ``t``, with the time means.

        Parameters
        ----------
        t : float
            Simulation time.

        Returns
        -------
        tuple or None
            ``(spectra, time)``: ``spectra[(field, panel)] = (n, E(t), mean E)`` for the
            ``nz`` and ``h`` panels; None if the run has no spectra.
        """
        if self.time is None:
            return None
        i = int(np.argmin(np.abs(self.time - t)))
        return {k: (n, E[i], self.mean[k]) for k, (n, E) in self.series.items()}, float(self.time[i])


def compose(info, wr, o, path, spec=None):
    """Render, add colour bars (and spectra), and save one frame.

    Parameters
    ----------
    info : rbc.vis.fields.RunInfo
        The run.
    wr : rbc.vis.fields.Write
        The snapshot write.
    o : argparse.Namespace
        Options.
    path : pathlib.Path
        PNG to write.
    spec : tuple, optional
        This frame's spectra (:meth:`Spectra.at`); drawn underneath with ``o.spectra``.

    Returns
    -------
    pathlib.Path
    """
    img, t, meta = render(info, wr, o)
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
        bars = [(r"$\theta/\theta_{\rm rms}$ (surfaces at $\pm" + f"{o.iso_T:g}" + "$)", "RdBu_r",
                 Normalize(-3, 3), None),
                (r"$w/w_{\rm rms}$ on $|\omega|/\omega_{\rm rms}=" + f"{o.iso_omega:g}" + "$", "PuOr_r",
                 Normalize(-3, 3), None)]
    else:
        lo, hi = meta["omega_range"]
        bars = [(r"$T$", o.cmap_T, Normalize(0, 1), opacity_T()),
                (r"$|\omega|^2$ (asinh scale)", o.cmap_omega, AsinhNorm(ASINH_WIDTH, lo, hi),
                 opacity_omega(*asinh_scale([lo, hi])))]
    for side, (lab, cmap, norm, op) in zip(("left", "right"), bars):
        cax = ax.inset_axes([0.08 if side == "left" else 0.58, 0.07, 0.34, 0.025])
        cb = fig.colorbar(ScalarMappable(norm, cmap), cax=cax, orientation="horizontal")
        cb.set_label(lab + (" (line: opacity)" if op is not None else ""), fontsize=7)
        cb.ax.tick_params(labelsize=6)
        if isinstance(norm, AsinhNorm):     # 0 and whole decades only (the default crowds near 0)
            top = int(np.floor(np.log10(max(abs(norm.vmin), abs(norm.vmax), 1.0))))
            ticks = [v for k in range(top + 1) for v in (-10.0 ** k, 10.0 ** k)] + [0.0]
            cb.set_ticks(sorted(v for v in ticks if norm.vmin <= v <= norm.vmax))
            cb.ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}" if abs(v) < 1e4 else f"{v:.0e}"))
            cb.ax.xaxis.set_minor_locator(NullLocator())
        if op is not None:          # the opacity transfer function, drawn across the bar
            u = np.linspace(0, 1, op.size)
            cax.plot(norm.inverse(u), 0.05 + 0.9 * op, color="k", lw=0.7)
            cax.set_ylim(0, 1)
    if o.spectra:
        if spec is not None:
            spec, ts = spec
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
    tmp = path.with_name(path.stem + ".part.png")     # renamed when complete
    fig.savefig(tmp, dpi=o.dpi)
    plt.close(fig)
    os.replace(tmp, path)
    return path


def _job(args):
    """Worker for the movie: one frame (workers are fresh interpreters, so set the style)."""
    style.use("diag")
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
    p.add_argument("--keep-negative", action="store_true",
                   help="keep points where the squared vorticity is negative (aliasing): marked in iso mode, "
                        "and opaque at the low end in volume mode, as originally")
    p.add_argument("--spectra", action="store_true", help="add the spectra at the same time underneath")
    p.add_argument("--write", type=int, default=-1, help="snapshot index (-1 = last)")
    p.add_argument("--movie", action="store_true", help="render every snapshot and assemble an MP4 (needs ffmpeg)")
    p.add_argument("--jobs", type=int, default=0, help="parallel processes (0 = as many as cores and memory allow)")
    p.add_argument("--redo", action="store_true", help="re-render frames that already exist")
    p.add_argument("--retries", type=int, default=2, help="extra rounds for frames that failed")
    p.add_argument("--frames-per-process", type=int, default=4,
                   help="frames each worker renders before it is replaced (Python 3.11+)")
    p.add_argument("--stride", type=int, default=1, help="horizontal subsampling (memory/speed)")
    p.add_argument("--nz-uniform", type=int, default=0, help="volume mode z resample (default 2 nz)")
    p.add_argument("--stretch-z", action="store_true",
                   help="volume mode: place the Chebyshev z points evenly, as the original script did")
    p.add_argument("--cmap-T", default="jet", help="volume mode colour map of T")
    p.add_argument("--cmap-omega", default="jet", help="volume mode colour map of the squared vorticity")
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
    w, h = a.window
    info = fl.RunInfo.load(run)
    writes = fl.writes(run)                              # listed once, not in every frame
    if not writes:
        raise SystemExit(f"{run}: no snapshots")
    spectra = Spectra(run) if a.spectra else None

    def spec(wr):
        return spectra.at(wr.time) if spectra is not None else None

    with VirtualDisplay((max(w, 1024), max(h, 768)), out / "xvfb.log") as display:
        probe()
        if not a.movie:
            wr = writes[a.write]
            print(compose(info, wr, a, out / f"render_{tag}.png", spec(wr)))
            return
        fdir = out / f"frames_{tag}"
        n = len(writes)
        skip = {"run", "movie", "jobs", "redo", "retries", "frames_per_process", "fps", "out"}
        parallel.prepare_frames(fdir, {k: v for k, v in vars(a).items() if k not in skip}, n, a.redo)
        items = [(fdir / f"frame_{k:05d}.png", (info, wr, a, fdir / f"frame_{k:05d}.png", spec(wr)))
                 for k, wr in enumerate(writes)]
        nx, ny, nz = fl.grid_shape(writes[0])
        points = (nx // a.stride) * (ny // a.stride) * nz
        jobs = parallel.auto_jobs(n, FRAME_BASE + BYTES_PER_POINT[a.mode] * points, a.jobs)
        # software OpenGL (llvmpipe) would otherwise start one thread per core in every worker
        os.environ.setdefault("LP_NUM_THREADS", str(parallel.worker_threads(jobs)))
        lost = parallel.render_all(_job, items, jobs, retries=a.retries, spawn=True,
                                   tasks_per_child=a.frames_per_process, before_retry=display.ensure)
    if lost:
        raise SystemExit(f"{len(lost)} of {n} frames failed (first: {lost[0].name}); no movie assembled. "
                         "Rerun the same command to render only those.")
    movie = parallel.make_movie(fdir, "frame_%05d.png", out / f"movie_{tag}.mp4", a.fps)
    if movie:
        print(movie)


if __name__ == "__main__":
    main()
