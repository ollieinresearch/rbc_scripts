"""Snapshot access and the normalization shared by the field visualizations.

Fields are normalized so that colour scales are fixed across frames and across Pr:

* temperature: theta = T - Tbar(z), shown as theta / theta_rms(z);
* velocity: w / w_rms(z) (and u, v by the horizontal rms per component);
* vorticity: ``|ω|`` / omega_rms with omega_rms = sqrt((Nu - 1) / Pr), since in
  free-fall units the exact energy balance gives <``|∇u|²``> = <``|ω|²``> =
  (Nu - 1)/Pr (no-slip plates, periodic sides).

The profiles come from ``outputs/profiles.npz`` (time means); without it, statistics of
the plane itself are used. ``|ω|²`` is written as a pointwise product truncated to
the output grid, so where it is small it can dip below zero; those points are clipped
unless ``keep_negative``.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import h5py
import numpy as np
from scipy.interpolate import BarycentricInterpolator

from rbc import io
from rbc.cases import parse_run_path


@dataclass
class Write:
    """One snapshot write: the file, the index within it, its time and write number."""
    file: Path
    index: int
    time: float
    number: int


def writes(run: Path, folder="snapshots") -> list[Write]:
    """Every write of a file handler across sets, restart overlaps removed.

    Parameters
    ----------
    run : pathlib.Path
        Run directory.
    folder : str, optional
        File-handler folder.

    Returns
    -------
    list of Write
        In time order.
    """
    out = []
    for f, a, b in io.plan_reads(io.set_files(run / folder, folder)):
        with h5py.File(f, "r") as h:
            t = h["scales/sim_time"][a:b]
            n = h["scales/write_number"][a:b]
        out += [Write(f, a + i, float(t[i]), int(n[i])) for i in range(b - a)]
    return out


def _task(h, name):
    """Dataset and component index of a field, for current and older file layouts.

    Parameters
    ----------
    h : h5py.File
        Snapshot file.
    name : {"T", "u", "v", "w", "vort"}
        Field.

    Returns
    -------
    tuple
        ``(dataset, component)``, component None for scalars.
    """
    t = h["tasks"]
    if name == "T":
        return t["temperature"], None
    if name == "vort":
        return t["vorticity"], None
    comp = {"u": 0, "v": 1, "w": 2}[name]
    if name in t and t[name].ndim == 4:
        return t[name], None
    return t["u"], comp


def grids(h, name="T", L=(2.0, 2.0)):
    """x, y, z grids of a snapshot field.

    Parameters
    ----------
    h : h5py.File
        Snapshot file.
    name : str, optional
        Field whose dimension scales are read.
    L : tuple of float, optional
        Box size, for uniform x, y grids when no scale is attached.

    Returns
    -------
    list of numpy.ndarray
        ``[x, y, z]``.

    Raises
    ------
    ValueError
        If the file has no z scale.
    """
    d, comp = _task(h, name)
    off = d.ndim - 3
    ax = []
    for k in range(3):
        dims = d.dims[off + k]
        if len(dims):
            ax.append(dims[0][:])
        elif k < 2:
            ax.append(np.linspace(0, L[k], d.shape[off + k], endpoint=False))
        else:
            raise ValueError("no z scale in snapshot file")
    return ax


def plane(h, name, idx, z0, k=8):
    """Horizontal plane of a field at any height.

    Barycentric interpolation over the ``k`` nearest grid levels, which agrees with the
    Chebyshev interpolant to O(h^k) and reads only those levels from disk.

    Parameters
    ----------
    h : h5py.File
        Snapshot file.
    name : str
        Field.
    idx : int
        Write index within the file.
    z0 : float
        Height in [0, 1].
    k : int, optional
        Number of levels used.

    Returns
    -------
    numpy.ndarray
        (nx, ny) values.
    """
    d, comp = _task(h, name)
    z = grids(h, name)[2]
    j = np.argsort(np.abs(z - z0))[:k]
    lo, hi = j.min(), j.max() + 1
    sel = (idx,) + ((comp,) if comp is not None else ()) + (slice(None), slice(None), slice(lo, hi))
    block = d[sel].astype(np.float64)
    return BarycentricInterpolator(z[lo:hi], block, axis=-1)(z0)


def volume(h, name, idx, dtype=np.float32):
    """A whole 3D field.

    Parameters
    ----------
    h : h5py.File
        Snapshot file.
    name : str
        Field.
    idx : int
        Write index.
    dtype : numpy.dtype, optional
        Output type (single precision halves the memory).

    Returns
    -------
    numpy.ndarray
        (nx, ny, nz) values.
    """
    d, comp = _task(h, name)
    sel = (idx,) + ((comp,) if comp is not None else ())
    return d[sel].astype(dtype)


@dataclass
class RunInfo:
    """What a visualization needs to know about a run.

    Attributes
    ----------
    path : pathlib.Path
        Run directory.
    Ra, Pr : float
        Parameters (from the path).
    Nu : float
        From ``outputs/nusselt.json`` (NaN if absent).
    prof : dict
        Time-mean profiles from ``outputs/profiles.npz``, z in [0, 1].
    """
    path: Path
    Ra: float
    Pr: float
    Nu: float
    prof: dict = field(default_factory=dict)

    @classmethod
    def load(cls, run: Path) -> "RunInfo":
        """Read a run's parameters, Nu and profiles.

        Parameters
        ----------
        run : pathlib.Path
            Run directory.

        Returns
        -------
        RunInfo
        """
        run = Path(run).resolve()
        p = parse_run_path(run) or {}
        nj = run / "outputs" / "nusselt.json"
        Nu = json.loads(nj.read_text())["Nu"] if nj.exists() else np.nan
        prof = {}
        pz = run / "outputs" / "profiles.npz"
        if pz.exists():
            with np.load(pz) as d:
                prof = {k: d[k] for k in d.files}
                prof["z"] = prof["z"] + 0.5                    # stored centred; use [0, 1]
        return cls(run, 10.0 ** p.get("log_ra", np.nan), 10.0 ** p.get("log_pr", np.nan), Nu, prof)

    @property
    def omega_rms(self) -> float:
        """float: rms vorticity sqrt((Nu - 1) / Pr) in free-fall units."""
        return float(np.sqrt((self.Nu - 1) / self.Pr))

    @property
    def delta_T(self) -> float:
        """float: thermal boundary-layer thickness 1/(2 Nu)."""
        return 1 / (2 * self.Nu)

    def height(self, h) -> float:
        """A named or numeric height.

        Parameters
        ----------
        h : str or float
            ``"mid"`` (1/2), ``"bl"`` (delta_T/2 = 1/(4 Nu)), ``"bl-edge"`` (delta_T), or
            a number.

        Returns
        -------
        float
        """
        if isinstance(h, (int, float)):
            return float(h)
        return {"mid": 0.5, "bl": self.delta_T / 2, "bl-edge": self.delta_T}.get(h) or float(h)

    def profile(self, key, z):
        """A time-mean profile interpolated to height(s) z (None if unavailable).

        Parameters
        ----------
        key : str
            Profile name in ``profiles.npz`` (``T``, ``w_rms``, ``uh_rms``, ``theta_rms``).
        z : float or numpy.ndarray
            Heights in [0, 1].

        Returns
        -------
        numpy.ndarray or None
        """
        return np.interp(z, self.prof["z"], self.prof[key]) if key in self.prof else None

    def label(self) -> str:
        """str: ``Ra = 10^..., Pr = 10^...`` for titles."""
        from rbc.style import pow10
        return rf"$Ra={pow10(np.log10(self.Ra))}$, $Pr={pow10(np.log10(self.Pr))}$"


def normalized(info: RunInfo, name, data, z, keep_negative=False):
    """Apply the normalization described in the module docstring.

    Parameters
    ----------
    info : RunInfo
        The run.
    name : {"T", "u", "v", "w", "vort"}
        Field.
    data : numpy.ndarray
        A plane (nx, ny) or volume (nx, ny, nz).
    z : float or numpy.ndarray
        Height of the plane, or heights of the volume's last axis.
    keep_negative : bool, optional
        For vorticity: keep ``|ω|²`` < 0 points (signed) instead of clipping them.

    Returns
    -------
    numpy.ndarray
        Normalized field.
    """
    if name == "T":
        tbar = info.profile("T", z)
        th_rms = info.profile("theta_rms", z)
        if tbar is None:
            tbar = data.mean(axis=(0, 1))
        theta = data - tbar
        if th_rms is None:
            th_rms = np.sqrt((theta ** 2).mean(axis=(0, 1)))
        return theta / th_rms
    if name in ("u", "v", "w"):
        key = "w_rms" if name == "w" else "uh_rms"
        r = info.profile(key, z)
        if r is not None and name != "w":
            r = r / np.sqrt(2)
        if r is None:
            r = np.sqrt((data ** 2).mean(axis=(0, 1)))
        return data / r
    if name == "vort":
        d = data if keep_negative else np.maximum(data, 0)
        return np.sign(d) * np.sqrt(np.abs(d)) / info.omega_rms
    raise ValueError(name)
