"""Reading Dedalus file-handler output across restarts.

Every job writes new sets (``analysis_s1.h5``, ``analysis_s2.h5``, ...). A restart
resumes from the last *state* write, which is earlier than the last analysis write, so
consecutive sets overlap in time. The later set is the trajectory that continues, so it
wins: earlier sets are truncated at the first time of each later set. Only the
surviving rows are read from disk.
"""
from __future__ import annotations

import re
from pathlib import Path

import h5py
import numpy as np

_SET_RE = re.compile(r"_s(\d+)\.h5$")


def set_files(folder: Path | str, name: str | None = None) -> list[Path]:
    """The set files of one file handler, in set order.

    Parameters
    ----------
    folder : str or pathlib.Path
        File-handler folder, e.g. ``RUN/analysis``.
    name : str, optional
        Handler name (default: the folder's name).

    Returns
    -------
    list of pathlib.Path
        ``<name>_s<N>.h5`` files sorted by N (not alphabetically: s10 after s9).
    """
    folder = Path(folder)
    name = name or folder.name
    files = [p for p in folder.glob(f"{name}_s*.h5") if _SET_RE.search(p.name)]
    return sorted(files, key=lambda p: int(_SET_RE.search(p.name).group(1)))


def plan_reads(files, t_min=-np.inf, t_max=np.inf):
    """Which rows of which files make one monotonic time series.

    Parameters
    ----------
    files : list of pathlib.Path
        Set files in set order (:func:`set_files`).
    t_min, t_max : float, optional
        Time window.

    Returns
    -------
    list of tuple
        ``(file, i0, i1)``: rows ``i0:i1`` of each file, later sets winning where they
        overlap.
    """
    segs = []
    for f in files:
        with h5py.File(f, "r") as h:
            t = h["scales/sim_time"][:]
        if t.size == 0:
            continue
        kept = []
        for g, i0, i1, tg in segs:
            i1 = min(i1, int(np.searchsorted(tg, t[0], "left")))
            if i1 > i0:
                kept.append((g, i0, i1, tg))
        segs = kept + [(f, 0, t.size, t)]
    out = []
    for f, i0, i1, t in segs:
        a = max(i0, int(np.searchsorted(t, t_min, "left")))
        b = min(i1, int(np.searchsorted(t, t_max, "right")))
        if b > a:
            out.append((f, a, b))
    return out


def _squeeze(a: np.ndarray) -> np.ndarray:
    """Drop singleton axes after the first (time) axis."""
    keep = [0] + [i for i in range(1, a.ndim) if a.shape[i] != 1]
    return a.reshape([a.shape[i] for i in keep])


def load_series(folder, tasks, t_min=-np.inf, t_max=np.inf, name=None, optional=()):
    """Time series of file-handler tasks across all sets.

    Parameters
    ----------
    folder : str or pathlib.Path
        File-handler folder.
    tasks : list of str
        Task names to read.
    t_min, t_max : float, optional
        Time window.
    name : str, optional
        Handler name (default: the folder's name).
    optional : iterable of str, optional
        Tasks that may be missing from the files (skipped instead of an error).

    Returns
    -------
    t : numpy.ndarray
        Sample times.
    data : dict
        Task -> array of shape ``(nt, ...)`` with singleton spatial axes removed.

    Raises
    ------
    FileNotFoundError
        If no data falls in the window.
    KeyError
        If a required task is missing.
    """
    plan = plan_reads(set_files(folder, name), t_min, t_max)
    if not plan:
        raise FileNotFoundError(f"no data in {folder} for t in [{t_min}, {t_max}]")
    with h5py.File(plan[0][0], "r") as h:
        present = set(h["tasks"])
    missing = [k for k in tasks if k not in present and k not in optional]
    if missing:
        raise KeyError(f"{folder}: missing tasks {missing}")
    tasks = [k for k in tasks if k in present]

    ts, chunks = [], {k: [] for k in tasks}
    for f, a, b in plan:
        with h5py.File(f, "r") as h:
            ts.append(h["scales/sim_time"][a:b])
            for k in tasks:
                chunks[k].append(h["tasks"][k][a:b])
    return np.concatenate(ts), {k: _squeeze(np.concatenate(v)) for k, v in chunks.items()}


def task_grid(folder, task, axis=-1, name=None) -> np.ndarray | None:
    """Grid that Dedalus attached to a task as an HDF5 dimension scale.

    Parameters
    ----------
    folder : str or pathlib.Path
        File-handler folder.
    task : str
        Task name.
    axis : int, optional
        Axis of the task dataset (default: the last, z for a profile).
    name : str, optional
        Handler name.

    Returns
    -------
    numpy.ndarray or None
        The grid, or None if no scale is attached.
    """
    files = set_files(folder, name)
    with h5py.File(files[0], "r") as h:
        d = h["tasks"][task]
        dims = d.dims[axis % d.ndim]
        return dims[0][:] if len(dims) else None


def time_mean(t: np.ndarray, x: np.ndarray, axis: int = 0) -> np.ndarray:
    """Simpson-rule time average.

    Parameters
    ----------
    t : numpy.ndarray
        Sample times.
    x : numpy.ndarray
        Samples, time along ``axis``.
    axis : int, optional
        Time axis of ``x``.

    Returns
    -------
    numpy.ndarray
        Average over time.
    """
    from scipy.integrate import simpson
    return simpson(x, x=t, axis=axis) / (t[-1] - t[0])
