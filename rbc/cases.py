"""Parameter registry and naming conventions shared by every script.

``data/cases.csv`` (edited by hand, one row per (dim, Ra, Pr) point) says what *should*
exist and its status. ``data/results.csv`` (written by :mod:`rbc.collect`, one row per
run directory) holds the computed numbers. They join on ``(dim, log_ra, log_pr)``.

Everything is keyed by log10 exponents on exact half-decades, never by floats like
3.16e-3 or by folder strings like ``pr3e-3``.

The data directory is ``data/`` in the repository, or ``$RBC_DATA`` if set (the
documentation's example gallery uses this to point at synthetic data).
"""
from __future__ import annotations

import os
import re
from pathlib import Path

import numpy as np
import pandas as pd

#: Slices are labelled by their Pr at Ra = 10**RA_REF.
RA_REF = 5.0
#: Directory holding ``cases.csv`` and ``results.csv``.
DATA_DIR = Path(os.environ.get("RBC_DATA", Path(__file__).resolve().parents[1] / "data"))
CASES_CSV = DATA_DIR / "cases.csv"
RESULTS_CSV = DATA_DIR / "results.csv"

#: Allowed values of the ``status`` column of ``cases.csv``.
STATUSES = ("finished", "running", "planned")

#: Slice slopes a in Pr = Pr_0 (Ra / 10**RA_REF)**a.
FAMILIES = (-0.5, -1.0)


def slice_c(log_ra, log_pr, a=-0.5):
    """Label of the slice through a point: c = log10 Pr_0.

    Along Pr = Pr_0 (Ra / 10**RA_REF)**a, the quantity log10 Pr - a (log10 Ra - RA_REF)
    is constant; it is rounded to the nearest half-decade.

    Parameters
    ----------
    log_ra, log_pr : float or array_like
        log10 Ra and log10 Pr of the point(s).
    a : float, optional
        Slice slope (-0.5 for Pr ∝ Ra^-1/2, -1.0 for Pr ∝ Ra^-1).

    Returns
    -------
    numpy.ndarray or float
        c for each point.
    """
    return np.round((np.asarray(log_pr) - a * (np.asarray(log_ra) - RA_REF)) * 2) / 2


def slices(df: pd.DataFrame, a=-0.5, min_points=3) -> dict[float, pd.DataFrame]:
    """All slices of slope ``a`` present in a table.

    Parameters
    ----------
    df : pandas.DataFrame
        Rows with ``log_ra`` and ``log_pr`` columns (cases or results).
    a : float, optional
        Slice slope.
    min_points : int, optional
        A slice needs at least this many distinct Ra values to be returned.

    Returns
    -------
    dict
        c = log10 Pr_0 -> the rows on that slice, sorted by Ra; largest c first.
    """
    c = slice_c(df["log_ra"], df["log_pr"], a)
    out = {}
    for cv in sorted(set(c), reverse=True):
        g = df[c == cv].sort_values("log_ra")
        if g["log_ra"].nunique() >= min_points:
            out[float(cv)] = g
    return out


def half_decade(x: float) -> float:
    """log10 of a folder value, rounded to the nearest half-decade.

    Parameters
    ----------
    x : float or str
        Value as in a folder name, e.g. ``3e-3``.

    Returns
    -------
    float
        e.g. -2.5.
    """
    return round(np.log10(float(x)) * 2) / 2


def folder_value(e: float) -> str:
    """Folder-name form of a half-decade exponent (inverse of :func:`half_decade`).

    Parameters
    ----------
    e : float
        log10 of the value, a multiple of 0.5.

    Returns
    -------
    str
        e.g. ``-2.5 -> '3e-3'``, ``8 -> '1e8'``.
    """
    e = float(e)
    k = int(np.floor(e))
    mant = 1 if e == k else 3
    return f"{mant}e{k}"


_RUN_RE = re.compile(r"ra(?P<ra>[\d.]+e[-+]?\d+)/pr(?P<pr>[\d.]+e[-+]?\d+)/(?P<res>[^/]+)")


def parse_run_path(path) -> dict | None:
    """Parameters encoded in a run directory's path.

    Parameters
    ----------
    path : str or pathlib.Path
        e.g. ``.../3d/ra1e8/pr3e-3/384_Gam2_2``.

    Returns
    -------
    dict or None
        ``dim`` (2 if ``/2d/`` is in the path, else 3), ``log_ra``, ``log_pr``,
        ``nz`` (first number of the label), ``nx`` (second number, unless the label is
        of the ``<nz>_Gam2_2`` form, which encodes the box rather than nx) and
        ``label``; None if the path does not match.
    """
    p = Path(path).as_posix()
    m = _RUN_RE.search(p)
    if not m:
        return None
    dim = 2 if "/2d/" in p else 3
    res = m["res"]
    nums = [int(t) for t in res.split("_") if t.isdigit()]
    nz = nums[0] if nums else None
    nx = nums[1] if len(nums) > 1 and "Gam" not in res else None
    return dict(dim=dim, log_ra=half_decade(m["ra"]), log_pr=half_decade(m["pr"]),
                nz=nz, nx=nx, label=res)


def load_cases(path: Path | str = CASES_CSV, dim: int | None = None) -> pd.DataFrame:
    """Read and validate ``cases.csv``.

    Parameters
    ----------
    path : str or pathlib.Path, optional
        The registry file. Lines starting with ``#`` are comments.
    dim : int, optional
        Keep only this dimension.

    Returns
    -------
    pandas.DataFrame
        The registry with ``res_check`` as bool and added columns ``c`` (slice label for
        a = -1/2), ``Ra`` and ``Pr``.

    Raises
    ------
    ValueError
        On an unknown status or a duplicated (dim, Ra, Pr) point.
    """
    df = pd.read_csv(path, comment="#", skipinitialspace=True)
    df.columns = df.columns.str.strip()
    df["status"] = df["status"].str.strip().str.lower()
    bad = set(df["status"]) - set(STATUSES)
    if bad:
        raise ValueError(f"{path}: unknown status {bad}; use one of {STATUSES}")
    df["res_check"] = df["res_check"].fillna(0).astype(bool)
    df["c"] = slice_c(df["log_ra"], df["log_pr"])
    df["Ra"] = 10.0 ** df["log_ra"]
    df["Pr"] = 10.0 ** df["log_pr"]
    if df.duplicated(["dim", "log_ra", "log_pr"]).any():
        dup = df[df.duplicated(["dim", "log_ra", "log_pr"], keep=False)]
        raise ValueError(f"{path}: duplicate parameter points\n{dup}")
    return df if dim is None else df[df["dim"] == dim].reset_index(drop=True)


#: Column names accepted for the three columns of a simple results file.
SIMPLE_COLUMNS = {"ra": "Ra", "pr": "Pr", "nusselt": "Nu", "nu": "Nu"}


def _first_row_is_header(path: Path | str) -> bool:
    """Whether the first non-comment line of a CSV contains non-numeric fields."""
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                try:
                    [float(v) for v in line.split(",") if v.strip()]
                    return False
                except ValueError:
                    return True
    return False


def is_simple_results(path: Path | str) -> bool:
    """Whether a results file is the simple three-column ``Ra,Pr,Nusselt`` format.

    A file counts as simple if it has exactly three columns and either no header or a
    header of the names in :data:`SIMPLE_COLUMNS` (any case).

    Parameters
    ----------
    path : str or pathlib.Path
        The file.

    Returns
    -------
    bool
    """
    if not _first_row_is_header(path):
        return pd.read_csv(path, header=None, comment="#", nrows=1).shape[1] == 3
    cols = [c.strip().lower() for c in pd.read_csv(path, comment="#", nrows=0).columns]
    return len(cols) == 3 and all(c in SIMPLE_COLUMNS for c in cols)


def _snap_half_decade(logv: np.ndarray, name: str, tol: float = 0.05) -> np.ndarray:
    """Round log10 values to the nearest half-decade when within ``tol`` of it, so that
    3e-3 and 3.16e-3 both become -2.5; values further away are kept and reported."""
    snapped = np.round(logv * 2) / 2
    off = np.abs(logv - snapped) > tol
    if off.any():
        print(f"note: {off.sum()} {name} value(s) are not on half-decades and are kept as given: "
              f"{sorted({float(v) for v in np.round(10 ** logv[off], 6)})}")
    return np.where(off, logv, snapped)


def read_simple_results(path: Path | str, dim: int = 3) -> pd.DataFrame:
    """Read a three-column ``Ra,Pr,Nusselt`` file into the ``results.csv`` layout.

    The columns are Ra, Pr and Nu, in that order, with or without a header line. Values
    are the numbers themselves (1e5, 0.01); if every Ra is below 100, all values are taken
    as log10 instead (convection needs Ra > 1708, so this is unambiguous). log10 Ra and
    log10 Pr within 0.05 of a half-decade are rounded to it. Every point is treated as a
    finished production run; there are no standard errors, so fits are unweighted.

    Parameters
    ----------
    path : str or pathlib.Path
        The file.
    dim : int, optional
        Dimension assigned to every row.

    Returns
    -------
    pandas.DataFrame
        Columns ``dim, log_ra, log_pr, Ra, Pr, Nu, status, res_check, production``.

    Raises
    ------
    ValueError
        On non-positive or missing values, or a (Ra, Pr) point given twice.
    """
    header = _first_row_is_header(path)
    raw = pd.read_csv(path, comment="#", header=0 if header else None, skipinitialspace=True)
    if header:
        raw.columns = [SIMPLE_COLUMNS[c.strip().lower()] for c in raw.columns]
        raw = raw[["Ra", "Pr", "Nu"]]
    else:
        raw.columns = ["Ra", "Pr", "Nu"]
    raw = raw.apply(pd.to_numeric, errors="coerce")
    bad = raw.isna().any(axis=1)
    if bad.any():
        raise ValueError(f"{path}: non-numeric or missing values in rows {list(raw.index[bad] + 1 + header)}")
    ra, pr = raw["Ra"].to_numpy(float), raw["Pr"].to_numpy(float)
    if (ra < 100).all():
        print(f"note: every Ra in {path} is below 100, so Ra and Pr are read as log10 values")
        lra, lpr = ra, pr
    else:
        if (ra <= 0).any() or (pr <= 0).any():
            raise ValueError(f"{path}: Ra and Pr must be positive")
        lra, lpr = np.log10(ra), np.log10(pr)
    df = pd.DataFrame(dict(dim=dim, log_ra=_snap_half_decade(lra, "Ra"), log_pr=_snap_half_decade(lpr, "Pr"),
                           Nu=raw["Nu"].to_numpy(float)))
    dup = df.duplicated(["log_ra", "log_pr"], keep=False)
    if dup.any():
        raise ValueError(f"{path}: these (log10 Ra, log10 Pr) points appear more than once; keep one row each:\n"
                         f"{df.loc[dup, ['log_ra', 'log_pr', 'Nu']].to_string()}")
    df["Ra"], df["Pr"] = 10.0 ** df["log_ra"], 10.0 ** df["log_pr"]
    df["status"], df["res_check"], df["production"] = "finished", False, True
    print(f"read {len(df)} points from simple results file {path}")
    return df


def load_results(path: Path | str = RESULTS_CSV, dim: int | None = 3, production=True,
                 converged_only=False) -> pd.DataFrame:
    """Read ``results.csv`` and join the status columns of ``cases.csv``.

    A three-column ``Ra,Pr,Nusselt`` file is also accepted and detected automatically
    (:func:`read_simple_results`); its points are all treated as finished, and
    ``cases.csv`` is not consulted.

    Parameters
    ----------
    path : str or pathlib.Path, optional
        The results table.
    dim : int or None, optional
        Keep only this dimension (None: all).
    production : bool, optional
        Keep only the production run of each point (see :func:`rbc.collect.mark_production`).
    converged_only : bool, optional
        Keep only rows with ``converged`` true.

    Returns
    -------
    pandas.DataFrame
        Rows with added ``Ra``, ``Pr``, and, when ``cases.csv`` exists, ``status`` and
        ``res_check``.
    """
    if is_simple_results(path):
        return read_simple_results(path, dim if dim is not None else 3)
    df = pd.read_csv(path)
    if dim is not None:
        df = df[df["dim"] == dim]
    if production and "production" in df:
        df = df[df["production"].astype(bool)]
    if converged_only and "converged" in df:
        df = df[df["converged"].astype(bool)]
    df = df.copy()
    df["Ra"], df["Pr"] = 10.0 ** df["log_ra"], 10.0 ** df["log_pr"]
    if CASES_CSV.exists():
        c = load_cases()[["dim", "log_ra", "log_pr", "status", "res_check"]]
        df = df.merge(c, on=["dim", "log_ra", "log_pr"], how="left")
    return df.reset_index(drop=True)


def require_columns(df: pd.DataFrame, cols, what: str) -> None:
    """Stop with a clear message if a figure needs columns the results do not have.

    Parameters
    ----------
    df : pandas.DataFrame
        Loaded results.
    cols : iterable of str
        Needed columns.
    what : str
        Name of the figure, for the message.

    Raises
    ------
    SystemExit
        If any column is missing.
    """
    missing = [c for c in cols if c not in df]
    if missing:
        raise SystemExit(f"{what} needs {', '.join(missing)}, which these results do not have "
                         "(a simple Ra,Pr,Nusselt file only supports the Nu figures)")