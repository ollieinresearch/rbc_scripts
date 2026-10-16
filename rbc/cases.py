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


def load_results(path: Path | str = RESULTS_CSV, dim: int | None = 3, production=True,
                 converged_only=False) -> pd.DataFrame:
    """Read ``results.csv`` and join the status columns of ``cases.csv``.

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
