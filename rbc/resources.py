"""Run history and compute cost of one run, from its own files.

Every job since this version appends a row to ``outputs/jobs.csv`` (written by the
simulation script): ranks, wall time, core-hours and the simulation time covered.
For older jobs the cost is reconstructed from the Dedalus output itself:

* each job starts a new ``analysis`` set (file handlers use ``mode="append"``), so
  one set = one job;
* the ``wall_time`` scale is seconds since the solver was created, so the largest value
  in a set is the job's run time (a lower bound: Python start-up and module imports
  before the solver existed are missing);
* the ``iteration`` and ``timestep`` scales give the number of steps and the timestep.

The number of ranks of old jobs is not stored anywhere in the output, so the estimate
uses ``MESHX * MESHY * C`` from the run's current ``params.conf``; it is wrong for jobs
that ran on a different layout. ``core_hours_source`` says which applies.

    python -m rbc.resources RUN_DIR [RUN_DIR ...]
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from rbc import io


def job_table(run: Path) -> pd.DataFrame | None:
    """The per-job log ``outputs/jobs.csv``.

    Parameters
    ----------
    run : pathlib.Path
        Run directory.

    Returns
    -------
    pandas.DataFrame or None
        One row per job (columns as ``sim/rayleigh_benard.py`` ``JOB_FIELDS``), or
        None if the run predates the log.
    """
    p = Path(run) / "outputs" / "jobs.csv"
    return pd.read_csv(p) if p.exists() else None


def segments(run: Path, folder: str = "analysis") -> pd.DataFrame:
    """One row per output set of ``folder``, i.e. per job.

    Parameters
    ----------
    run : pathlib.Path
        Run directory.
    folder : str, optional
        File-handler folder to read (``analysis`` writes most often).

    Returns
    -------
    pandas.DataFrame
        Columns ``set, t0, t1, iter0, iter1, writes, wall_s, dt_mean, dt_min``. Sets
        overlap in simulation time when a job restarted from a checkpoint older than
        the previous job's last analysis write; the overlap was computed twice and is
        counted twice here, since it cost compute time.
    """
    rows = []
    for f in io.set_files(Path(run) / folder, folder):
        try:
            with h5py.File(f, "r") as h:
                sc = h["scales"]
                t = sc["sim_time"][:]
                if t.size == 0:
                    continue
                it = sc["iteration"][:] if "iteration" in sc else np.full(t.size, -1)
                wt = sc["wall_time"][:] if "wall_time" in sc else np.full(t.size, np.nan)
                ts = sc["timestep"][:] if "timestep" in sc else np.full(t.size, np.nan)
        except OSError:
            continue
        rows.append(dict(set=int(f.stem.rsplit("_s", 1)[1]), t0=float(t[0]), t1=float(t[-1]),
                         iter0=int(it[0]), iter1=int(it[-1]), writes=int(t.size),
                         wall_s=float(np.nanmax(wt)) if np.isfinite(wt).any() else np.nan,
                         dt_mean=float(np.nanmean(ts)) if np.isfinite(ts).any() else np.nan,
                         dt_min=float(np.nanmin(ts)) if np.isfinite(ts).any() else np.nan))
    return pd.DataFrame(rows, columns=["set", "t0", "t1", "iter0", "iter1", "writes", "wall_s",
                                       "dt_mean", "dt_min"])


def params(run: Path) -> dict:
    """Numeric ``export NAME=value`` lines of a run's ``params.conf``.

    Shell arithmetic such as ``$((MESHX * MESHY))`` is evaluated for products and sums
    of variables defined earlier in the file; anything else is skipped.

    Parameters
    ----------
    run : pathlib.Path
        Run directory.

    Returns
    -------
    dict
        Variable name -> float.
    """
    conf = Path(run) / "params.conf"
    out: dict[str, float] = {}
    if not conf.exists():
        return out
    for m in re.finditer(r"^\s*export\s+(\w+)=([^\s#]+)", conf.read_text(), re.M):
        name, val = m.group(1), m.group(2).strip("\"'")
        a = re.fullmatch(r"\$\(\((.+)\)\)", val)
        if a:
            expr = re.sub(r"\$?\{?([A-Za-z_]\w*)\}?", lambda v: repr(out.get(v.group(1), "x")), a.group(1))
            if re.fullmatch(r"[\d.\s*+\-/()]+", expr):
                out[name] = float(eval(expr))       # digits and operators only
            continue
        try:
            out[name] = float(val)
        except ValueError:
            pass
    return out


def grid_shape(run: Path) -> tuple[int, int, int] | None:
    """Coefficient shape (Nx, Ny, Nz) of a run.

    Taken from the spectra files' attributes, else from the newest coefficient-space
    checkpoint, else from ``params.conf`` (RES_HOR, RES_Z).

    Parameters
    ----------
    run : pathlib.Path
        Run directory.

    Returns
    -------
    tuple of int or None
    """
    run = Path(run)
    for f in sorted((run / "spectra").glob("spectra_s*.h5")):
        try:
            with h5py.File(f, "r") as h:
                return tuple(int(h.attrs[k]) for k in ("Nx", "Ny", "Nz"))
        except (OSError, KeyError):
            continue
    states = io.set_files(run / "state", "state")
    if states:
        from rbc import regrid
        try:
            s = regrid.coeff_shape(states[-1])
            if s:
                return s
        except (OSError, KeyError):
            pass
    p = params(run)
    if "RES_Z" in p and "RES_HOR" in p:
        return int(p["RES_HOR"]), int(p["RES_HOR"]), int(p["RES_Z"])
    return None


def summary(run: Path) -> dict:
    """History and cost of one run, for ``results.csv``.

    Parameters
    ----------
    run : pathlib.Path
        Run directory.

    Returns
    -------
    dict
        ``n_jobs``, ``t_run`` (simulation time computed, overlaps included), ``iterations``,
        ``dt_mean`` (time-weighted mean timestep), ``wall_hours``, ``ranks`` (latest),
        ``core_hours``, ``core_hours_source`` (``jobs.csv``, ``estimate`` or a mix),
        ``nx``, ``ny``, ``nz``.
    """
    run = Path(run)
    seg = segments(run)
    jobs = job_table(run)
    p = params(run)
    ranks_now = p.get("N") or (p.get("MESHX", np.nan) * p.get("MESHY", np.nan))
    cpr = p.get("C", 1.0)
    out = dict(n_jobs=int(len(seg)))
    if len(seg):
        span = (seg["t1"] - seg["t0"]).clip(lower=0)
        steps = (seg["iter1"] - seg["iter0"]).clip(lower=0)
        out.update(t_run=float(span.sum()), iterations=int(steps.sum()),
                   dt_mean=float(span.sum() / steps.sum()) if steps.sum() else np.nan,
                   wall_hours=float(seg["wall_s"].sum() / 3600))
    ch, est, logged = 0.0, 0, 0
    for _, s in seg.iterrows():
        hit = None
        if jobs is not None and len(jobs):
            tol = 1e-6 + 2 * (s["dt_mean"] if np.isfinite(s["dt_mean"]) else 0)
            m = (jobs["t_start"] <= s["t0"] + tol) & (jobs["t_end"] >= s["t1"] - tol)
            if m.any():
                hit = jobs[m].iloc[-1]
        if hit is not None:
            ch += float(hit["core_hours"])
            logged += 1
        elif np.isfinite(s["wall_s"]) and np.isfinite(ranks_now):
            ch += s["wall_s"] * ranks_now * cpr / 3600
            est += 1
    if logged or est:
        out["core_hours"] = ch
        out["core_hours_source"] = ("jobs.csv" if not est else "estimate" if not logged
                                    else f"jobs.csv ({logged}) + estimate ({est})")
    if jobs is not None and len(jobs):
        out["ranks"] = int(jobs["ranks"].iloc[-1])
    elif np.isfinite(ranks_now):
        out["ranks"] = int(ranks_now)
    shape = grid_shape(run)
    if shape:
        out.update(nx=shape[0], ny=shape[1], nz=shape[2])
    return out


def get_parser() -> argparse.ArgumentParser:
    """Command-line options.

    Returns
    -------
    argparse.ArgumentParser
    """
    p = argparse.ArgumentParser(prog="python -m rbc.resources",
                                description="Print the job history and cost of runs.")
    p.add_argument("runs", type=Path, nargs="+", help="run directories")
    p.add_argument("--jobs", action="store_true", help="also print every job (analysis set)")
    return p


def main(argv=None) -> None:
    """Print :func:`summary` (and optionally :func:`segments`) for each run.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments (default: ``sys.argv[1:]``).
    """
    a = get_parser().parse_args(argv)
    with pd.option_context("display.width", 160, "display.max_columns", 20):
        for run in a.runs:
            print(run)
            for k, v in summary(run).items():
                print(f"  {k:<18} {v}")
            if a.jobs:
                print(segments(run).to_string(index=False))


if __name__ == "__main__":
    main()
