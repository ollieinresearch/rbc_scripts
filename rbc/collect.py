"""Collect every run's results into ``data/results.csv`` (one row per run directory).

    python -m rbc.collect [ROOT] [--recalculate] [--jobs N] [--no-plots]

ROOT defaults to ``$DATA_ROOT``, else ``$SCRATCH/ollie_rb_data/3d``. It can be any copy
of the run tree, for example the backup on project space:

    python -m rbc.collect /project/def-goluskin/ollie/ollie_rb_data/3d

Runs are directories ``ra*/pr*/<label>`` with an ``analysis/`` folder; paths matching
``--exclude`` (default: ``_old`` or ``test``) are skipped.

For each run, ``outputs/{nusselt,profiles,spectra}.json`` are reused if they are newer
than the data; ``--recalculate`` redoes them, with the averaging window starting at
AVG_TIME from the run's ``params.conf`` (a number or ``auto``) unless ``--t0`` is given.
Runs are processed in parallel, one process per run.

Each row is tagged with the cluster it was collected on. Rows from other clusters
already in the output are kept, so copying another cluster's ``results.csv`` in with
``--merge`` before collecting combines them. ``results.parquet`` (same table, column
types preserved) is written too when pyarrow is available.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import socket
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

from rbc import cases as rc

CLUSTER = os.environ.get("CC_CLUSTER") or socket.gethostname().split(".")[0].rstrip("0123456789-")
#: Identities whose diagnostics enter the convergence columns (as in rbc.analysis.nusselt).
NU_KEYS = ("avg_wT", "avg_grad_T_sq")


def default_root() -> Path:
    """The run tree to collect when none is given: ``$DATA_ROOT`` or
    ``$SCRATCH/ollie_rb_data/3d``.

    Returns
    -------
    pathlib.Path
    """
    if os.environ.get("DATA_ROOT"):
        return Path(os.environ["DATA_ROOT"])
    return Path(os.environ.get("SCRATCH", ".")) / "ollie_rb_data" / "3d"


def find_runs(root: Path, exclude: str) -> list[Path]:
    """Every run directory below ``root``.

    Parameters
    ----------
    root : pathlib.Path
        Top of the run tree.
    exclude : str
        Regular expression; runs whose path relative to ``root`` matches are skipped.

    Returns
    -------
    list of pathlib.Path
        Sorted run directories (those with an ``analysis/`` folder and a parseable path).
    """
    pat = re.compile(exclude) if exclude else None
    runs = []
    for a in root.rglob("analysis"):
        run = a.parent
        if a.is_dir() and rc.parse_run_path(run) and not (pat and pat.search(run.relative_to(root).as_posix())):
            runs.append(run)
    return sorted(runs)


def avg_time(run: Path) -> float | str:
    """Averaging start (AVG_TIME) from a run's ``params.conf``.

    Parameters
    ----------
    run : pathlib.Path
        Run directory.

    Returns
    -------
    float or str
        The value, ``"auto"``, or 0.0 if unset.
    """
    conf = run / "params.conf"
    if conf.exists():
        m = re.search(r"^\s*export\s+AVG_TIME=([-+\w.]+)", conf.read_text(), re.M)
        if m:
            return "auto" if m.group(1) == "auto" else float(m.group(1))
    return 0.0


def _stale(out: Path, data_dir: Path) -> bool:
    """True if ``out`` is missing or older than any ``.h5`` file in ``data_dir``."""
    if not out.exists():
        return True
    newest = max((p.stat().st_mtime for p in data_dir.glob("*.h5")), default=0)
    return newest > out.stat().st_mtime


def process_run(run: Path, recalc: bool, t0, plots: bool) -> tuple[dict | None, str]:
    """Bring one run's ``outputs/*.json`` up to date and return its row.

    Parameters
    ----------
    run : pathlib.Path
        Run directory.
    recalc : bool
        Recompute even if the outputs are newer than the data.
    t0 : float, "auto" or None
        Averaging start; None means the run's AVG_TIME.
    plots : bool
        Also write the per-run diagnostic PDFs.

    Returns
    -------
    row : dict or None
        The ``results.csv`` row, or None if the run failed.
    message : str
        One line for the log (or a traceback).
    """
    from rbc.analysis import nusselt, profiles

    out = run / "outputs"
    out.mkdir(exist_ok=True)
    t0 = avg_time(run) if t0 is None else t0
    log = []
    try:
        if recalc or _stale(out / "nusselt.json", run / "analysis"):
            res, series = nusselt.compute(run, t0)
            (out / "nusselt.json").write_text(json.dumps(res, indent=1))
            nusselt.write_info(res, out / "info.txt")
            if plots:
                nusselt.plot(res, series, out / "info.pdf")
            log.append("nusselt")
        t_win = json.loads((out / "nusselt.json").read_text())["t0"]
        if (run / "horizontal_analysis").exists() and (recalc or _stale(out / "profiles.json", run / "horizontal_analysis")):
            prof = profiles.compute(run, t_win)
            np.savez(out / "profiles.npz", **{k: v for k, v in prof.items() if k != "diag"})
            (out / "profiles.json").write_text(json.dumps(prof["diag"], indent=1))
            if plots:
                profiles.plot(prof, out / "profiles.pdf")
            log.append("profiles")
        if (run / "spectra").exists() and (recalc or _stale(out / "spectra.json", run / "spectra")):
            from rbc.spectra import cumulative
            if plots:
                cumulative.run(run, "auto", t_win, 0.0)
            else:
                _spectra_json_only(run, t_win)
            log.append("spectra")
    except Exception:
        return None, f"FAILED {run}\n{traceback.format_exc(limit=3)}"
    return row(run), f"{run}: {' '.join(log) or 'up to date'}"


def _spectra_json_only(run: Path, t0: float) -> None:
    """Write ``outputs/spectra.json`` without the spectra figures."""
    from rbc.spectra import core as cs, cumulative as cu
    time, modes, data, attrs = cs.load_spectra(cs.spectra_files(run))
    time, data = cu.select_window(time, data, t0, 0.0)
    cu.write_json(run / "outputs" / "spectra.json", time, cu.summarize(time, modes, data, attrs), attrs)


def _read(p: Path) -> dict:
    """A JSON file as a dict ({} if missing)."""
    return json.loads(p.read_text()) if p.exists() else {}


def row(run: Path) -> dict:
    """The ``results.csv`` row of one run, from its ``outputs/*.json`` and history.

    Parameters
    ----------
    run : pathlib.Path
        Run directory.

    Returns
    -------
    dict
        See the "results.csv columns" table in the documentation.
    """
    from rbc import resources

    p = rc.parse_run_path(run)
    o = run / "outputs"
    nu, pr, sp = _read(o / "nusselt.json"), _read(o / "profiles.json"), _read(o / "spectra.json")
    r = dict(cluster=CLUSTER, path=str(run), dim=p["dim"], log_ra=p["log_ra"], log_pr=p["log_pr"],
             label=p["label"], nz=p["nz"], nx=sp.get("Nx") or p["nx"])
    if nu:
        r.update(t0=nu["t0"], t1=nu["t1"], T_avg=nu["T_avg"], n_samples=nu["n_samples"],
                 Nu=nu["Nu"], Nu_se=nu["Nu_se"], Nu_gradT=nu["avg_grad_T_sq"]["Nu"],
                 Nu_vort=nu.get("avg_vorticity_sq", {}).get("Nu", np.nan),
                 vort_diff_pct=nu.get("vort_diff_pct", np.nan),
                 section_spread_pct=max(nu[k]["spread_pct"] for k in NU_KEYS),
                 qoi_spread_pct=nu["qoi_spread_pct"], converged=nu["converged"])
        if "se_tau" in nu["avg_wT"]:
            w = nu["avg_wT"]
            r.update(Nu_se_tau=w["se_tau"], tau_int=w["tau_int"], n_eff=w["n_eff"],
                     p_sections=min(nu[k]["p_sections"] for k in NU_KEYS),
                     drift_z=max((nu[k]["drift_z"] for k in NU_KEYS), key=abs),
                     t0_mser=nu.get("t0_mser"), converged_stat=nu.get("converged_stat"))
        if "plates" in nu:
            r.update(Nu_plates=0.5 * (nu["plates"]["Nu_bot"] + nu["plates"]["Nu_top"]),
                     plate_imbalance_z=nu["plates"]["imbalance_z"])
        if "Re" in nu:
            r.update(Re=nu["Re"]["Re"], Re_se=nu["Re"]["se"])
    for k in ("delta_T", "delta_T_rms", "delta_u", "delta_u_slope", "n_delta_T", "n_delta_T_rms",
              "n_delta_u", "n_delta_u_slope", "Nu_flux_spread_pct", "Re_rms"):
        if k in pr:
            r[k] = pr[k]
    r.update({k: v for k, v in sp.items() if k.startswith(("D_", "Dmin_"))})
    hist = resources.summary(run)
    for k in ("nx", "ny", "nz"):
        if hist.get(k):
            r[k] = hist.pop(k)
        hist.pop(k, None)
    r.update(hist)
    return r


def mark_production(df: pd.DataFrame, cases: pd.DataFrame | None) -> pd.DataFrame:
    """Flag one production run per (dim, Ra, Pr).

    The run with the ``nz`` pinned in ``cases.csv`` wins; else the highest-``nz`` run
    marked converged; else the highest-``nz`` run.

    Parameters
    ----------
    df : pandas.DataFrame
        All rows.
    cases : pandas.DataFrame or None
        ``cases.csv`` as loaded by :func:`rbc.cases.load_cases`.

    Returns
    -------
    pandas.DataFrame
        Copy of ``df`` with a boolean ``production`` column.
    """
    df = df.copy()
    df["production"] = False
    pinned = {}
    if cases is not None and "nz" in cases:
        for _, c in cases.dropna(subset=["nz"]).iterrows():
            pinned[(c["dim"], c["log_ra"], c["log_pr"])] = int(c["nz"])
    for key, g in df.groupby(["dim", "log_ra", "log_pr"]):
        if key in pinned and (g["nz"] == pinned[key]).any():
            idx = g.index[g["nz"] == pinned[key]][0]
        else:
            ok = g[g.get("converged", pd.Series(False, index=g.index)).fillna(False).astype(bool)]
            idx = (ok if len(ok) else g)["nz"].idxmax()
        df.loc[idx, "production"] = True
    return df


def report(df: pd.DataFrame, cases: pd.DataFrame) -> None:
    """Print disagreements between ``cases.csv`` statuses and the collected results.

    Parameters
    ----------
    df : pandas.DataFrame
        Collected rows with the ``production`` flag.
    cases : pandas.DataFrame
        The registry.
    """
    prod = df[df["production"]]
    have = set(zip(prod["dim"], prod["log_ra"], prod["log_pr"]))
    for _, c in cases[cases["dim"].isin(df["dim"].unique())].iterrows():
        key = (c["dim"], c["log_ra"], c["log_pr"])
        if c["status"] == "finished" and key not in have:
            print(f"  cases.csv says finished but no results: {key}")
        if c["status"] != "finished" and key in have:
            r = prod[(prod["dim"] == key[0]) & (prod["log_ra"] == key[1]) & (prod["log_pr"] == key[2])].iloc[0]
            if bool(r.get("converged", False)):
                print(f"  converged results but cases.csv says {c['status']}: {key}")


def write_tables(df: pd.DataFrame, out: Path) -> list[Path]:
    """Write ``results.csv`` and, if pyarrow is importable, ``results.parquet``.

    Parameters
    ----------
    df : pandas.DataFrame
        The table.
    out : pathlib.Path
        Path of the CSV; the Parquet file goes next to it.

    Returns
    -------
    list of pathlib.Path
        Files written.
    """
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False, float_format="%.10g")
    written = [out]
    try:
        import pyarrow  # noqa: F401
    except ImportError:
        return written
    df.to_parquet(out.with_suffix(".parquet"), index=False)
    written.append(out.with_suffix(".parquet"))
    return written


def _t0_arg(v: str):
    """argparse type: a number or ``"auto"``."""
    return v if v == "auto" else float(v)


def get_parser() -> argparse.ArgumentParser:
    """Command-line options (see the module docstring).

    Returns
    -------
    argparse.ArgumentParser
    """
    ap = argparse.ArgumentParser(prog="python -m rbc.collect",
                                 description="Collect run results into results.csv.")
    ap.add_argument("root", type=Path, nargs="?", default=None,
                    help="top of the run tree (default: $DATA_ROOT or $SCRATCH/ollie_rb_data/3d)")
    ap.add_argument("--recalculate", action="store_true", help="redo every run's analysis")
    ap.add_argument("--t0", type=_t0_arg, help="averaging start for every run (number or 'auto')")
    ap.add_argument("--jobs", type=int, default=min(32, os.cpu_count() or 1), help="parallel processes")
    ap.add_argument("--no-plots", action="store_true", help="skip the per-run diagnostic PDFs")
    ap.add_argument("--exclude", default=r"(_old|test)", help="regex of run paths to skip")
    ap.add_argument("--merge", type=Path, nargs="+", default=[],
                    help="results.csv files from other clusters to merge in (their rows are kept)")
    ap.add_argument("--out", type=Path, default=rc.RESULTS_CSV, help="output CSV")
    return ap


def main(argv=None) -> None:
    """Collect every run under the root and write the results tables.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments (default: ``sys.argv[1:]``).
    """
    a = get_parser().parse_args(argv)
    root = (a.root or default_root()).resolve()
    runs = find_runs(root, a.exclude)
    print(f"{len(runs)} runs under {root} (cluster {CLUSTER}), {a.jobs} workers")
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    rows = []
    with ProcessPoolExecutor(max_workers=max(1, a.jobs)) as ex:
        futs = [ex.submit(process_run, r, a.recalculate, a.t0, not a.no_plots) for r in runs]
        for f in as_completed(futs):
            r, msg = f.result()
            print(" ", msg)
            if r is not None:
                rows.append(r)
    new = pd.DataFrame(rows)
    old = [pd.read_csv(p) for p in ([a.out] if a.out.exists() else []) + list(a.merge)]
    if old:
        old = pd.concat(old, ignore_index=True)
        old = old[old["cluster"] != CLUSTER].drop_duplicates(subset=["cluster", "path"], keep="last")
        new = pd.concat([old, new], ignore_index=True)
    cases = rc.load_cases() if rc.CASES_CSV.exists() else None
    new = mark_production(new, cases).sort_values(["dim", "log_ra", "log_pr", "nz"]).reset_index(drop=True)
    for p in write_tables(new, a.out):
        print(f"wrote {len(new)} rows to {p}")
    if cases is not None:
        report(new, cases)


if __name__ == "__main__":
    main()