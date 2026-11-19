"""Check runs without a resolution study against the calibrated thresholds.

For each run (anything with ``outputs/nusselt.json``) and each measure in
``calibration.json`` (:mod:`testing.calibrate`): the run's value, the threshold and
whether the run is beyond it (``pass``). ``margin`` is value / threshold for measures
where larger is better and threshold / value otherwise, so ``margin > 1`` passes.

A run outside the (Ra, Pr) range of the calibration is flagged: its thresholds are an
extrapolation. Writes ``check.{txt,csv,pdf}`` to ``--out``::

    python -m testing.check PATH [PATH ...] --calibration CAL.json [--threshold strict|safe|bad] [--out DIR]

A path may be a run directory or any folder above runs (every run below it is checked).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from rbc import cases as rc
from testing import convergence as cv


def find_runs(paths) -> list[Path]:
    """Run directories (with ``outputs/nusselt.json``) at or below the paths.

    Parameters
    ----------
    paths : list of path-like

    Returns
    -------
    list of pathlib.Path
    """
    runs = set()
    for p in map(Path, paths):
        for nj in [p / "outputs" / "nusselt.json", *p.rglob("outputs/nusselt.json")]:
            if nj.exists() and rc.parse_run_path(nj.parent.parent) is not None:
                runs.add(nj.parent.parent.resolve())
    return sorted(runs)


def check(runs, cal: dict, which: str = "safe") -> pd.DataFrame:
    """Every run against every calibrated threshold.

    Parameters
    ----------
    runs : list of pathlib.Path
        Run directories.
    cal : dict
        Contents of ``calibration.json``.
    which : {"strict", "safe", "bad"}, optional
        Which threshold (:mod:`testing.calibrate`).

    Returns
    -------
    pandas.DataFrame
        One row per run: :func:`testing.convergence.run_row` fields, ``in_range`` and per
        measure ``<m>_threshold``, ``<m>_margin``, ``<m>_pass``; ``n_pass`` and
        ``best_pass`` (the verdict of the best calibrated measure).
    """
    ra0, ra1 = cal["log_ra_range"]
    pr0, pr1 = cal["log_pr_range"]
    rows = []
    for run in runs:
        r = cv.run_row(run)
        if r is None:
            continue
        r["in_range"] = bool(ra0 <= r["log_ra"] <= ra1 and pr0 <= r["log_pr"] <= pr1)
        n = 0
        for m, t in cal["overall"].items():
            thr = t[f"t_{which}"]
            v = float(r.get(m, np.nan))
            if thr is None or not np.isfinite(v) or v == 0:
                r[f"{m}_threshold"], r[f"{m}_margin"], r[f"{m}_pass"] = thr, np.nan, None
                continue
            margin = v / thr if t["larger_better"] else thr / v
            r[f"{m}_threshold"], r[f"{m}_margin"], r[f"{m}_pass"] = thr, margin, bool(margin > 1)
            n += margin > 1
        r["n_pass"] = int(n)
        best = next((m for m in cal["overall"] if r.get(f"{m}_pass") is not None), None)
        r["best_measure"] = best
        r["best_pass"] = r.get(f"{best}_pass") if best else None
        rows.append(r)
    return pd.DataFrame(rows)


def write_txt(df: pd.DataFrame, cal: dict, path: Path, which: str) -> None:
    """Write ``check.txt``.

    Parameters
    ----------
    df : pandas.DataFrame
        Output of :func:`check`.
    cal : dict
        Calibration.
    path : pathlib.Path
        File to write.
    which : str
        Threshold used.
    """
    ms = [m for m in cal["overall"] if f"{m}_margin" in df]
    L = [f"Runs against the calibrated thresholds (t_{which}; tol {cal['tol']:g}% of Nu-1). "
         f"Calibrated on log Ra {cal['log_ra_range']}, log Pr {cal['log_pr_range']}.",
         "margin > 1 passes; measures in calibration order (best first). '*' = outside the calibrated range.", ""]
    L.append(f"  {'run':<34}{'Nu':>9}" + "".join(f"{m:>15}" for m in ms) + "   verdict (best measure)")
    for r in df.itertuples():
        name = "/".join(Path(r.path).parts[-3:]) + ("" if r.in_range else " *")
        cells = "".join(f"{cv._g(getattr(r, f'{m}_margin'), '.2f'):>15}" for m in ms)
        verdict = "-" if r.best_pass is None else ("pass" if r.best_pass else "FAIL")
        L.append(f"  {name:<34}{r.Nu:9.4f}{cells}   {verdict} ({r.best_measure})")
    L += ["", "Thresholds:"] + [f"  {m:<16}{'>' if cal['overall'][m]['larger_better'] else '<'} "
                                f"{cal['overall'][m][f't_{which}']}" for m in ms]
    path.write_text("\n".join(L) + "\n")


def plot(df: pd.DataFrame, cal: dict, path: Path, which: str) -> None:
    """``check.pdf``: each run's margin per measure (one row per run).

    Parameters
    ----------
    df : pandas.DataFrame
        Output of :func:`check`.
    cal : dict
        Calibration.
    path : pathlib.Path
        File to write.
    which : str
        Threshold used.
    """
    import matplotlib.pyplot as plt
    from rbc import style

    style.use("diag")
    ms = [m for m in cal["overall"] if f"{m}_margin" in df and np.isfinite(df[f"{m}_margin"].astype(float)).any()]
    names = ["/".join(Path(r.path).parts[-3:]) + ("" if r.in_range else " *") for r in df.itertuples()]
    y = np.arange(len(df))
    fig, axes = plt.subplots(1, len(ms), figsize=(3.2 * len(ms) + 2.5, 0.28 * len(df) + 1.8), sharey=True,
                             layout="constrained", squeeze=False)
    for ax, m in zip(axes[0], ms):
        v = df[f"{m}_margin"].astype(float).to_numpy()
        ok = v > 1
        ax.scatter(v[ok], y[ok], color="#4EB265", s=18, zorder=3)
        ax.scatter(v[~ok], y[~ok], color="#DC050C", marker="X", s=26, zorder=3)
        ax.axvline(1, color="k", lw=0.8)
        if np.nanmin(v) > 0:
            ax.set_xscale("log")
        ax.set_title(m, fontsize=9)
        ax.set_xlabel("margin (> 1 passes)", fontsize=8)
    axes[0][0].set_yticks(y, names, fontsize=7)
    axes[0][0].invert_yaxis()
    fig.suptitle(f"Runs against calibrated thresholds (t_{which}); * outside the calibrated range")
    fig.savefig(path)
    plt.close(fig)


def get_parser() -> argparse.ArgumentParser:
    """Command-line options (see the module docstring).

    Returns
    -------
    argparse.ArgumentParser
    """
    p = argparse.ArgumentParser(prog="python -m testing.check", description=__doc__.split("\n")[0])
    p.add_argument("paths", type=Path, nargs="+", help="run directories, or folders above them")
    p.add_argument("--calibration", type=Path, required=True, help="calibration.json from testing.calibrate")
    p.add_argument("--threshold", choices=("strict", "safe", "bad"), default="safe",
                   help="which calibrated threshold (default t_safe)")
    p.add_argument("--out", type=Path, default=Path("."), help="output folder")
    p.add_argument("--no-plot", action="store_true", help="skip check.pdf")
    return p


def main(argv=None) -> None:
    """Command-line entry point.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments (default: ``sys.argv[1:]``).
    """
    a = get_parser().parse_args(argv)
    cal = json.loads(a.calibration.read_text())
    df = check(find_runs(a.paths), cal, a.threshold)
    if df.empty:
        print("no runs with outputs/nusselt.json found")
        return
    a.out.mkdir(parents=True, exist_ok=True)
    df.to_csv(a.out / "check.csv", index=False)
    write_txt(df, cal, a.out / "check.txt", a.threshold)
    if not a.no_plot:
        plot(df, cal, a.out / "check.pdf", a.threshold)
    print((a.out / "check.txt").read_text())


if __name__ == "__main__":
    main()
