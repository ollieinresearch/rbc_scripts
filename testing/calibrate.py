"""Which resolution measure predicts convergence of Nu? Calibration over every multi-resolution point.

Every parameter folder under the given roots with two or more analysed runs goes through
the plateau test (:mod:`testing.convergence`), which labels each run ``good`` or ``bad``
(shown within or beyond ``tol`` of the converged Nu), or, when the noise is too large to
tell, ``likely good`` or ``likely bad`` by its point estimate. Then, for each resolution
measure (k_max η, k_max η_T, decades, decades of n²E, the vorticity identity), with the
best value of the measure among a set of runs as threshold:

``t_strict``
    over every run not shown to be good. Beyond it, every run was shown to be good.
``t_safe``
    over bad and likely bad runs (the default). Beyond it, no run was, or looked, beyond
    tol.
``t_bad``
    over bad runs only. Beyond it, no run was shown to be beyond tol.
``ok_below``
    Good and likely good runs that ``t_safe`` would still reject. A measure that predicts convergence
    well separates the runs with few of these, and gives similar thresholds at every
    (Ra, Pr) (``per_point``): then one threshold can be applied to runs without a
    resolution study (:mod:`testing.check`).

Writes ``calibration.{txt,csv,pdf,json}`` to ``--out`` (default ``<first root>/outputs/testing``)
and each folder's ``outputs/testing/convergence.*``::

    python -m testing.calibrate ROOT [ROOT ...] [--tol 1] [--alpha 0.05] [--out DIR]

A root may be the data root (``.../3d``), a Ra folder or one parameter folder.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from rbc import cases as rc
from testing import convergence as cv


def param_dirs(roots) -> list[Path]:
    """Parameter folders under the roots that hold at least one analysed run.

    Parameters
    ----------
    roots : list of path-like
        Data root, Ra folders, parameter folders or run directories.

    Returns
    -------
    list of pathlib.Path
        Sorted, without duplicates.
    """
    found = set()
    for root in map(Path, roots):
        for nj in [root / "outputs" / "nusselt.json", *root.rglob("outputs/nusselt.json")]:
            run = nj.parent.parent
            if nj.exists() and rc.parse_run_path(run) is not None:
                found.add(run.parent.resolve())
    return sorted(found)


def expense_key(log_ra: float, log_pr: float) -> float:
    """Ordering of parameter points by cost: log10 of N⁴ with N ∝ Ra^(1/4) Pr^(-1/2).

    Parameters
    ----------
    log_ra, log_pr : float

    Returns
    -------
    float
    """
    return log_ra - 2 * log_pr


def thresholds(df: pd.DataFrame, col: str, larger_better: bool) -> dict:
    """Empirical thresholds of one measure.

    Parameters
    ----------
    df : pandas.DataFrame
        Calibration rows with ``status`` and the measure.
    col : str
        Measure.
    larger_better : bool
        Direction of the measure.

    Returns
    -------
    dict
        ``larger_better``, ``n`` (runs with the measure), counts per status (``n_good``,
        ``n_likely_good``, ``n_likely_bad``, ``n_bad``), ``t_strict``, ``t_safe``,
        ``t_bad`` (None if there is no such run), and ``ok_above``, ``ok_below``: good or
        likely good runs beyond and not beyond ``t_safe``.
    """
    d = df[np.isfinite(df[col].astype(float))]
    sign = 1 if larger_better else -1
    v = sign * d[col].astype(float)
    st = d["status"]
    ok = st.isin(["good", "likely good"])

    def best(mask):
        return float(sign * v[mask].max()) if mask.any() else None

    t_safe = best(~ok)
    cut = sign * t_safe if t_safe is not None else -np.inf
    return dict(larger_better=larger_better, n=int(len(d)),
                **{f"n_{k.replace(' ', '_')}": int((st == k).sum()) for k in ("good", "likely good", "likely bad", "bad")},
                t_strict=best(st != "good"), t_safe=t_safe, t_bad=best(st == "bad"),
                ok_above=int((ok & (v > cut)).sum()), ok_below=int((ok & (v <= cut)).sum()))


def calibrate(df: pd.DataFrame) -> tuple[dict, dict]:
    """Thresholds of every measure, overall and per parameter point.

    Parameters
    ----------
    df : pandas.DataFrame
        Calibration rows.

    Returns
    -------
    overall : dict
        Measure -> :func:`thresholds`, for measures with data, best first (smallest
        fraction of good and likely good runs rejected by ``t_safe``).
    per_point : dict
        Measure -> {"(log Ra, log Pr)": (t_safe, t_bad)}.
    """
    overall, per_point = {}, {}
    for col, (_, hi, _) in cv.PREDICTORS.items():
        if col not in df or not np.isfinite(df[col].astype(float)).any():
            continue
        overall[col] = thresholds(df, col, hi)
        per_point[col] = {f"({a:g}, {b:g})": (t["t_safe"], t["t_bad"])
                          for (a, b), g in df.groupby(["log_ra", "log_pr"])
                          for t in [thresholds(g, col, hi)]}
    order = sorted(overall, key=lambda c: (overall[c]["ok_below"] / max(overall[c]["ok_below"] + overall[c]["ok_above"], 1),
                                           -overall[c]["n"]))
    return {c: overall[c] for c in order}, per_point


def _t(v, fmt=".3g") -> str:
    """A threshold for the report."""
    return "-" if v is None else format(v, fmt)


def write_txt(df: pd.DataFrame, summaries: dict, overall: dict, per_point: dict, path: Path, tol: float,
              alpha: float) -> None:
    """Write ``calibration.txt``.

    Parameters
    ----------
    df : pandas.DataFrame
        Calibration rows.
    summaries : dict
        "(log Ra, log Pr)" -> summary of :func:`testing.convergence.analyse`.
    overall, per_point : dict
        Outputs of :func:`calibrate`.
    path : pathlib.Path
        File to write.
    tol, alpha : float
        Tolerance (% of Nu-1) and plateau significance level.
    """
    L = [f"Resolution calibration: tol {tol:g}% of Nu-1, plateau alpha {alpha:g}", "",
         "Parameter points (log Ra, log Pr):"]
    for key, s in summaries.items():
        g = df[(df["key"] == key)]
        plat = (f"plateau from {s['plateau_start']} ({s['plateau_n']} runs, p = {s['p']:.2g})"
                if s["plateau_start"] else f"NO PLATEAU (two finest disagree, p = {s['p']:.2g})")
        L.append(f"  {key:<14} {len(g)} runs, {plat}, Nu = {s['Nu_ref']:.5f} ± {s['Nu_ref_se']:.2g}; "
                 + ", ".join(f"{r.label} {r.status}" for r in g.itertuples()))
    L += ["", "Thresholds, best measure first ('>' larger is better, '<' smaller is). Beyond t_strict every run",
          "was shown good; beyond t_safe none was or looked bad; beyond t_bad none was shown bad.",
          "ok_below: good or likely good runs that t_safe would still reject.", ""]
    L.append(f"  {'measure':<16}{'dir':>4}{'runs':>6}{'good':>6}{'l.good':>7}{'l.bad':>6}{'bad':>5}"
             f"{'t_strict':>10}{'t_safe':>9}{'t_bad':>9}{'ok_above':>10}{'ok_below':>10}")
    for c, t in overall.items():
        L.append(f"  {c:<16}{'>' if t['larger_better'] else '<':>4}{t['n']:>6}{t['n_good']:>6}"
                 f"{t['n_likely_good']:>7}{t['n_likely_bad']:>6}{t['n_bad']:>5}{_t(t['t_strict']):>10}"
                 f"{_t(t['t_safe']):>9}{_t(t['t_bad']):>9}{t['ok_above']:>10}{t['ok_below']:>10}")
    L += ["", "Per-point t_safe (t_bad) — a measure that predicts convergence gives similar values at every point:"]
    for c in overall:
        L.append(f"  {c:<16}" + "  ".join(f"{k} {_t(a)} ({_t(b)})" for k, (a, b) in per_point[c].items()))
    L += ["", "Runs:", f"  {'point':<14}{'label':<22}{'delta%':>8}{'sigma%':>8}{'bound%':>8}  {'status':<12}"
          + "".join(f"{c:>14}" for c in overall)]
    for r in df.itertuples():
        L.append(f"  {r.key:<14}{r.label:<22}{r.delta_pct:+8.2f}{r.sigma_pct:8.2f}{r.bound_pct:8.2f}  {r.status:<12}"
                 + "".join(f"{cv._g(getattr(r, c), '.3g'):>14}" for c in overall))
    path.write_text("\n".join(L) + "\n")


def plot(df: pd.DataFrame, overall: dict, path: Path, tol: float) -> None:
    """``calibration.pdf``: the difference from the converged Nu against each measure.

    Points are coloured by parameter point, red (cheapest) to purple (dearest); markers as in
    :data:`testing.convergence.STATUSES`. Dashed line: ``t_safe``; dotted: ``t_bad``.

    Parameters
    ----------
    df : pandas.DataFrame
        Calibration rows.
    overall : dict
        Output of :func:`calibrate`.
    path : pathlib.Path
        File to write.
    tol : float
        Tolerance, drawn as a band.
    """
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from rbc import style

    style.use("diag")
    keys = sorted(df["key"].unique(), key=lambda k: expense_key(*df.loc[df["key"] == k, ["log_ra", "log_pr"]].iloc[0]))
    col = dict(zip(keys, style.ranked_colours(len(keys))))
    cols = list(overall)
    nc = 3
    nr = int(np.ceil(len(cols) / nc))
    fig, axes = plt.subplots(nr, nc, figsize=(5.2 * nc, 4.2 * nr), layout="constrained", squeeze=False)
    for ax, c in zip(axes.flat, cols):
        label, hi, logx = cv.PREDICTORS[c]
        d = df[np.isfinite(df[c].astype(float))]
        ax.axhspan(-tol, tol, color="0.9", lw=0, zorder=0)
        ax.axhline(0, color="0.6", lw=0.8)
        for r in d.itertuples():
            m, filled = cv.STATUSES[r.status]
            ax.errorbar(getattr(r, c), r.delta_pct, yerr=2 * r.sigma_pct, fmt=m, ms=5, capsize=2,
                        ecolor=col[r.key], **style.marker_kw(col[r.key], filled))
        t = overall[c]
        if t["t_safe"] is not None:
            ax.axvline(t["t_safe"], color="k", ls="--", lw=0.9)
        if t["t_bad"] is not None:
            ax.axvline(t["t_bad"], color="k", ls=":", lw=0.9)
        if logx and (d[c].astype(float) > 0).all():
            ax.set_xscale("log")
        ax.set_xlabel(label + (" (larger is better)" if hi else " (smaller is better)"), fontsize=8)
        ax.set_ylabel("Nu − Nu_converged  [% of Nu−1] ± 2σ", fontsize=8)
        ax.set_title(f"{c}: t_safe {_t(t['t_safe'])}, ok rejected {t['ok_below']}/{t['ok_below'] + t['ok_above']}",
                     fontsize=9)
    for ax in list(axes.flat)[len(cols):]:
        ax.axis("off")
    handles = [Line2D([], [], color=col[k], marker="o", ls="") for k in keys]
    handles += [Line2D([], [], color="k", marker=m, ls="", mfc="k" if f else "none") for m, f in cv.STATUSES.values()]
    fig.legend(handles, [f"(log Ra, log Pr) = {k}" for k in keys] + list(cv.STATUSES), loc="outside right upper",
               fontsize=8)
    fig.suptitle(f"Resolution calibration (grey: ±{tol:g}% of Nu−1)")
    fig.savefig(path)
    plt.close(fig)


def run(roots, out: Path | None = None, tol: float = cv.TOL, alpha: float = cv.ALPHA,
        include_variants: bool = False, plots: bool = True) -> dict | None:
    """Run the plateau test on every multi-resolution point and calibrate the measures.

    Parameters
    ----------
    roots : list of path-like
        Where to look for runs.
    out : pathlib.Path, optional
        Output folder (default ``<first root>/outputs/testing``).
    tol, alpha : float, optional
        Tolerance (% of Nu-1) and plateau significance level.
    include_variants : bool, optional
        Keep timestepper and safety-factor variants.
    plots : bool, optional
        Write the PDFs.

    Returns
    -------
    dict or None
        The contents of ``calibration.json``; None if no point has two runs.
    """
    rows, summaries = [], {}
    for pdir in param_dirs(roots):
        r = cv.run(pdir, tol, alpha, include_variants, plots)
        if r is None:
            continue
        df, s = r
        key = f"({df['log_ra'].iloc[0]:g}, {df['log_pr'].iloc[0]:g})"
        summaries[key] = s
        rows.append(df.assign(key=key))
    if not rows:
        return None
    df = pd.concat(rows, ignore_index=True)
    overall, per_point = calibrate(df)
    out = Path(out) if out else Path(roots[0]).resolve() / "outputs" / "testing"
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / "calibration.csv", index=False)
    write_txt(df, summaries, overall, per_point, out / "calibration.txt", tol, alpha)
    if plots:
        plot(df, overall, out / "calibration.pdf", tol)
    rec = dict(tol=tol, alpha=alpha, overall=overall, per_point=per_point,
               log_ra_range=[float(df["log_ra"].min()), float(df["log_ra"].max())],
               log_pr_range=[float(df["log_pr"].min()), float(df["log_pr"].max())],
               points=list(summaries))
    (out / "calibration.json").write_text(json.dumps(rec, indent=1))
    return rec


def get_parser() -> argparse.ArgumentParser:
    """Command-line options (see the module docstring).

    Returns
    -------
    argparse.ArgumentParser
    """
    p = argparse.ArgumentParser(prog="python -m testing.calibrate", description=__doc__.split("\n")[0])
    p.add_argument("roots", type=Path, nargs="+", help="data root, Ra folders or parameter folders")
    p.add_argument("--out", type=Path, help="output folder (default <first root>/outputs/testing)")
    p.add_argument("--tol", type=float, default=cv.TOL, help="tolerance, percent of Nu-1")
    p.add_argument("--alpha", type=float, default=cv.ALPHA, help="significance level of the plateau test")
    p.add_argument("--include-variants", action="store_true", help="keep timestepper / safety-factor variants")
    p.add_argument("--no-plot", action="store_true", help="skip the PDFs")
    return p


def main(argv=None) -> None:
    """Command-line entry point.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments (default: ``sys.argv[1:]``).
    """
    a = get_parser().parse_args(argv)
    rec = run(a.roots, a.out, a.tol, a.alpha, a.include_variants, not a.no_plot)
    if rec is None:
        print("no parameter point with two or more analysed runs")
        return
    for c, t in rec["overall"].items():
        print(f"  {c:<16} t_safe {_t(t['t_safe'])}, ok runs rejected {t['ok_below']}/{t['ok_below'] + t['ok_above']}")


if __name__ == "__main__":
    main()
