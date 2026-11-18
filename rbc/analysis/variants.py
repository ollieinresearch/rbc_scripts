"""Compare the runs of one (Ra, Pr): resolution, timestepper and CFL safety factor.

Every run directory under one parameter folder (``.../ra1e7/pr1e-2/``) is a variant
of the same flow. This module pairs them up and turns each difference into a bound on the
systematic error of Nu. Writes ``PARAM_DIR/outputs/variants.txt``, ``variants.pdf`` and
``variants.json`` (read by :mod:`rbc.collect`). Usage::

    python -m rbc.analysis.variants PARAM_DIR [--tol 0.5]

Each run needs ``outputs/nusselt.json`` (and, for the decades panel,
``outputs/resolution.json``). The configuration of a run comes from its ``params.conf``
and its label (:func:`rbc.cases.run_config`).

Resolution pairs
    Runs with the same timestepper and safety factor but different (N_z, N_x). Each run
    is compared with the finest of its group; the finest with the next coarser one.
Timestep pairs
    Runs with the same (N_z, N_x) but a different timestepper or safety factor, each
    compared with the baseline (the run whose label has no suffix).

For a pair, with Q = Nu - 1 from ⟨wT⟩ (⟨|∇T|²⟩ is reported alongside)::

    delta = Q_b - Q_a,   sigma = sqrt(se_a² + se_b²),   bound = |delta| + 2 sigma

``bound`` is a 95% upper bound on |delta|, used as the bias bound of both runs. For the
finest run this is conservative: with exponential convergence its own error is a small
fraction of the difference to a coarser run. A pair ``agrees`` if |delta| ≤ max(2 sigma,
tol (Nu - 1)). Per run: ``B_res`` (from its resolution pair), ``B_dt`` (largest bound of
its timestep pairs) and ``B_sys = B_res + B_dt``, added linearly because they are biases,
not random errors. All bounds are in Nu units (the same for Nu and Nu - 1).

Inherited bounds
    A timestep variant without a resolution pair takes ``B_res`` of its baseline (same
    grid). A run without a timestep pair takes the largest ``B_dt`` measured with the
    same stepper and safety factor on the same or a coarser grid: with a CFL-limited
    step, dt shrinks as the grid is refined, so the time error does not grow. Inherited
    values are marked ``inherited from <label>``.

:func:`fill_bounds` estimates a bound for runs without a measured one, for
:mod:`rbc.collect`.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from rbc import cases as rc

#: Identities compared: key in nusselt.json -> label.
QOIS = {"avg_wT": "<wT>", "avg_grad_T_sq": "<|grad T|^2>"}
#: Default agreement tolerance, in percent of Nu - 1.
TOL_PCT = 0.5
#: Colours of the pair kinds in variants.pdf.
PAIR_COLOURS = {"resolution": "#4477AA", "timestep": "#EE7733"}
#: Implicit solves per step, for the cost estimate (Dedalus IMEX schemes).
STAGES = {"RK111": 1, "RK222": 2, "RK443": 4, "RKSMR": 3, "SBDF1": 1, "SBDF2": 1, "SBDF3": 1, "SBDF4": 1,
          "CNAB1": 1, "CNAB2": 1, "MCNAB2": 1, "CNLF2": 1}


def cost(nz, nx, ny, stepper: str, safety) -> float:
    """Relative cost of a run per unit simulated time.

    Grid points times solves per step divided by the CFL safety factor (the step is
    proportional to it). Only used to order runs and colour them by expense.

    Parameters
    ----------
    nz, nx, ny : int
        Resolution.
    stepper : str
        Timestepper name.
    safety : float or None
        CFL safety factor (None: taken as 1).

    Returns
    -------
    float
    """
    if not (nz and nx and ny):
        return np.nan
    return float(nz * nx * ny) * STAGES.get(str(stepper).upper(), 1) / (safety or 1.0)


def _read(p: Path) -> dict:
    """A JSON file as a dict ({} if missing)."""
    return json.loads(p.read_text()) if p.exists() else {}


def load_runs(param_dir: Path) -> pd.DataFrame:
    """Every analysed run of one parameter folder.

    Parameters
    ----------
    param_dir : pathlib.Path
        e.g. ``.../ra1e7/pr1e-2``.

    Returns
    -------
    pandas.DataFrame
        One row per run with ``label``, ``path``, ``nz``, ``nx``, ``stepper``,
        ``safety``, ``variant``, ``grid`` (N_z N_x N_y), ``cost`` (:func:`cost`), ``mismatch``, Nu and standard
        error per identity (``Nu_<qoi>``, ``se_<qoi>``) and ``D_run`` when the run has
        ``resolution.json``. Sorted by cost.
    """
    rows = []
    for run in sorted(p for p in Path(param_dir).iterdir() if p.is_dir()):
        nu = _read(run / "outputs" / "nusselt.json")
        if not nu or rc.parse_run_path(run) is None:
            continue
        res = _read(run / "outputs" / "resolution.json")
        cfg = rc.run_config(run)
        nz = res.get("Nz") or cfg["nz"]
        nx = res.get("Nx") or cfg["nx"] or (2 * nz if nz else None)
        ny = res.get("Ny") or nx
        r = dict(label=run.name, path=str(run), nz=nz, nx=nx, ny=ny, stepper=cfg["stepper"],
                 safety=cfg["safety"], variant=cfg["variant"], mismatch="; ".join(cfg["mismatch"]),
                 grid=float(nz * nx * ny) if nz and nx else np.nan,
                 cost=cost(nz, nx, ny, cfg["stepper"], cfg["safety"]), D_run=res.get("D_run", np.nan),
                 D_run_at=res.get("D_run_at", ""))
        for q in QOIS:
            r[f"Nu_{q}"] = nu.get(q, {}).get("Nu", np.nan)
            r[f"se_{q}"] = nu.get(q, {}).get("se", np.nan)
        for k in ("D_vel_x", "D_vel_y", "D_vel_z", "D_temp_x", "D_temp_y", "D_temp_z"):
            r[k] = res.get(k, np.nan)
        rows.append(r)
    df = pd.DataFrame(rows)
    return df.sort_values(["cost", "label"]).reset_index(drop=True) if len(df) else df


def compare(a: pd.Series, b: pd.Series, kind: str, tol_pct: float) -> dict:
    """Difference of Nu between two runs, per identity.

    Parameters
    ----------
    a, b : pandas.Series
        Rows of :func:`load_runs`; ``b`` is the finer (or the variant) run.
    kind : {"resolution", "timestep"}
        Which pair this is.
    tol_pct : float
        Agreement tolerance in percent of Nu - 1.

    Returns
    -------
    dict
        ``a``, ``b``, ``kind`` and, per identity q, ``delta_<q>``, ``sigma_<q>``,
        ``z_<q>``, ``pct_<q>`` (delta in % of Nu_b - 1), ``bound_<q>`` and
        ``agree_<q>``.
    """
    out = dict(a=a["label"], b=b["label"], kind=kind)
    for q in QOIS:
        d = b[f"Nu_{q}"] - a[f"Nu_{q}"]
        s = float(np.hypot(a[f"se_{q}"], b[f"se_{q}"]))
        q1 = b[f"Nu_{q}"] - 1
        out.update({f"delta_{q}": d, f"sigma_{q}": s, f"z_{q}": d / s if s > 0 else np.nan,
                    f"pct_{q}": 100 * d / q1, f"bound_{q}": abs(d) + 2 * s,
                    f"agree_{q}": bool(abs(d) <= max(2 * s, tol_pct / 100 * q1))})
    return out


def analyse(df: pd.DataFrame, tol_pct: float = TOL_PCT) -> tuple[pd.DataFrame, list[dict]]:
    """Pair the runs and assign bounds.

    Parameters
    ----------
    df : pandas.DataFrame
        Output of :func:`load_runs`.
    tol_pct : float, optional
        Agreement tolerance in percent of Nu - 1.

    Returns
    -------
    runs : pandas.DataFrame
        ``df`` with ``B_res``, ``B_res_vs``, ``B_dt``, ``B_dt_vs`` and ``B_sys`` (NaN
        where neither measured nor inherited). ``*_vs`` names the other run of the pair,
        or ``inherited from <label>``.
    pairs : list of dict
        Every comparison (:func:`compare`).
    """
    df = df.copy()
    df["B_res"], df["B_dt"] = np.nan, np.nan
    df["B_res_vs"], df["B_dt_vs"] = "", ""
    pairs = []
    key_dt = df["stepper"].astype(str) + "|" + df["safety"].astype(str)
    for _, g in df.groupby(key_dt, sort=False):
        g = g[np.isfinite(g["grid"])].sort_values("grid")
        if g["grid"].nunique() < 2:
            continue
        fine = g.iloc[-1]
        for i in g.index[:-1]:
            if g.loc[i, "grid"] == fine["grid"]:
                continue
            p = compare(g.loc[i], fine, "resolution", tol_pct)
            pairs.append(p)
            df.loc[i, ["B_res", "B_res_vs"]] = [p["bound_avg_wT"], fine["label"]]
        coarser = g[g["grid"] < fine["grid"]].iloc[-1]
        p = next(q for q in pairs if q["a"] == coarser["label"] and q["b"] == fine["label"])
        df.loc[fine.name, ["B_res", "B_res_vs"]] = [p["bound_avg_wT"], coarser["label"]]
    for _, g in df.groupby(["nz", "nx"], sort=False):
        base = g[g["variant"] == ""]
        if len(g) < 2 or base.empty:
            continue
        b0 = base.iloc[0]
        for i in g.index:
            if i == b0.name or (g.loc[i, "stepper"] == b0["stepper"] and g.loc[i, "safety"] == b0["safety"]):
                continue
            p = compare(b0, g.loc[i], "timestep", tol_pct)
            pairs.append(p)
            for j, other in ((i, b0["label"]), (b0.name, g.loc[i, "label"])):
                if not (df.loc[j, "B_dt"] >= p["bound_avg_wT"]):
                    df.loc[j, ["B_dt", "B_dt_vs"]] = [p["bound_avg_wT"], other]
            if np.isfinite(df.loc[b0.name, "B_res"]) and not np.isfinite(df.loc[i, "B_res"]):
                df.loc[i, ["B_res", "B_res_vs"]] = [df.loc[b0.name, "B_res"], f"inherited from {b0['label']}"]
    measured = df[df["B_dt_vs"].astype(bool)]
    for i in df.index[~df["B_dt_vs"].astype(bool)]:
        src = measured[(measured["grid"] <= df.loc[i, "grid"]) & (key_dt[measured.index] == key_dt[i])]
        if len(src):
            k = src["B_dt"].idxmax()
            df.loc[i, ["B_dt", "B_dt_vs"]] = [src.loc[k, "B_dt"], f"inherited from {src.loc[k, 'label']}"]
    df["B_sys"] = df[["B_res", "B_dt"]].sum(axis=1, min_count=1)
    return df, pairs


def _g(v, fmt=".3g") -> str:
    """A number for the text report, ``-`` if missing."""
    return format(v, fmt) if v is not None and np.isfinite(v) else "-"


def write_txt(df: pd.DataFrame, pairs: list[dict], path: Path, title: str, tol_pct: float) -> None:
    """Write ``variants.txt``.

    Parameters
    ----------
    df, pairs
        Outputs of :func:`analyse`.
    path : pathlib.Path
        File to write.
    title : str
        First line (the parameter point).
    tol_pct : float
        Agreement tolerance used.
    """
    L = [title, "", "Runs, cheapest first (Nu from <wT>, 1 standard error; cost relative to the cheapest):"]
    L.append(f"  {'label':<30}{'N_z':>5}{'N_x':>6}{'stepper':>9}{'safety':>9}{'cost':>7}{'Nu':>11}{'se':>9}"
             f"{'D_run':>7}  limit")
    c0 = df["cost"].min()
    for _, r in df.iterrows():
        L.append(f"  {r['label']:<30}{r['nz']:>5}{r['nx']:>6}{r['stepper']:>9}{_g(r['safety'], '.4g'):>9}"
                 f"{_g(r['cost'] / c0, '.3g'):>7}{r['Nu_avg_wT']:11.5f}{_g(r['se_avg_wT'], '.2g'):>9}"
                 f"{_g(r['D_run'], '.2f'):>7}  {r['D_run_at']}")
        if r["mismatch"]:
            L.append(f"    WARNING: {r['mismatch']}")
    L += ["", f"Pairs: delta = Nu(b) - Nu(a); bound = |delta| + 2 sigma; agree if |delta| <= max(2 sigma, "
              f"{tol_pct:g}% of Nu-1)"]
    for p in pairs:
        L.append(f"  {p['kind']:<11}{p['a']} -> {p['b']}")
        for q, name in QOIS.items():
            L.append(f"    {name:<14} delta {p[f'delta_{q}']:+.5f} ({p[f'pct_{q}']:+.3f}% of Nu-1), "
                     f"sigma {p[f'sigma_{q}']:.2g}, z {p[f'z_{q}']:+.1f}, bound {p[f'bound_{q}']:.3g}, "
                     f"{'agree' if p[f'agree_{q}'] else 'DIFFER'}")
    if not pairs:
        L.append("  none: no two runs differ only in resolution, or only in timestepper / safety factor")
    L += ["", "Systematic bounds per run (Nu units, '-' = not measured; B_sys = B_res + B_dt):"]
    for _, r in df.iterrows():
        pct = 100 * r["B_sys"] / (r["Nu_avg_wT"] - 1)
        L.append(f"  {r['label']:<30} B_res {_g(r['B_res'])} ({r['B_res_vs'] or '-'}), "
                 f"B_dt {_g(r['B_dt'])} ({r['B_dt_vs'] or '-'}), B_sys {_g(r['B_sys'])} = {_g(pct, '.3f')}% of Nu-1")
    path.write_text("\n".join(L) + "\n")


def plot(df: pd.DataFrame, pairs: list[dict], path: Path, title: str, tol_pct: float) -> None:
    """``variants.pdf``: Nu of every run, the pair differences and the decades.

    Runs are ordered and coloured by cost (red cheapest to purple dearest).

    Parameters
    ----------
    df, pairs
        Outputs of :func:`analyse`.
    path : pathlib.Path
        File to write.
    title : str
        Figure title.
    tol_pct : float
        Agreement tolerance, drawn as a band.
    """
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from rbc import style

    style.use("diag")
    colours = dict(zip(df["label"], style.ranked_colours(len(df))))
    marks = {"avg_wT": ("o", True), "avg_grad_T_sq": ("D", False)}
    fig = plt.figure(figsize=(15, 9.5), layout="constrained")
    gs = fig.add_gridspec(2, 2)
    ax = fig.add_subplot(gs[0, 0])
    x = np.arange(len(df))
    for i, r in df.iterrows():
        c = colours[r["label"]]
        if np.isfinite(r["B_sys"]):
            ax.add_patch(plt.Rectangle((i - 0.3, r["Nu_avg_wT"] - r["B_sys"]), 0.6, 2 * r["B_sys"],
                                       color=c, alpha=0.25, lw=0))
        for k, (q, (mk, filled)) in enumerate(marks.items()):
            ax.errorbar(i + (k - 0.5) * 0.24, r[f"Nu_{q}"], yerr=2 * r[f"se_{q}"], fmt=mk, ms=5, capsize=3,
                        ecolor=c, **style.marker_kw(c, filled))
    ax.set_xticks(x, df["label"], rotation=25, ha="right", fontsize=8)
    ax.set_xlim(-0.6, len(df) - 0.4)
    ax.set_ylabel("Nu")
    ax.set_title("Nu ± 2 SE; shaded: ± B_sys", fontsize=10)
    ax.legend([Line2D([], [], marker="o", color="k", ls=""), Line2D([], [], marker="D", mfc="none", color="k", ls="")],
              ["⟨wT⟩", "⟨|∇T|²⟩"], fontsize=8)

    ax = fig.add_subplot(gs[0, 1])
    res = df[np.isfinite(df["D_run"])]
    if len(res):
        for col, ls in (("x", "-"), ("y", "--"), ("z", ":")):
            for fld, c in (("vel", "k"), ("temp", "0.6")):
                ax.plot(res.index, res[f"D_{fld}_{col}"], ls=ls, marker=".", color=c, label=f"{fld} {col}")
        ax.axhline(4, color="0.4", lw=0.8, ls="--")
        ax.legend(fontsize=7, ncol=2)
    ax.set_xticks(x, df["label"], rotation=25, ha="right", fontsize=8)
    ax.set_xlim(-0.6, len(df) - 0.4)
    ax.set_ylabel("decades")
    ax.set_title("Decades of the time-averaged marginal spectra (screening only)", fontsize=10)

    ax = fig.add_subplot(gs[1, :])
    ax.axhspan(-tol_pct, tol_pct, color="0.88", lw=0, zorder=0)
    ax.axhline(0, color="0.5", lw=0.8)
    for k, p in enumerate(pairs):
        nb = df.loc[df["label"] == p["b"], "Nu_avg_wT"].iloc[0] - 1
        c = PAIR_COLOURS[p["kind"]]
        for j, (q, (mk, filled)) in enumerate(marks.items()):
            ax.errorbar(k + (j - 0.5) * 0.2, p[f"pct_{q}"], yerr=200 * p[f"sigma_{q}"] / nb, fmt=mk, ms=5,
                        capsize=3, color=c, ecolor=c, mfc=c if filled else "none")
    ax.set_xticks(range(len(pairs)), [f"{p['a']}\n→ {p['b']}" for p in pairs], fontsize=8)
    ax.set_xlim(-0.6, max(len(pairs), 1) - 0.4)
    ax.set_ylabel("Δ(Nu − 1)  [% of Nu − 1]")
    ax.set_title(f"Pair differences ± 2σ; grey: ±{tol_pct:g}% tolerance", fontsize=10)
    ax.legend([Line2D([], [], color=c, marker="o", ls="") for c in PAIR_COLOURS.values()],
              ["resolution pair", "timestepper / safety pair"], fontsize=8)
    fig.suptitle(title)
    fig.savefig(path)
    plt.close(fig)


def run(param_dir, tol_pct: float = TOL_PCT, plots: bool = True) -> dict | None:
    """Compare every run of one parameter folder and write the outputs.

    Parameters
    ----------
    param_dir : str or pathlib.Path
        e.g. ``.../ra1e7/pr1e-2`` (a run directory is accepted and its parent used).
    tol_pct : float, optional
        Agreement tolerance in percent of Nu - 1.
    plots : bool, optional
        Write ``variants.pdf``.

    Returns
    -------
    dict or None
        The contents of ``variants.json``; None with fewer than two analysed runs.
    """
    param_dir = Path(param_dir).resolve()
    if rc.parse_run_path(param_dir) is not None:
        param_dir = param_dir.parent
    df = load_runs(param_dir)
    if len(df) < 2:
        return None
    df, pairs = analyse(df, tol_pct)
    out = param_dir / "outputs"
    out.mkdir(exist_ok=True)
    title = "/".join(param_dir.parts[-2:])
    write_txt(df, pairs, out / "variants.txt", title, tol_pct)
    if plots:
        plot(df, pairs, out / "variants.pdf", title, tol_pct)
    keep = ["label", "nz", "nx", "stepper", "safety", "variant", "B_res", "B_res_vs", "B_dt", "B_dt_vs", "B_sys"]
    rec = dict(tol_pct=tol_pct, runs={r["label"]: {k: (None if isinstance(r[k], float) and not np.isfinite(r[k])
                                                       else r[k]) for k in keep}
                                       for _, r in df.iterrows()}, pairs=pairs)
    (out / "variants.json").write_text(json.dumps(rec, indent=1, default=lambda v: v.item()))
    return rec


def fill_bounds(df: pd.DataFrame) -> pd.DataFrame:
    """Systematic bound for every production row: measured, or estimated from the decades.

    A row without a measured ``B_sys`` gets ``B_sys_rel = max over references of
    B_ref_rel 10^((D_ref - D_run) / 2)``, the references being the production rows with
    a measured bound, and relative meaning divided by Nu - 1. This scales a measured
    bound by the "D decades -> 10^(-D/2)" rule of thumb, so it is only as good as that
    rule; the source is recorded so that figures and tables can tell the two apart.

    Parameters
    ----------
    df : pandas.DataFrame
        Results rows with ``Nu``, ``D_run``, ``B_sys`` (NaN if not measured) and
        ``production``.

    Returns
    -------
    pandas.DataFrame
        Copy with ``Nu_sys`` (bound in Nu units) and ``Nu_sys_source`` (``measured``,
        ``estimated`` or ``""``).
    """
    df = df.copy()
    df["Nu_sys"] = df.get("B_sys", pd.Series(np.nan, index=df.index))
    df["Nu_sys_source"] = np.where(np.isfinite(df["Nu_sys"]), "measured", "")
    if "D_run" not in df:
        return df
    prod = df["production"].astype(bool) if "production" in df else pd.Series(True, index=df.index)
    ref = df[prod & np.isfinite(df["Nu_sys"]) & np.isfinite(df["D_run"])]
    if ref.empty:
        return df
    rel = (ref["Nu_sys"] / (ref["Nu"] - 1)).to_numpy()
    Dref = ref["D_run"].to_numpy()
    need = prod & ~np.isfinite(df["Nu_sys"]) & np.isfinite(df["D_run"])
    for i in df.index[need]:
        b = np.max(rel * 10 ** ((Dref - df.loc[i, "D_run"]) / 2))
        df.loc[i, ["Nu_sys", "Nu_sys_source"]] = [b * (df.loc[i, "Nu"] - 1), "estimated"]
    return df


def get_parser() -> argparse.ArgumentParser:
    """Command-line options (see the module docstring).

    Returns
    -------
    argparse.ArgumentParser
    """
    p = argparse.ArgumentParser(prog="python -m rbc.analysis.variants", description=__doc__.split("\n")[0])
    p.add_argument("param_dir", type=Path, help="folder holding every run of one (Ra, Pr), or one of its runs")
    p.add_argument("--tol", type=float, default=TOL_PCT, help="agreement tolerance in percent of Nu-1")
    p.add_argument("--no-plot", action="store_true", help="skip variants.pdf")
    return p


def main(argv=None) -> None:
    """Command-line entry point.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments (default: ``sys.argv[1:]``).
    """
    a = get_parser().parse_args(argv)
    rec = run(a.param_dir, a.tol, not a.no_plot)
    if rec is None:
        print(f"{a.param_dir}: fewer than two analysed runs, nothing to compare")
        return
    for lab, r in rec["runs"].items():
        print(f"  {lab}: B_sys = {r['B_sys']}")


if __name__ == "__main__":
    main()
