"""Resolution convergence of Nu at one (Ra, Pr): the plateau test.

The runs of one parameter folder (``.../ra1e6/pr1e-2``) at different resolutions are
sorted by grid size. The **plateau** is the longest run of finest resolutions whose Nu
are consistent with one common value: for every suffix of the sorted runs (coarsest
dropped first), a χ² test of "all equal" against their standard errors; the plateau
starts at the coarsest run from which every suffix passes at level ``alpha``.

Each run is then compared with the plateau mean, leaving the run itself out when it is
part of the plateau (so a run is never compared with itself)::

    delta = Nu_run - Nu_ref,   sigma = sqrt(se_run² + se_ref²),   bound = |delta| + 2 sigma

all in percent of Nu - 1 (Nu from ⟨wT⟩). ``bound`` is an approximately 95% upper bound on
the run's resolution error. A run is

* ``good`` if bound <= tol: shown to be within tol of the converged value;
* ``bad`` if |delta| - 2 sigma > tol: shown to be further than tol from it;
* otherwise the noise is too large to tell, and the point estimate decides between
  ``likely good`` (|delta| <= tol) and ``likely bad``.

Without a plateau (the two finest runs disagree), the finest run is the reference for the
others and is itself compared with the next coarser run.

Runs that fail the statistical convergence rule of :mod:`rbc.analysis.nusselt`
(``converged_stat``) are shown and compared but never used for the plateau or the
reference, and are classed ``unconverged``: a short average cannot define the converged
value. Only baseline runs (no timestepper or safety-factor suffix in the label) with the
folder's most common aspect ratio (``Gam2_2``) are used, unless ``--include-variants``. Writes ``PARAM_DIR/outputs/testing/convergence.{txt,pdf,json}``::

    python -m testing.convergence PARAM_DIR [--tol 1] [--alpha 0.05]

Each run needs ``outputs/nusselt.json`` (:mod:`rbc.analysis.nusselt`); the decades come
from ``outputs/resolution.json`` (:mod:`rbc.analysis.resolution`) where it exists.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from rbc import cases as rc
from rbc.analysis.resolution import grid_spacings, scale_tests

#: Default tolerance, percent of Nu - 1.
TOL = 1.0
#: Default significance level of the plateau test.
ALPHA = 0.05
#: Run classes, best first, with their markers (filled: shown, open: point estimate only).
STATUSES = {"good": ("o", True), "likely good": ("o", False), "likely bad": ("s", False), "bad": ("X", True),
            "unconverged": ("d", False)}
#: Resolution measures: column -> (label, True if larger is better, log axis).
PREDICTORS = {
    "kmax_eta_min": ("k_max η (smallest direction)", True, True),
    "kmax_eta_x": ("k_max η (horizontal)", True, True),
    "kmax_eta_z": ("k_max η (vertical, mid-plane)", True, True),
    "kmax_etaT_min": ("k_max η_T (smallest direction)", True, True),
    "vort_abs_pct": ("|Nu_ω − Nu| (% of Nu−1)", False, True),
    "D_run": ("decades, smallest direction", True, False),
    "D_vel": ("decades, velocity", True, False),
    "D_temp": ("decades, temperature", True, False),
    "Dk2_vel": ("decades of n²E, velocity", True, False),
    "Dk2_temp": ("decades of n²E, temperature", True, False),
}


def _read(p: Path) -> dict:
    """A JSON file as a dict ({} if missing)."""
    return json.loads(p.read_text()) if p.exists() else {}


def _num(prm: dict, key: str, default: float) -> float:
    """A number from ``params.conf`` (``default`` if absent or not a number)."""
    try:
        return float(prm[key])
    except (KeyError, ValueError):
        return default


def run_row(run: Path) -> dict | None:
    """Nu, resolution and resolution measures of one run.

    Parameters
    ----------
    run : pathlib.Path
        Run directory with ``outputs/nusselt.json``.

    Returns
    -------
    dict or None
        ``label``, ``path``, ``log_ra``, ``log_pr``, ``Ra``, ``Pr``, ``nz``, ``nx``,
        ``ny``, ``grid``, ``stepper``, ``safety``, ``variant``, ``Nu``, ``se`` (⟨wT⟩),
        ``Nu_gradT``, ``Re``, ``T_avg``, ``converged``, ``vort_diff_pct``,
        ``vort_abs_pct``, the measures of :data:`PREDICTORS` (NaN where unknown) and the
        grid spacings; None without ``nusselt.json`` or a parsable path.
    """
    run = Path(run)
    p = rc.parse_run_path(run)
    nu = _read(run / "outputs" / "nusselt.json")
    if p is None or not nu:
        return None
    res = _read(run / "outputs" / "resolution.json")
    cfg = rc.run_config(run)
    prm = rc.read_params(run)
    nz = res.get("Nz") or cfg["nz"]
    nx = res.get("Nx") or cfg["nx"] or 2 * nz
    ny = res.get("Ny") or nx
    Lx, Ly = _num(prm, "LX", 2.0), _num(prm, "LY", 2.0)
    Ra, Pr = nu.get("Ra", 10 ** p["log_ra"]), nu.get("Pr", 10 ** p["log_pr"])
    sp = grid_spacings(nx, ny, nz, Lx, Ly, 1.0)
    sc = scale_tests(Ra, Pr, nu["Nu"], sp)
    vort = nu.get("vort_diff_pct", np.nan)
    r = dict(label=run.name, path=str(run), log_ra=p["log_ra"], log_pr=p["log_pr"], Ra=Ra, Pr=Pr,
             nz=nz, nx=nx, ny=ny, grid=float(nz * nx * ny), stepper=cfg["stepper"], safety=cfg["safety"],
             variant=cfg["variant"], gam=(re.search(r"Gam[^_]*_[^_]*", run.name) or [""])[0], Nu=nu["Nu"], se=nu.get("Nu_se", np.nan),
             Nu_gradT=nu.get("avg_grad_T_sq", {}).get("Nu", np.nan), Re=nu.get("Re", {}).get("Re", np.nan),
             T_avg=nu.get("T_avg", np.nan), converged=nu.get("converged_stat", nu.get("converged")),
             vort_diff_pct=vort, vort_abs_pct=abs(vort), eta=sc["eta"], eta_T=sc["eta_T"],
             dx=sp["x"], dy=sp["y"], dz=sp["z"],
             kmax_eta_x=min(sc["kmax_eta_x"], sc["kmax_eta_y"]), kmax_eta_z=sc["kmax_eta_z"],
             kmax_eta_min=min(sc[f"kmax_eta_{d}"] for d in "xyz"),
             kmax_etaT_min=min(sc[f"kmax_etaT_{d}"] for d in "xyz"),
             D_run=res.get("D_run", np.nan))
    for f in ("vel", "temp"):
        r[f"D_{f}"] = min((res.get(f"D_{f}_{d}", np.nan) for d in "xyz"), default=np.nan)
        r[f"Dk2_{f}"] = min((res.get(f"Dk2_{f}_{d}", np.nan) for d in "xyz"), default=np.nan)
    return r


def load_point(param_dir: Path, include_variants: bool = False) -> pd.DataFrame:
    """Every analysed run of one parameter folder, coarsest first.

    Parameters
    ----------
    param_dir : pathlib.Path
        e.g. ``.../ra1e6/pr1e-2``.
    include_variants : bool, optional
        Keep runs with a timestepper or safety-factor suffix.

    Returns
    -------
    pandas.DataFrame
        Rows of :func:`run_row`, sorted by grid size.
    """
    rows = [r for d in sorted(Path(param_dir).iterdir()) if d.is_dir() and (r := run_row(d))]
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    if not include_variants:
        df = df[df["variant"] == ""]
    if len(df) and df["gam"].nunique() > 1:
        gams = df["gam"].value_counts()
        keep = "Gam2_2" if "Gam2_2" in gams else gams.index[0]
        print(f"{param_dir}: using {keep} runs only, skipping " + ", ".join(df.loc[df["gam"] != keep, "label"]))
        df = df[df["gam"] == keep]
    return df.sort_values(["grid", "label"]).reset_index(drop=True)


def wmean(nu, se) -> tuple[float, float]:
    """Inverse-variance weighted mean and its standard error.

    Parameters
    ----------
    nu, se : array_like
        Values and standard errors.

    Returns
    -------
    tuple of float
    """
    w = 1 / np.asarray(se, float) ** 2
    return float(np.sum(w * nu) / np.sum(w)), float(1 / np.sqrt(np.sum(w)))


def chi2_p(nu, se) -> tuple[float, int, float]:
    """χ² test that values with standard errors ``se`` share one mean.

    Parameters
    ----------
    nu, se : array_like
        Values and standard errors.

    Returns
    -------
    chi2, dof, p : float, int, float
    """
    nu, se = np.asarray(nu, float), np.asarray(se, float)
    m, _ = wmean(nu, se)
    c = float(np.sum(((nu - m) / se) ** 2))
    dof = len(nu) - 1
    return c, dof, float(stats.chi2.sf(c, dof)) if dof > 0 else 1.0


def plateau(df: pd.DataFrame, alpha: float = ALPHA) -> int | None:
    """Index of the coarsest run of the plateau.

    Parameters
    ----------
    df : pandas.DataFrame
        Output of :func:`load_point` (coarsest first).
    alpha : float, optional
        Significance level: a suffix is consistent if its χ² p-value is at least alpha.

    Returns
    -------
    int or None
        Smallest k such that the runs k, k+1, ..., k', ... to the finest are consistent
        for every starting point from k on; None if the two finest disagree.
    """
    n = len(df)
    start = None
    for k in range(n - 2, -1, -1):
        sub = df.iloc[k:]
        if chi2_p(sub["Nu"], sub["se"])[2] < alpha:
            break
        start = k
    return start


def analyse(df: pd.DataFrame, tol: float = TOL, alpha: float = ALPHA) -> tuple[pd.DataFrame, dict]:
    """Compare every run with the converged value.

    Parameters
    ----------
    df : pandas.DataFrame
        Output of :func:`load_point`.
    tol : float, optional
        Tolerance in percent of Nu - 1.
    alpha : float, optional
        Significance level of the plateau test.

    Returns
    -------
    runs : pandas.DataFrame
        ``df`` with ``in_plateau``, ``ref`` (what the run was compared with), ``delta_pct``,
        ``sigma_pct``, ``z``, ``bound_pct`` and ``status`` (:data:`STATUSES`).
    summary : dict
        ``plateau_start`` (label or None), ``plateau_n``, ``chi2``, ``dof``, ``p``,
        ``Nu_ref`` and ``Nu_ref_se`` (plateau mean, or the finest run), ``tol``, ``alpha``.
    """
    df = df.copy()
    df["converged"] = [v is not False and v != False for v in df["converged"]]
    c = df[df["converged"]]
    if len(c) < 2:
        raise ValueError("fewer than two converged runs")
    k = plateau(c, alpha)
    df["in_plateau"] = False
    if k is not None:
        df.loc[c.index[k:], "in_plateau"] = True
        P = df[df["in_plateau"]]
        c2, dof, p = chi2_p(P["Nu"], P["se"])
        mu, se_mu = wmean(P["Nu"], P["se"])
    else:
        c2, dof, p = chi2_p(c["Nu"].iloc[-2:], c["se"].iloc[-2:])
        mu, se_mu = float(c["Nu"].iloc[-1]), float(c["se"].iloc[-1])
    out = []
    for i, r in df.iterrows():
        if k is not None and r["in_plateau"]:
            others = df[df["in_plateau"] & (df.index != i)]
            m, s = wmean(others["Nu"], others["se"])
            ref = "plateau without this run"
        elif k is not None:
            m, s, ref = mu, se_mu, "plateau"
        elif i == c.index[-1]:
            m, s, ref = float(c["Nu"].iloc[-2]), float(c["se"].iloc[-2]), c["label"].iloc[-2]
        else:
            m, s, ref = mu, se_mu, c["label"].iloc[-1]
        q1 = m - 1
        d = 100 * (r["Nu"] - m) / q1
        sig = 100 * np.hypot(r["se"], s) / q1
        b = abs(d) + 2 * sig
        if not r["converged"]:
            status = "unconverged"
        elif b <= tol:
            status = "good"
        elif abs(d) - 2 * sig > tol:
            status = "bad"
        else:
            status = "likely good" if abs(d) <= tol else "likely bad"
        out.append(dict(ref=ref, delta_pct=d, sigma_pct=sig, z=d / sig, bound_pct=b, status=status))
    df = pd.concat([df, pd.DataFrame(out, index=df.index)], axis=1)
    summary = dict(plateau_start=None if k is None else c["label"].iloc[k],
                   plateau_n=0 if k is None else int(df["in_plateau"].sum()), chi2=c2, dof=dof, p=p,
                   Nu_ref=mu, Nu_ref_se=se_mu, tol=tol, alpha=alpha)
    return df, summary


def _g(v, fmt=".3g") -> str:
    """A number for the text report, ``-`` if missing."""
    try:
        return format(v, fmt) if np.isfinite(v) else "-"
    except (TypeError, ValueError):
        return "-"


def write_txt(df: pd.DataFrame, s: dict, path: Path, title: str) -> None:
    """Write ``convergence.txt``.

    Parameters
    ----------
    df, s
        Outputs of :func:`analyse`.
    path : pathlib.Path
        File to write.
    title : str
        First line.
    """
    L = [title, ""]
    if s["plateau_start"] is None:
        L.append(f"No plateau: the two finest runs disagree (chi2 {s['chi2']:.2f}, p = {s['p']:.3g}). "
                 f"Reference: the finest run, Nu = {s['Nu_ref']:.5f} ± {s['Nu_ref_se']:.2g}.")
    else:
        L.append(f"Plateau from {s['plateau_start']} on ({s['plateau_n']} runs): chi2 = {s['chi2']:.2f} on "
                 f"{s['dof']} dof, p = {s['p']:.3g} (alpha {s['alpha']:g}).")
        L.append(f"Converged Nu (weighted plateau mean) = {s['Nu_ref']:.5f} ± {s['Nu_ref_se']:.2g} "
                 f"({100 * s['Nu_ref_se'] / (s['Nu_ref'] - 1):.3f}% of Nu-1).")
    L += ["", f"delta, sigma, bound in % of Nu-1; bound = |delta| + 2 sigma; good: bound <= {s['tol']:g}%, "
              f"bad: |delta| - 2 sigma > {s['tol']:g}%; otherwise likely good / likely bad by |delta|", ""]
    L.append(f"{'label':<22}{'Nu':>10}{'se%':>7}{'T_avg':>8}{'delta%':>9}{'sigma%':>8}{'z':>6}{'bound%':>8}"
             f"  {'status':<12}{'kη_x':>6}{'kη_z':>6}{'kηT':>7}{'Nu_ω%':>8}{'D_run':>6}{'Dk2_v':>7}  plateau")
    for _, r in df.iterrows():
        L.append(f"{r['label']:<22}{r['Nu']:10.5f}{_g(100 * r['se'] / (r['Nu'] - 1), '.2f'):>7}"
                 f"{_g(r['T_avg'], '.0f'):>8}{r['delta_pct']:+9.2f}{r['sigma_pct']:8.2f}{r['z']:+6.1f}"
                 f"{r['bound_pct']:8.2f}  {r['status']:<12}{_g(r['kmax_eta_x'], '.2f'):>6}"
                 f"{_g(r['kmax_eta_z'], '.2f'):>6}{_g(r['kmax_etaT_min'], '.1f'):>7}"
                 f"{_g(r['vort_diff_pct'], '+.1f'):>8}{_g(r['D_run'], '.2f'):>6}{_g(r['Dk2_vel'], '.2f'):>7}"
                 f"  {'yes' if r['in_plateau'] else ''}")
    L += ["", "kη: k_max η = π η / Δ (η from the exact dissipation, Δ the largest spacing); kηT: the same with "
              "the thermal scale; Nu_ω%: vorticity identity minus <wT>, % of Nu-1.",
          "Each run in the plateau is compared with the plateau mean of the other plateau runs."]
    path.write_text("\n".join(L) + "\n")


def plot(df: pd.DataFrame, s: dict, path: Path, title: str) -> None:
    """``convergence.pdf``: Nu against resolution, the differences and the vorticity identity.

    Runs are coloured by grid size, red (coarsest) to purple (finest).

    Parameters
    ----------
    df, s
        Outputs of :func:`analyse`.
    path : pathlib.Path
        File to write.
    title : str
        Figure title.
    """
    import matplotlib.pyplot as plt
    from rbc import style

    style.use("diag")
    col = style.ranked_colours(len(df))
    fig, axes = plt.subplots(1, 3, figsize=(16, 5), layout="constrained")
    q1 = s["Nu_ref"] - 1
    tol = s["tol"]

    ax = axes[0]
    x = np.arange(len(df))
    ax.axhspan(s["Nu_ref"] - tol / 100 * q1, s["Nu_ref"] + tol / 100 * q1, color="0.9", lw=0, zorder=0)
    ax.axhspan(s["Nu_ref"] - 2 * s["Nu_ref_se"], s["Nu_ref"] + 2 * s["Nu_ref_se"], color="0.7", lw=0, zorder=0)
    for i, r in df.iterrows():
        m, filled = STATUSES[r["status"]]
        ax.errorbar(x[i], r["Nu"], yerr=2 * r["se"], fmt=m, capsize=3, ecolor=col[i], **style.marker_kw(col[i], filled))
    ax.set_xticks(x, df["label"], rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("Nu (⟨wT⟩) ± 2 SE")
    ax.set_title(f"dark: plateau mean ± 2 SE; light: ± {tol:g}% of Nu−1", fontsize=9)

    ax = axes[1]
    ax.axhspan(-tol, tol, color="0.9", lw=0, zorder=0)
    ax.axhline(0, color="0.5", lw=0.8)
    for i, r in df.iterrows():
        m, filled = STATUSES[r["status"]]
        ax.errorbar(r["kmax_eta_min"], r["delta_pct"], yerr=2 * r["sigma_pct"], fmt=m, capsize=3,
                    ecolor=col[i], **style.marker_kw(col[i], filled))
    ax.set_xscale("log")
    ax.set_xlabel("k_max η (smallest direction)")
    ax.set_ylabel("Nu − Nu_ref  [% of Nu−1] ± 2σ")
    ax.set_title("filled o good, open o likely good, open □ likely bad, X bad", fontsize=9)

    ax = axes[2]
    v = df["vort_abs_pct"]
    for i, r in df.iterrows():
        if np.isfinite(r["vort_abs_pct"]) and r["vort_abs_pct"] > 0:
            m, filled = STATUSES[r["status"]]
            ax.plot(r["kmax_eta_min"], r["vort_abs_pct"], m, **style.marker_kw(col[i], filled))
    if np.isfinite(v).any():
        ax.set_xscale("log")
        ax.set_yscale("log")
    ax.axhline(tol, color="0.5", ls="--", lw=0.8)
    ax.set_xlabel("k_max η (smallest direction)")
    ax.set_ylabel("|Nu_ω − Nu|  [% of Nu−1]")
    ax.set_title("vorticity identity (reported, not a test)", fontsize=9)
    fig.suptitle(title)
    fig.savefig(path)
    plt.close(fig)


def run(param_dir, tol: float = TOL, alpha: float = ALPHA, include_variants: bool = False,
        plots: bool = True) -> tuple[pd.DataFrame, dict] | None:
    """Plateau test of one parameter folder; writes ``outputs/testing/convergence.*``.

    Parameters
    ----------
    param_dir : str or pathlib.Path
        e.g. ``.../ra1e6/pr1e-2``.
    tol : float, optional
        Tolerance in percent of Nu - 1.
    alpha : float, optional
        Significance level of the plateau test.
    include_variants : bool, optional
        Keep timestepper and safety-factor variants.
    plots : bool, optional
        Write ``convergence.pdf``.

    Returns
    -------
    tuple or None
        Outputs of :func:`analyse`; None with fewer than two converged runs (the point is
        then skipped, and ``convergence.txt`` says why).
    """
    param_dir = Path(param_dir).resolve()
    df = load_point(param_dir, include_variants)
    out = param_dir / "outputs" / "testing"
    title = "/".join(param_dir.parts[-2:])
    conv = [v is not False and v != False for v in df["converged"]] if len(df) else []
    if sum(conv) < 2:
        if len(df):
            out.mkdir(parents=True, exist_ok=True)
            (out / "convergence.txt").write_text(
                f"{title}\n\nSkipped: {sum(conv)} of {len(df)} runs pass the statistical convergence rule; "
                f"two are needed.\n" + "".join(f"  {l}: {'converged' if c else 'not converged'}\n"
                                                for l, c in zip(df["label"], conv)))
        return None
    df, s = analyse(df, tol, alpha)
    out.mkdir(parents=True, exist_ok=True)
    write_txt(df, s, out / "convergence.txt", title)
    if plots:
        plot(df, s, out / "convergence.pdf", title)
    rec = dict(summary=s, runs=json.loads(df.to_json(orient="records")))
    (out / "convergence.json").write_text(json.dumps(rec, indent=1))
    return df, s


def get_parser() -> argparse.ArgumentParser:
    """Command-line options (see the module docstring).

    Returns
    -------
    argparse.ArgumentParser
    """
    p = argparse.ArgumentParser(prog="python -m testing.convergence", description=__doc__.split("\n")[0])
    p.add_argument("param_dirs", type=Path, nargs="+", help="parameter folders (e.g. .../ra1e6/pr1e-2)")
    p.add_argument("--tol", type=float, default=TOL, help="tolerance, percent of Nu-1")
    p.add_argument("--alpha", type=float, default=ALPHA, help="significance level of the plateau test")
    p.add_argument("--include-variants", action="store_true", help="keep timestepper / safety-factor variants")
    p.add_argument("--no-plot", action="store_true", help="skip convergence.pdf")
    return p


def main(argv=None) -> None:
    """Command-line entry point.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments (default: ``sys.argv[1:]``).
    """
    a = get_parser().parse_args(argv)
    for d in a.param_dirs:
        r = run(d, a.tol, a.alpha, a.include_variants, not a.no_plot)
        if r is None:
            print(f"{d}: fewer than two converged runs, skipped")
            continue
        df, s = r
        print(f"{d}: plateau from {s['plateau_start']} (p = {s['p']:.3g}), Nu = {s['Nu_ref']:.5f} ± "
              f"{s['Nu_ref_se']:.2g}; " + ", ".join(f"{l} {st}" for l, st in zip(df["label"], df["status"])))


if __name__ == "__main__":
    main()
