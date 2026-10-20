"""Resolution studies (figure H). Every (Ra, Pr) with at least 3 resolutions in results.csv.

``study``
    per point, Nu from the three exact formulas relative to the finest run's <wT> value,
    and the spectral decades D, against N_z.
``decades``
    at one point, the runs whose D is nearest each target (default 4, 5, 6): (a) their
    time-averaged spectra, (b) the relative Nu error against D for all runs, with the rule
    "D decades -> D/2 significant figures" drawn as 10^(-D/2). Also prints the N_z needed
    for each target from a linear fit D(N_z), which holds for exponential (spectral)
    convergence.

D of a run = smallest of D[<E>] over {velocity, temperature} x {Chebyshev n_z,
horizontal k}. Usage::

    python -m rbc.plots.resolution study
    python -m rbc.plots.resolution decades --point "6 -2" --targets 4 5 6
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from rbc import cases as rc
from rbc import style
from rbc.spectra import compare

QOI = [("Nu", "o", r"$\langle wT\rangle$"), ("Nu_gradT", "D", r"$\langle|\nabla T|^2\rangle$"),
       ("Nu_vort", "*", r"$\langle|\omega|^2\rangle$")]
D_COLS = ["D_vel_nz", "D_vel_h", "D_temp_nz", "D_temp_h"]


def studies(df, min_runs=3):
    """Every (Ra, Pr) point with several resolutions.

    Parameters
    ----------
    df : pandas.DataFrame
        All results rows (not just production).
    min_runs : int, optional
        Fewest distinct N_z.

    Returns
    -------
    dict
        (log_ra, log_pr) -> rows sorted by N_z, with ``D_run`` added.
    """
    df = df.copy()
    df["D_run"] = df[[c for c in D_COLS if c in df]].min(axis=1)
    out = {}
    for key, g in df.groupby(["log_ra", "log_pr"]):
        if g["nz"].nunique() >= min_runs:
            out[key] = g.sort_values("nz")
    return out


def study_fig(st, out, png):
    """Figure ``resolution_study``: Nu error and decades against N_z, one column per point.

    Parameters
    ----------
    st : dict
        Output of :func:`studies`.
    out : pathlib.Path
        Output folder.
    png : bool
        Also write a PNG.
    """
    n = len(st)
    fig, axes = plt.subplots(2, n, figsize=style.figsize(1.0, 0.62 if n > 1 else 1.2), squeeze=False,
                             sharex="col")
    for j, ((lra, lpr), g) in enumerate(st.items()):
        ref = g.iloc[-1]["Nu"]
        ax = axes[0, j]
        ax.axhline(0, color="0.5", lw=0.6, ls="--")
        for col, mk, lab in QOI:
            ax.plot(g["nz"], 100 * (g[col] / ref - 1), marker=mk, ls="-", lw=0.7,
                    mfc="none" if mk == "D" else None, ms=4, label=lab)
        ax.set_title(rf"$Ra={style.pow10(lra)}$, $Pr={style.pow10(lpr)}$", fontsize=8)
        ax2 = axes[1, j]
        for c, ls in zip(D_COLS, ["-", "--", "-", "--"]):
            if c in g:
                fld = "vel" if "vel" in c else "temp"
                ax2.plot(g["nz"], g[c], ls=ls, marker=".", color="C0" if fld == "vel" else "C2",
                         label={"D_vel_nz": r"$u$, $n_z$", "D_vel_h": r"$u$, $k_h$",
                                "D_temp_nz": r"$b$, $n_z$", "D_temp_h": r"$b$, $k_h$"}[c])
        prod = g[g["production"].astype(bool)]
        for a_ in (ax, ax2):
            for nz in prod["nz"]:
                a_.axvline(nz, color="0.75", lw=0.8, zorder=0)
        ax2.set_xlabel(style.SYM["Nz"])
    axes[0, 0].set_ylabel(r"$Nu/Nu_{\rm ref}-1$ (\%)" if plt.rcParams["text.usetex"] else r"$Nu/Nu_{\rm ref}-1$ (%)")
    axes[1, 0].set_ylabel(r"decades $D$")
    axes[0, 0].legend(fontsize=6.5)
    axes[1, 0].legend(fontsize=6.5, ncol=2)
    print(style.save(fig, out / "resolution_study", png=png))


def decades_fig(g, targets, out, png, lra, lpr):
    """Figure ``resolution_decades``: spectra at target decades, and Nu error against D.

    Also prints the N_z needed for each target from a linear fit D(N_z), which holds
    for exponential (spectral) convergence.

    Parameters
    ----------
    g : pandas.DataFrame
        Runs of one point (an entry of :func:`studies`).
    targets : list of float
        Decades to show.
    out : pathlib.Path
        Output folder.
    png : bool
        Also write a PNG.
    lra, lpr : float
        log10 Ra and log10 Pr of the point.
    """
    D, nz = g["D_run"].values, g["nz"].values
    ref = g.iloc[-1]
    m, b = np.polyfit(nz, D, 1)
    print(f"Ra=1e{lra:g}, Pr=1e{lpr:g}: D ≈ {m:.4f} N_z {b:+.3f}")
    for t in targets:
        print(f"  D = {t}: N_z ≈ {(t - b) / m:.0f}" + ("" if D.min() <= t <= D.max() else "  (extrapolated)"))
    pick = []
    for t in targets:
        i = int(np.argmin(np.abs(D - t)))
        if abs(D[i] - t) > 0.5:
            print(f"  no run within 0.5 decades of D = {t} (nearest D = {D[i]:.2f}, N_z = {nz[i]})")
        pick.append(i)

    fig, (a1, a2) = plt.subplots(1, 2, figsize=style.figsize(1.0, 0.42))
    from matplotlib import colormaps
    cm = colormaps["viridis"]
    arrays = []
    for k, i in enumerate(dict.fromkeys(pick)):
        r = g.iloc[i]
        sp = compare.mean_spectra(Path(r["path"]), r.get("t0"))
        if sp is None:
            continue
        n, E, Dv = sp[("vel", "nz")]
        msk = (n >= 1) & (E > 0)
        a1.loglog(n[msk], E[msk], color=cm(0.85 * k / max(len(pick) - 1, 1)),
                  label=rf"$N_z={int(r['nz'])}$, $D={r['D_run']:.1f}$")
        arrays.append(E[msk])
    a1.set(xlabel=r"Chebyshev degree $n_z$", ylabel=r"$E_u(n_z)$")
    a1.legend(fontsize=6.5)

    others = g.iloc[:-1]
    for col, mk, lab in QOI:
        err = np.abs(others[col] / ref["Nu"] - 1)
        a2.semilogy(others["D_run"], err.clip(lower=1e-12), marker=mk, ls="none", ms=4,
                    mfc="none" if mk == "D" else None, label=lab)
    dd = np.linspace(min(D.min(), min(targets)) - 0.3, max(D.max(), max(targets)) + 0.3, 50)
    a2.semilogy(dd, 10 ** (-dd / 2), "k--", lw=0.8, label=r"$10^{-D/2}$")
    for i in pick:
        a2.axvline(D[i], color="0.8", lw=0.6, zorder=0)
    a2.set(xlabel=r"decades $D$", ylabel=r"$|Nu-Nu_{\rm ref}|/Nu_{\rm ref}$")
    a2.legend(fontsize=6.5)
    tag = f"ra{lra:g}_pr{lpr:g}"
    print(style.save(fig, out / f"resolution_decades_{tag}", png=png))


def get_parser() -> argparse.ArgumentParser:
    """Command-line options (see the module docstring).

    Returns
    -------
    argparse.ArgumentParser
    """
    p = argparse.ArgumentParser(prog="python -m rbc.plots.resolution", description=__doc__.split("\n")[0])
    p.add_argument("figure", choices=("study", "decades"), help="study: Nu and decades against N_z per point; decades: spectra at target decades")
    p.add_argument("--point", help='"logRa logPr" for decades')
    p.add_argument("--targets", type=float, nargs="+", default=[4, 5, 6], help="decades D to compare")
    p.add_argument("--min-runs", type=int, default=3, help="fewest resolutions for a point to count as a study")
    p.add_argument("--results", type=Path, default=rc.RESULTS_CSV, help="results table (default: data/results.csv)")
    p.add_argument("--out", type=Path, default=Path("figures"), help="output folder")
    p.add_argument("--png", action="store_true", help="also write PNG previews")
    p.add_argument("--no-tex", action="store_true", help="never use LaTeX for text (default: use it where it is installed)")
    return p


def main(argv=None) -> None:
    """Command-line entry point.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments (default: ``sys.argv[1:]``).
    """
    a = get_parser().parse_args(argv)
    style.use("paper", usetex=False if a.no_tex else None)
    df = rc.load_results(a.results, production=False)
    rc.require_columns(df, ["nz", "Nu_gradT", "Nu_vort"], "the resolution figure")
    st = studies(df, a.min_runs)
    if not st:
        raise SystemExit("no (Ra, Pr) point has enough resolutions")
    if a.figure == "study":
        study_fig(st, a.out, a.png)
    else:
        key = tuple(float(v) for v in a.point.split()) if a.point else next(iter(st))
        decades_fig(st[key], a.targets, a.out, a.png, *key)


if __name__ == "__main__":
    main()