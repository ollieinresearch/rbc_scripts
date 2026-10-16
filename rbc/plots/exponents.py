"""Local exponents alpha, beta in Q ~ Pr^alpha Ra^beta, estimated from the data
(see rbc.analysis.exponents), plotted against Ra and Pr.

Panels: alpha vs Ra | alpha vs Pr  /  beta vs Ra | beta vs Pr. Colour is the other
parameter. A weighted line exponent = c0 + c1 log10(param) is fitted in each panel;
c1 ~ 0 with chi2_nu(flat) ~ 1 means the exponent is constant there. Dotted lines mark
the reference values (default 1/2 and 1/3). Usage::

    python -m rbc.plots.exponents                       # Q = Nu - 1, two-point secants + cells
    python -m rbc.plots.exponents --method lines        # one regression per row/column
    python -m rbc.plots.exponents --surface 2           # overlay a quadratic surface fit
    python -m rbc.plots.exponents --q Nu --refs 0.5 0.3333 0.25

Writes ``figures/exponents_<q>_<method>.pdf``, a CSV of every estimate and a CSV of the
cell-centre estimates, and prints the surface-fit table, the single global power law,
the flatness fits and the slope along the slices at the cell centres.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import colormaps
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

from rbc import cases as rc
from rbc import style
from rbc.analysis import exponents as ex
from rbc.plots.common import decade_ticks

SYMB = {"alpha": r"\alpha", "beta": r"\beta"}
COORD = {"x": "Ra", "y": "Pr"}


def frac(v: float) -> str:
    """A number as a simple fraction (0.3333 -> ``1/3``), for the reference labels.

    Parameters
    ----------
    v : float

    Returns
    -------
    str
    """
    f = Fraction(v).limit_denominator(12)
    return f"{f.numerator}/{f.denominator}" if f.denominator > 1 else str(f.numerator)


def gather(d, method, use_cells):
    """All exponent estimates in one long table.

    Parameters
    ----------
    d : pandas.DataFrame
        Output of :func:`rbc.analysis.exponents.prepare`.
    method : {"pairs", "lines"}
        Two-point secants or row/column regressions.
    use_cells : bool
        Add the cell-centre estimates.

    Returns
    -------
    pandas.DataFrame
        Columns ``kind``, ``x``, ``y``, ``value``, ``err``, ``source`` (and ``lo``,
        ``hi`` for lines).
    """
    parts = [ex.edges(d) if method == "pairs" else ex.lines(d)]
    if use_cells:
        c = ex.cells(d)
        if len(c):
            for k in ("alpha", "beta"):
                parts.append(pd.DataFrame(dict(kind=k, x=c.x, y=c.y, value=c[k], err=c[f"{k}_err"],
                                               source="cell")))
    return pd.concat(parts, ignore_index=True)


def panel(ax, est, kind, coord, cmap, norm, refs, fit_source, surf=None, grid=None):
    """One exponent against one parameter, coloured by the other.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes.
    est : pandas.DataFrame
        Output of :func:`gather`.
    kind : {"alpha", "beta"}
        Exponent.
    coord : {"x", "y"}
        Against log10 Ra or log10 Pr.
    cmap : matplotlib.colors.Colormap
        Colour map for the other parameter.
    norm : matplotlib.colors.Normalize
        Its normalization.
    refs : list of float
        Reference values drawn as dotted lines.
    fit_source : {"edge", "line", "cell"}
        Estimates used for the flatness fit.
    surf : rbc.analysis.exponents.Surface, optional
        Surface fit to overlay along each row or column (within the data range).
    grid : pandas.DataFrame, optional
        The data points (needed with ``surf``).

    Returns
    -------
    dict
        Output of :func:`rbc.analysis.exponents.flatness` (empty if too few points).
    """
    other = "y" if coord == "x" else "x"
    e = est[est.kind == kind]
    for src, mk, ms, fill in (("edge", "o", 3.5, True), ("line", "o", 4, True), ("cell", "D", 2.6, False)):
        g = e[e.source == src]
        if g.empty:
            continue
        cols = cmap(norm(g[other]))
        xerr = None
        if src == "line" and "lo" in g and coord == ("y" if kind == "alpha" else "x"):
            xerr = np.vstack([10 ** g[coord] - 10 ** g["lo"], 10 ** g["hi"] - 10 ** g[coord]])
        for i, (_, r) in enumerate(g.iterrows()):
            ax.errorbar(10 ** r[coord], r["value"], yerr=r["err"] if np.isfinite(r["err"]) else None,
                        xerr=None if xerr is None else xerr[:, [i]], fmt=mk, ms=ms, color=cols[i],
                        mfc=cols[i] if fill else "white", mew=0.7, elinewidth=0.6, capsize=1.2, zorder=3)
    for v in refs:
        ax.axhline(v, color="0.45", lw=0.7, ls=":", zorder=1)
        ax.annotate(frac(v), (1.0, v), xycoords=("axes fraction", "data"), xytext=(2, 0),
                    textcoords="offset points", va="center", fontsize=7, color="0.35",
                    annotation_clip=False)
    fs = e[e.source == fit_source]
    stats = ex.flatness(fs[coord], fs["value"], fs["err"]) if len(fs) > 2 else {}
    if stats:
        t = np.linspace(fs[coord].min(), fs[coord].max(), 50)
        ax.plot(10 ** t, stats["mean"] + 0 * t, color="0.2", lw=0.6, ls="--", zorder=2)
        ax.plot(10 ** t, stats["slope"] * t + stats["intercept"], color="k", lw=1, zorder=2)
        chi = (rf"$\chi^2_\nu = {stats['chi2_flat']:.2g}$" if np.isfinite(stats["chi2_flat"])
               else "errors unknown")
        txt = (rf"$d{SYMB[kind]}/d\log_{{10}}{COORD[coord]} = {stats['slope']:+.3f}\pm{stats['slope_err']:.3f}$" "\n"
               rf"$\langle{SYMB[kind]}\rangle = {stats['mean']:.3f}\pm{stats['mean_err']:.3f}$, " + chi)
        ax.text(0.03, 0.97, txt, transform=ax.transAxes, va="top", fontsize=6.5,
                bbox=dict(fc="white", ec="none", alpha=0.8, pad=1.5))
    if surf is not None:      # drawn only across the data range of each row/column (no extrapolation)
        for o in np.sort(grid[other].unique()):
            near = grid[np.abs(grid[other] - o) < 1e-9][coord]
            if len(near) < 2:
                continue
            t = np.linspace(near.min(), near.max(), 40)
            xx, yy = (t, np.full_like(t, o)) if coord == "x" else (np.full_like(t, o), t)
            val, _ = surf.exponent(xx, yy, kind)
            ax.plot(10 ** t, val, color=cmap(norm(o)), lw=0.7, ls="-", alpha=0.6, zorder=1)
    ax.set_xscale("log")
    decade_ticks(ax, "x")
    ax.set_xlabel(style.SYM[COORD[coord]])
    ax.set_ylabel(rf"${SYMB[kind]}$")
    return stats


def get_parser() -> argparse.ArgumentParser:
    """Command-line options (see the module docstring).

    Returns
    -------
    argparse.ArgumentParser
    """
    p = argparse.ArgumentParser(prog="python -m rbc.plots.exponents", description=__doc__.split("\n")[0])
    p.add_argument("--q", choices=("Nu-1", "Nu", "Re"), default="Nu-1", help="quantity Q in Q ~ Pr^alpha Ra^beta")
    p.add_argument("--method", choices=("pairs", "lines"), default="pairs",
                   help="pairs: two-point secants at midpoints; lines: one regression per row/column")
    p.add_argument("--no-cells", action="store_true", help="omit the cell-centre (diagonal) estimates")
    p.add_argument("--fit-source", choices=("edge", "line", "cell"), help="estimates used for the flatness fit")
    p.add_argument("--surface", type=int, default=0, help="overlay a polynomial surface fit of this degree")
    p.add_argument("--refs", type=float, nargs="*", default=[0.5, 1 / 3], help="reference exponents drawn as dotted lines")
    p.add_argument("--min-points", type=int, default=3, help="for --method lines")
    p.add_argument("--include-running", action="store_true", help="also use points whose cases.csv status is not 'finished'")
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

    df = rc.load_results(a.results)
    if not a.include_running:
        df = df[df["status"] == "finished"]
    d = ex.prepare(df, a.q)
    est = gather(d, a.method, not a.no_cells)
    fit_source = a.fit_source or ("edge" if a.method == "pairs" else "line")

    surf = None
    if a.surface:
        fits = ex.surfaces(d, range(1, max(a.surface, 3) + 1))
        print(f"surface fits of log10 {a.q} ({'weighted by errors' if fits and fits[0].weighted else 'errors unknown'}):")
        print(f"{'degree':>6} {'params':>6} {'rms log10':>10} {'chi2/dof':>9} {'AICc':>8}")
        for s in fits:
            chi = f"{s.chi2 / max(s.dof, 1):9.3g}" if s.weighted else f"{'-':>9}"
            print(f"{s.deg:>6} {len(s.terms):>6} {s.rms:>10.2e} {chi} {s.aic:>8.2f}")
        surf = next((s for s in fits if s.deg == a.surface), None)
    glob = ex.Surface(d, 1)
    (al, ale), (be, bee) = glob.exponent(0, 0, "alpha"), glob.exponent(0, 0, "beta")
    gof = f"chi2/dof = {glob.chi2 / max(glob.dof, 1):.3g}" if glob.weighted else f"rms log10 residual = {glob.rms:.2e}"
    print(f"single power law over all points: {a.q} ~ Pr^({al:.4f} ± {ale:.4f}) Ra^({be:.4f} ± {bee:.4f})   {gof}")

    w, _ = style.figsize(1.0)
    fig, axes = plt.subplots(2, 2, figsize=(w, 0.78 * w), sharey="row")
    cmap = colormaps["viridis"]
    norms = {"x": Normalize(d["y"].min(), d["y"].max()), "y": Normalize(d["x"].min(), d["x"].max())}
    summary = []
    for i, kind in enumerate(("alpha", "beta")):
        for j, coord in enumerate(("x", "y")):
            st = panel(axes[i, j], est, kind, coord, cmap, norms[coord], a.refs, fit_source, surf, d)
            summary.append(dict(kind=kind, vs=COORD[coord], **st))
    for j, coord in enumerate(("x", "y")):
        other = "Pr" if coord == "x" else "Ra"
        cb = fig.colorbar(ScalarMappable(norms[coord], cmap), ax=axes[:, j], location="top", shrink=0.8,
                          aspect=30, pad=0.02)
        cb.set_label(rf"$\log_{{10}}{other}$", fontsize=7)
        cb.ax.tick_params(labelsize=6)
    qlab = {"Nu-1": "Nu-1"}.get(a.q, a.q)
    fig.suptitle(rf"local exponents of ${qlab}\sim Pr^{{\alpha}}Ra^{{\beta}}$ "
                 f"({'two-point secants' if a.method == 'pairs' else 'row/column regressions'}"
                 f"{', diamonds: cell centres' if not a.no_cells else ''})", fontsize=8)
    tag = f"{a.q.replace('-', 'm')}_{a.method}" + (f"_surf{a.surface}" if a.surface else "")
    print(style.save(fig, a.out / f"exponents_{tag}", png=a.png))
    a.out.mkdir(parents=True, exist_ok=True)
    est.assign(Ra=10 ** est.x, Pr=10 ** est.y).to_csv(a.out / f"exponents_{tag}.csv", index=False, float_format="%.6g")
    cells = ex.cells(d)
    if len(cells):
        cells.to_csv(a.out / f"exponents_{tag}_cells.csv", index=False, float_format="%.6g")
    print("\nflatness fits (exponent = mean, or intercept + slope*log10 param):")
    for r in summary:
        if "mean" in r:
            print(f"  {r['kind']:>5} vs {r['vs']}: mean {r['mean']:.4f} ± {r['mean_err']:.4f}, "
                  f"slope {r['slope']:+.4f} ± {r['slope_err']:.4f} per decade, chi2_nu(flat) {r['chi2_flat']:.3g}, n={r['n']}")
    with pd.option_context("display.width", 140, "display.precision", 4):
        if len(cells):
            print("\nslope along Pr ~ Ra^-1/2 at cell centres (s = beta - alpha/2; Shishkina-Lohse 1/4):")
            print(cells[["c", "x", "y", "alpha", "beta", "s_slice", "s_slice_err"]].sort_values(["c", "x"]).to_string(index=False))


if __name__ == "__main__":
    main()
