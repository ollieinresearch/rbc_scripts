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
from matplotlib.cm import ScalarMappable

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
        Columns ``kind``, ``x``, ``y``, ``value``, ``err``, ``Lmin``, ``source`` (and
        ``lo``, ``hi`` for lines).
    """
    parts = [ex.edges(d) if method == "pairs" else ex.lines(d)]
    if use_cells:
        c = ex.cells(d)
        if len(c):
            for k in ("alpha", "beta"):
                parts.append(pd.DataFrame(dict(kind=k, x=c.x, y=c.y, value=c[k], err=c[f"{k}_err"],
                                               Lmin=c.Lmin, source="cell")))
    return pd.concat(parts, ignore_index=True)


def panel(ax, est, kind, coord, cmap, norm, refs, fit_source, surf=None, grid=None, lmin=-np.inf,
          ref_labels=True):
    """One exponent against one parameter, coloured by the other.

    Two-point estimates along the same row (same value of the other parameter) are
    joined by a thin line in that row's colour, so each row's trend can be followed;
    cell-centre estimates are open diamonds. Estimates involving a run with
    log10 Q < ``lmin`` are drawn as grey crosses and left out of the trend fit and of
    the y range.

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
        Estimates used for the trend fit.
    surf : rbc.analysis.exponents.Surface, optional
        Surface fit whose exponents are drawn along each row or column (within the data).
    grid : pandas.DataFrame, optional
        The data points (needed with ``surf``).
    lmin : float, optional
        Smallest log10 Q of the runs an estimate may use.
    ref_labels : bool, optional
        Label the reference lines at the right edge.

    Returns
    -------
    dict
        Output of :func:`rbc.analysis.exponents.flatness` for the kept estimates (empty
        if too few).
    """
    other = "y" if coord == "x" else "x"
    e = est[est.kind == kind]
    ok = e["Lmin"] >= lmin if "Lmin" in e else pd.Series(True, index=e.index)
    for v in refs:
        ax.axhline(v, color="0.55", lw=0.6, ls=":", zorder=1)
        if ref_labels:
            ax.annotate(frac(v), (1.0, v), xycoords=("axes fraction", "data"), xytext=(2, 0),
                        textcoords="offset points", va="center", fontsize=7, color="0.35",
                        annotation_clip=False)
    if surf is not None:      # drawn only across the data range of each row/column (no extrapolation)
        for o in np.sort(grid[other].unique()):
            near = grid[np.abs(grid[other] - o) < 1e-9][coord]
            if len(near) < 2:
                continue
            t = np.linspace(near.min(), near.max(), 40)
            xx, yy = (t, np.full_like(t, o)) if coord == "x" else (np.full_like(t, o), t)
            val, _ = surf.exponent(xx, yy, kind)
            ax.plot(10 ** t, val, color=cmap(norm(o)), lw=0.6, ls="--", alpha=0.7, zorder=1)
    for src in ("edge", "line"):     # one thin line per row through its kept estimates
        g = e[(e.source == src) & ok]
        for o, r in g.groupby(other):
            r = r.sort_values(coord)
            if len(r) > 1:
                ax.plot(10 ** r[coord], r["value"], color=cmap(norm(o)), lw=0.7, alpha=0.8, zorder=2)
    for src, mk, ms, fill in (("edge", "o", 3.2, True), ("line", "o", 3.6, True), ("cell", "D", 2.4, False)):
        g = e[(e.source == src) & ok]
        for _, r in g.iterrows():
            col = cmap(norm(r[other]))
            ax.errorbar(10 ** r[coord], r["value"], yerr=r["err"] if np.isfinite(r["err"]) else None,
                        fmt=mk, ms=ms, color=col, mfc=col if fill else "white", mew=0.7, elinewidth=0.6,
                        capsize=1.2, zorder=3)
    bad = e[~ok]
    if len(bad):
        ax.plot(10 ** bad[coord], bad["value"], ls="none", marker="x", ms=3, mew=0.6, color="0.6", zorder=2)
        kept = e[ok]["value"]
        if len(kept):
            lo, hi = kept.min(), kept.max()
            pad = 0.12 * max(hi - lo, 0.05)
            lo, hi = min(lo - pad, min(refs, default=lo)), max(hi + pad, max(refs, default=hi))
            ax.set_ylim(lo - 0.02, hi + 0.02)
    fs = e[(e.source == fit_source) & ok]
    stats = ex.flatness(fs[coord], fs["value"], fs["err"]) if len(fs) > 2 else {}
    if stats:
        t = np.linspace(fs[coord].min(), fs[coord].max(), 50)
        ax.plot(10 ** t, stats["slope"] * t + stats["intercept"], color="k", lw=0.9, zorder=4)
        ax.set_title(rf"trend ${stats['slope']:+.3f}\pm{stats['slope_err']:.3f}$ per decade, "
                     rf"mean ${stats['mean']:.3f}$", fontsize=7, pad=3)
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
    p.add_argument("--draw-surface", action="store_true",
                   help="with --surface: also draw the surface's exponents along each row (dashed)")
    p.add_argument("--min-q", type=float, default=0.0,
                   help="leave out (grey crosses, not fitted) estimates using a run with Q below this, "
                        "e.g. 1 to drop the weakly convective corner where Nu-1 < 1")
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
    rc.require_columns(df, ["Nu" if a.q == "Nu-1" else a.q], "the exponent figure")
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
    fig, axes = plt.subplots(2, 2, figsize=(w, 0.88 * w), sharey="row", sharex="col")
    # points against Ra are coloured by their Pr, points against Pr by their Ra, both from red
    # (cheap) to violet (expensive), with the same colour for a value as in every figure
    scales = {"x": style.expense_norm("pr"), "y": style.expense_norm("ra")}
    lmin = np.log10(a.min_q) if a.min_q > 0 else -np.inf
    summary = []
    for i, kind in enumerate(("alpha", "beta")):
        for j, coord in enumerate(("x", "y")):
            st = panel(axes[i, j], est, kind, coord, *scales[coord], a.refs, fit_source,
                       surf if a.draw_surface else None, d, lmin, ref_labels=(j == 1))
            summary.append(dict(kind=kind, vs=COORD[coord], **st))
        lims = [ax.get_ylim() for ax in axes[i]]
        axes[i, 0].set_ylim(min(l[0] for l in lims), max(l[1] for l in lims))
    for ax in axes[0]:
        ax.set_xlabel("")
    for j, coord in enumerate(("x", "y")):
        other = "Pr" if coord == "x" else "Ra"
        vals = np.sort(d["y" if coord == "x" else "x"].unique())
        cmap, norm = scales[coord]
        cb = fig.colorbar(ScalarMappable(norm, cmap), ax=axes[:, j], location="top", shrink=0.85,
                          aspect=32, pad=0.02, ticks=vals[::2] if len(vals) > 7 else vals)
        cb.set_label(rf"colour: $\log_{{10}}{other}$", fontsize=7)
        cb.ax.tick_params(labelsize=6)
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], marker="o", ls="-", lw=0.7, ms=3, color="0.3", label="two-point, joined by row"
                      if a.method == "pairs" else "row/column fit"),
               Line2D([], [], marker="D", ls="none", ms=2.6, mfc="white", color="0.3", label="cell centre"),
               Line2D([], [], color="k", lw=0.9, label="weighted trend" if d["sL"].notna().all() else "trend (unweighted)")]
    if np.isfinite(lmin):
        handles.append(Line2D([], [], marker="x", ls="none", ms=3, color="0.6",
                              label=rf"left out: ${'Nu-1' if a.q == 'Nu-1' else a.q} < {a.min_q:g}$"))
    if a.draw_surface and surf is not None:
        handles.append(Line2D([], [], ls="--", lw=0.6, color="0.3", label=f"degree-{a.surface} surface"))
    fig.legend(handles=handles, loc="outside lower center", ncol=len(handles), fontsize=6.5,
               handlelength=1.6, columnspacing=1.2)
    qlab = {"Nu-1": "Nu-1"}.get(a.q, a.q)
    fig.suptitle(rf"Local exponents of ${qlab}\sim Pr^{{\alpha}}Ra^{{\beta}}$", fontsize=8)
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