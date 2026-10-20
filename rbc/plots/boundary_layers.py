"""Boundary layers and their resolution (figure F).

ratio : delta_u/delta_T against Pr, both by the wall-tangent definition, with the
        laminar Prandtl-Blasius asymptotes of Shishkina et al. (2010):
        delta_T/delta_u = 0.588 Pr^-1/2 (Pr << 0.046), 0.982 Pr^-1/3 (Pr >> 0.046);
        plus the number of Chebyshev nodes inside each layer against Ra.
grid  : near-wall profiles of chosen runs with every Chebyshev node drawn and the
        two layers shaded, so the reader sees how many nodes resolve each layer.

Usage::

    python -m rbc.plots.boundary_layers ratio
    python -m rbc.plots.boundary_layers grid --points "5 0" "5 -3" "8 -3"
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colormaps

from rbc import cases as rc
from rbc import style
from rbc.analysis.profiles import cheb_grid
from rbc.plots.common import decade_ticks

#: Pr where the two Prandtl–Blasius asymptotes of delta_T/delta_u cross.
PR_STAR = 0.046


def ratio_fig(df, out, png):
    """Figure ``bl_ratio``: delta_u/delta_T against Pr, and nodes per layer against Ra.

    Parameters
    ----------
    df : pandas.DataFrame
        Finished production results with boundary-layer columns.
    out : pathlib.Path
        Output folder.
    png : bool
        Also write a PNG.
    """
    fig, (a1, a2) = plt.subplots(1, 2, figsize=style.figsize(1.0, 0.42))
    cmap = colormaps["viridis"]
    lo, hi = df["log_ra"].min(), df["log_ra"].max()
    for lra, g in df.groupby("log_ra"):
        col = cmap(0.9 * (lra - lo) / max(hi - lo, 1e-9))
        g = g.sort_values("Pr")
        a1.plot(g["Pr"], g["delta_u_slope"] / g["delta_T"], "o", color=col, ms=3.5,
                label=rf"$Ra={style.pow10(lra)}$")
    pr = np.logspace(np.log10(df["Pr"].min()) - 0.3, np.log10(df["Pr"].max()) + 0.3, 200)
    lowm, highm = pr < 3 * PR_STAR, pr > PR_STAR / 3
    a1.plot(pr[lowm], pr[lowm] ** 0.5 / 0.588, "k--", lw=0.8, label=r"PB, $Pr\ll Pr^*$")
    a1.plot(pr[highm], pr[highm] ** (1 / 3) / 0.982, "k:", lw=0.8, label=r"PB, $Pr\gg Pr^*$")
    a1.set(xscale="log", yscale="log", xlabel=style.SYM["Pr"], ylabel=r"$\delta_u/\delta_T$")
    decade_ticks(a1, "x")
    a1.legend(fontsize=6.5)

    for (lpr, g), mk in zip(df.groupby("log_pr"), "os^vDP<>X"):
        g = g.sort_values("Ra")
        kw = dict(marker=mk, ms=3.5, lw=0.6)
        col = colormaps["plasma"](0.85 * (lpr - df["log_pr"].min()) / max(np.ptp(df["log_pr"]), 1e-9))
        a2.plot(g["Ra"], g["n_delta_T"], color=col, ls="-", **kw, label=rf"$Pr={style.pow10(lpr)}$")
        a2.plot(g["Ra"], g["n_delta_u_slope"], color=col, ls="--", mfc="none", **kw)
    a2.set(xscale="log", xlabel=style.SYM["Ra"], ylabel="Chebyshev nodes in layer")
    decade_ticks(a2, "x")
    a2.plot([], [], "k-", label=r"in $\delta_T$")
    a2.plot([], [], "k--", label=r"in $\delta_u$")
    a2.legend(fontsize=6, ncol=2)
    print(style.save(fig, out / "bl_ratio", png=png))


def grid_fig(df, points, out, png, span=3.0):
    """Figure ``bl_grid``: near-wall profiles with every Chebyshev node drawn.

    Parameters
    ----------
    df : pandas.DataFrame
        Finished production results.
    points : list of str
        ``"logRa logPr"`` of each panel.
    out : pathlib.Path
        Output folder.
    png : bool
        Also write a PNG.
    span : float, optional
        Height shown, in boundary-layer thicknesses.
    """
    rows = []
    for pt in points:
        lra, lpr = (float(v) for v in pt.split())
        r = df[(df["log_ra"] == lra) & (df["log_pr"] == lpr)]
        if r.empty:
            raise SystemExit(f"no production run at Ra=1e{lra}, Pr=1e{lpr}")
        rows.append(r.iloc[0])
    fig, axes = plt.subplots(1, len(rows), figsize=style.figsize(1.0, 0.5), squeeze=False)
    cT, cu = style.SLICE_STYLE[-1.0]["color"], style.SLICE_STYLE[0.0]["color"]
    for ax, r in zip(axes[0], rows):
        prof = np.load(Path(r["path"]) / "outputs" / "profiles.npz")
        z = prof["z"] + 0.5                      # distance from the bottom wall
        top = span * max(r["delta_T"], r["delta_u_slope"])
        m = z <= top * 1.05
        ax.axhspan(0, r["delta_T"], color=cT, alpha=0.15, lw=0)
        ax.axhspan(0, r["delta_u_slope"], color=cu, alpha=0.10, lw=0)
        ax.plot(2 * (1 - prof["T"][m]), z[m], color=cT, label=r"$2(1-\overline{T})$")
        ax.plot(prof["uh_rms"][m] / prof["uh_rms"].max(), z[m], color=cu, label=r"$u_{h}/u_{h,\max}$")
        nodes = cheb_grid(int(r["nz"])) + 0.5
        nodes = nodes[nodes <= top]
        for zn in nodes:
            ax.axhline(zn, color="0.6", lw=0.3, zorder=0)
        ax.plot(np.zeros_like(nodes), nodes, ls="none", marker="o", ms=1.8, color="k", clip_on=False)
        ax.text(0.03, r["delta_T"], rf"$\delta_T$: {int(r['n_delta_T'])} nodes", color=cT, fontsize=6.5,
                ha="left", va="bottom")
        ax.text(0.97, r["delta_u_slope"], rf"$\delta_u$: {int(r['n_delta_u_slope'])} nodes", color=cu,
                fontsize=6.5, ha="right", va="bottom")
        ax.text(0.97, 0.97, rf"$Ra={style.pow10(r['log_ra'])}$" "\n" rf"$Pr={style.pow10(r['log_pr'])}$"
                "\n" rf"$N_z={int(r['nz'])}$", transform=ax.transAxes, ha="right", va="top", fontsize=7,
                bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.5))
        ax.tick_params(right=False, top=False, which="both")
        ax.set(xlim=(0, 1.0), ylim=(0, top), xlabel="normalized profile")
    axes[0, 0].set_ylabel(r"distance from wall $z$")
    axes[0, 0].legend(fontsize=6.5, loc="upper left")
    print(style.save(fig, out / "bl_grid", png=png))


def get_parser() -> argparse.ArgumentParser:
    """Command-line options (see the module docstring).

    Returns
    -------
    argparse.ArgumentParser
    """
    p = argparse.ArgumentParser(prog="python -m rbc.plots.boundary_layers", description=__doc__.split("\n")[0])
    p.add_argument("figure", choices=("ratio", "grid"), help="ratio: delta_u/delta_T and nodes per layer; grid: near-wall profiles with every node")
    p.add_argument("--points", nargs="+", default=["5 0", "7 -2"], help='"logRa logPr" for grid')
    p.add_argument("--span", type=float, default=3.0, help="grid: height shown, in layer thicknesses")
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
    rc.require_columns(df, ["delta_T", "delta_u_slope", "nz"], "the boundary-layer figure")
    df = df[df["status"] == "finished"].dropna(subset=["delta_T", "delta_u_slope"])
    (ratio_fig if a.figure == "ratio" else lambda d, o, pn: grid_fig(d, a.points, o, pn, a.span))(df, a.out, a.png)


if __name__ == "__main__":
    main()