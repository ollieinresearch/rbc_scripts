"""Time-averaged horizontal-mean profiles and boundary layers of one run.

Writes ``outputs/profiles.npz`` (profiles, for figures), ``outputs/profiles.json``
(scalar diagnostics, for :mod:`rbc.collect`) and ``outputs/profiles.pdf``. Usage::

    python -m rbc.analysis.profiles RUN_DIR [--t0 T] [--t1 T] [--no-plot]

Free-fall units, as in the simulation. Boundary-layer thicknesses (mean of bottom and
top): ``delta_T`` = 1/(2 Nu) from the wall gradient; ``delta_T_rms`` and ``delta_u``,
the distance from the wall to the peak of theta_rms and of the horizontal rms velocity;
``delta_u_slope``, where the wall tangent of the rms velocity reaches its peak value.
``n_*`` is the number of Chebyshev nodes of the simulation grid inside each layer.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from numpy.polynomial import Chebyshev

from rbc import io
from rbc.cases import parse_run_path

TASKS = ["avg_T", "avg_u_sq", "avg_v_sq", "avg_w_sq"]
OPTIONAL = ["avg_T_sq", "havg_wT"]


def cheb_grid(n):
    """Dedalus ChebyshevT grid (Gauss points) on z in [-1/2, 1/2], ascending.

    Parameters
    ----------
    n : int
        Number of modes (= nodes at scale 1).

    Returns
    -------
    numpy.ndarray
        Node heights relative to mid-height.
    """
    return -0.5 * np.cos(np.pi * (np.arange(n) + 0.5) / n)


def _peak_dist(fc, n=20001):
    """Distance from each wall to the maximum of ``fc`` in that half, averaged over
    the two walls.

    Parameters
    ----------
    fc : callable
        Profile as a function of z in [-1/2, 1/2].
    n : int, optional
        Points of the search grid in each half.

    Returns
    -------
    float
    """
    zb, zt = np.linspace(-0.5, 0, n), np.linspace(0, 0.5, n)
    return 0.5 * ((zb[np.argmax(fc(zb))] + 0.5) + (0.5 - zt[np.argmax(fc(zt))]))


def compute(run: Path, t0=0.0, t1=np.inf) -> dict:
    """Time-averaged profiles and boundary-layer diagnostics.

    Parameters
    ----------
    run : pathlib.Path
        Run directory with ``horizontal_analysis/``.
    t0, t1 : float, optional
        Averaging window.

    Returns
    -------
    dict
        Arrays ``z`` (centred, in [-1/2, 1/2]), ``T``, ``w_rms``, ``uh_rms``, ``u2``,
        ``theta_rms`` and ``nu_z`` (total heat flux, flat when converged) when the
        tasks exist; and ``diag``, the scalars written to ``profiles.json``.

    Raises
    ------
    ValueError
        With fewer than three writes in the window.
    """
    hdir = run / "horizontal_analysis"
    t, d = io.load_series(hdir, TASKS + OPTIONAL, t0, t1, optional=OPTIONAL)
    if t.size < 3:
        raise ValueError(f"{run}: only {t.size} writes after t={t0}")
    z = io.task_grid(hdir, "avg_T")
    if z is None:  # older files saved z as a task
        _, zd = io.load_series(hdir, ["z_an"], t[0], t[0])
        z = zd["z_an"][0]
    order = np.argsort(z)
    z = z[order] - 0.5
    m = {k: io.time_mean(t, v)[order] for k, v in d.items()}

    T = m["avg_T"]
    out = dict(z=z, T=T, w_rms=np.sqrt(m["avg_w_sq"]),
               uh_rms=np.sqrt(m["avg_u_sq"] + m["avg_v_sq"]),
               u2=m["avg_u_sq"] + m["avg_v_sq"] + m["avg_w_sq"])
    if "avg_T_sq" in m:   # fluctuations about the time- and horizontal-mean profile
        out["theta_rms"] = np.sqrt(np.maximum(m["avg_T_sq"] - T ** 2, 0))

    p = parse_run_path(run) or {}
    Ra, Pr = 10.0 ** p.get("log_ra", np.nan), 10.0 ** p.get("log_pr", np.nan)
    nz_sim = p.get("nz") or int(round(z.size / 1.5))

    # Grid points are interior Gauss nodes: differentiate/integrate the Chebyshev
    # interpolant, which is spectrally accurate and evaluates at the walls exactly.
    cheb = lambda f: Chebyshev.fit(z, f, deg=z.size - 1, domain=[-0.5, 0.5])
    Tc = cheb(T)
    dTdz = Tc.deriv()(z)
    g_bot, g_top = np.abs(Tc.deriv()(np.array([-0.5, 0.5])))
    nu_wall = 0.5 * (g_bot + g_top)
    uc = cheb(out["uh_rms"])
    u_slope = np.abs(uc.deriv()(np.array([-0.5, 0.5]))).mean()
    bl = dict(delta_T=0.5 / nu_wall,                      # tangent at wall meets T = 1/2: 1/(2 Nu)
              delta_u=_peak_dist(uc),                      # wall to peak of horizontal rms velocity
              delta_u_slope=out["uh_rms"].max() / u_slope)  # tangent at wall meets the peak value
    if "theta_rms" in out:
        bl["delta_T_rms"] = _peak_dist(cheb(out["theta_rms"]))
    zs = cheb_grid(nz_sim)
    nodes = {f"n_{k}": int(np.sum(zs + 0.5 < v)) for k, v in bl.items()}   # nodes in the bottom layer

    diag = dict(Ra=Ra, Pr=Pr, nz=nz_sim, t0=float(t[0]), t1=float(t[-1]), n_writes=int(t.size),
                Nu_wall_bot=float(g_bot), Nu_wall_top=float(g_top),
                Re_rms=float((Ra / Pr) ** 0.5 * np.sqrt(cheb(out["u2"]).integ(lbnd=-0.5)(0.5))),
                **{k: float(v) for k, v in bl.items()}, **nodes)
    if "havg_wT" in m:
        out["nu_z"] = (Ra * Pr) ** 0.5 * m["havg_wT"] - dTdz   # total heat flux; flat when converged
        diag["Nu_flux_mean"] = float(np.mean(out["nu_z"]))
        diag["Nu_flux_spread_pct"] = float(100 * np.ptp(out["nu_z"]) / np.mean(out["nu_z"]))
    out["diag"] = diag
    return out


def plot(prof: dict, path: Path) -> None:
    """Diagnostic figure ``outputs/profiles.pdf``: one panel per profile.

    Parameters
    ----------
    prof : dict
        Output of :func:`compute`.
    path : pathlib.Path
        File to write.
    """
    import matplotlib.pyplot as plt
    from rbc import style

    style.use("diag")
    S = style.SYM
    panels = [
        ("T", r"$\overline{T}$"),
        ("w_rms", r"$\sqrt{\overline{w^2}}$"),
        ("uh_rms", r"$\sqrt{\overline{u^2+v^2}}$"),
        ("u2", r"$\overline{u^2+v^2+w^2}$"),
    ]
    if "theta_rms" in prof:
        panels.insert(1, ("theta_rms", r"$\theta_{\rm rms}$"))
    if "nu_z" in prof:
        prof = dict(prof, nu_dev=100 * (prof["nu_z"] / prof["diag"]["Nu_flux_mean"] - 1))
        panels.append(("nu_dev", r"$Nu(z)/\langle Nu\rangle_z - 1$ (\%)"))
    fig, axes = plt.subplots(1, len(panels), figsize=(3.2 * len(panels), 4.2), sharey=True)
    for ax, (k, lab) in zip(axes, panels):
        ax.plot(prof[k], prof["z"])
        ax.set_xlabel(lab)
    axes[0].set_ylabel(S["z"])
    axes[0].set_ylim(-0.5, 0.5)
    axes[0].set_xlim(0, 1)
    dg = prof["diag"]
    if "nu_z" in prof:
        axes[-1].axvline(0, color="0.5", lw=0.8, ls="--")
        axes[-1].set_title(rf"$\langle Nu\rangle_z={dg['Nu_flux_mean']:.5g}$")
    fig.suptitle(f"$t\\in[{dg['t0']:.1f},{dg['t1']:.1f}]$, "
                 f"$\\delta_T={dg['delta_T']:.4f}$ ({dg['n_delta_T']} nodes), "
                 f"$\\delta_u={dg['delta_u']:.4f}$ ({dg['n_delta_u']} nodes), $Re={dg['Re_rms']:.4g}$")
    fig.savefig(path)
    plt.close(fig)


def get_parser() -> argparse.ArgumentParser:
    """Command-line options.

    Returns
    -------
    argparse.ArgumentParser
    """
    p = argparse.ArgumentParser(prog="python -m rbc.analysis.profiles",
                                description="Time-averaged horizontal profiles and boundary layers.")
    p.add_argument("run", type=Path, help="run directory")
    p.add_argument("--t0", "--time", type=float, default=0.0, dest="t0", help="start of the averaging window")
    p.add_argument("--t1", type=float, default=np.inf, help="end of the averaging window")
    p.add_argument("--no-plot", action="store_true", help="skip outputs/profiles.pdf")
    return p


def main(argv=None) -> None:
    """Compute and write one run's profiles.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments (default: ``sys.argv[1:]``).
    """
    a = get_parser().parse_args(argv)
    run = a.run.resolve()
    prof = compute(run, a.t0, a.t1)
    out = run / "outputs"
    out.mkdir(exist_ok=True)
    np.savez(out / "profiles.npz", **{k: v for k, v in prof.items() if k != "diag"})
    (out / "profiles.json").write_text(json.dumps(prof["diag"], indent=1))
    if not a.no_plot:
        plot(prof, out / "profiles.pdf")
    print(json.dumps(prof["diag"], indent=1))


if __name__ == "__main__":
    main()
