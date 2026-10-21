"""Nusselt and Reynolds numbers of one run, with statistical-convergence diagnostics.

Writes ``outputs/nusselt.json``, ``outputs/info.txt`` and ``outputs/info.pdf``. Usage::

    python -m rbc.analysis.nusselt RUN_DIR [--t0 T | --t0 auto] [--t1 T] [--no-plot]

In free-fall units, with ⟨·⟩ the volume average and an overline the time average over
the window [t0, t1], three exact identities give Nu in a statistically steady state::

    Nu = 1 + sqrt(Ra Pr) ⟨wT⟩ = ⟨|∇T|²⟩ = 1 + Pr ⟨|ω|²⟩,     Re = sqrt(Ra/Pr ⟨|u|²⟩).

The plate fluxes Nu_bot = -∂T/∂z(0) and Nu_top = -∂T/∂z(1), from the horizontally
averaged temperature profile, are a fourth estimate and a stationarity check.

Convergence is judged on Nu - 1 from ⟨wT⟩ and from ``⟨|∇T|²⟩`` only (``CONVERGENCE_QOI``).
The vorticity identity is computed and reported but not checked: ``⟨|ω|²⟩`` is dominated by
the smallest scales, which can be under-resolved while the velocity and temperature
spectra, and so the heat transport, are not. Its difference from ⟨wT⟩ is reported as
``vort_diff_pct``, a resolution diagnostic rather than a convergence one.

Convergence diagnostics (see the "Statistical convergence" page of the documentation):

``spread_pct``
    largest difference of Nu - 1 between the three equal sections of the window, in
    percent of Nu - 1 (the group's original rule: below 1% for ⟨wT⟩ and ``⟨|∇T|²⟩``).
``se``, ``se_tau``
    standard error of the time mean: batch means over 10 blocks, and from the
    integrated autocorrelation time ``tau_int`` (Sokal's automatic window). ``n_eff``
    = T / (2 tau_int) is the number of effectively independent samples.
``p_sections``
    chi-square test that the three section means differ only by sampling noise.
``drift_z``
    fitted linear trend over the window divided by its standard error.
``t0_mser``
    MSER-5 estimate of where the initial transient ends (White 1997), the later of the
    two checked identities.
``plate_imbalance_z``
    mean of Nu_bot - Nu_top over the window divided by its standard error. In a
    steady state it vanishes, since d⟨T⟩/dt = kappa (Nu_bot - Nu_top).
``converged``
    the original thirds-within-1% rule; ``converged_stat`` the proposed statistical
    rule (thresholds in ``CRITERIA``).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid, trapezoid

from rbc import io
from rbc.cases import parse_run_path

N_SECTIONS, N_BLOCKS = 3, 10
QOI = {  # task -> (Nu from time-mean of task, label)
    "avg_wT": (lambda m, Ra, Pr: 1 + np.sqrt(Ra * Pr) * m, r"$\langle wT\rangle$"),
    "avg_grad_T_sq": (lambda m, Ra, Pr: m, r"$\langle|\nabla T|^2\rangle$"),
    "avg_vorticity_sq": (lambda m, Ra, Pr: 1 + Pr * m, r"$\langle|\omega|^2\rangle$"),
}
#: Identities whose convergence is checked (``converged``, ``converged_stat``, the MSER
#: window start and ``qoi_spread_pct``). The vorticity identity is reported only.
CONVERGENCE_QOI = ("avg_wT", "avg_grad_T_sq")
#: Thresholds of the proposed statistical convergence rule (``converged_stat``).
CRITERIA = dict(rel_se=0.005, n_eff=50, p_sections=0.01, drift_z=3.0, plate_z=3.0)


def section_bounds(t: np.ndarray, n: int) -> list[tuple[int, int]]:
    """Split a time series into equal-duration sections.

    Parameters
    ----------
    t : numpy.ndarray
        Increasing sample times.
    n : int
        Number of sections.

    Returns
    -------
    list of tuple of int
        Inclusive ``(i0, i1)`` index pairs; neighbouring sections share an endpoint.
    """
    edges = t[0] + (t[-1] - t[0]) * np.arange(1, n) / n
    idx = np.searchsorted(t, edges)
    return list(zip(np.r_[0, idx], np.r_[idx, len(t) - 1]))


def tmean(t: np.ndarray, x: np.ndarray, i0: int = 0, i1: int | None = None) -> float:
    """Trapezoidal time average of ``x`` between samples ``i0`` and ``i1``.

    Parameters
    ----------
    t, x : numpy.ndarray
        Sample times and values.
    i0, i1 : int, optional
        Inclusive index range (default: the whole series).

    Returns
    -------
    float
        The time average.
    """
    i1 = len(t) - 1 if i1 is None else i1
    s = slice(i0, i1 + 1)
    return trapezoid(x[s], t[s]) / (t[i1] - t[i0])


def run_params(run: Path, f_tasks: dict) -> tuple[float, float]:
    """Ra and Pr of a run, from its directory name or else from the analysis tasks.

    Parameters
    ----------
    run : pathlib.Path
        Run directory (``.../ra1e8/pr1e0/108_Gam2_2``).
    f_tasks : dict
        Loaded analysis tasks; ``"Ra"`` and ``"Pr"`` are used if the path does not parse.

    Returns
    -------
    tuple of float
        ``(Ra, Pr)``.
    """
    p = parse_run_path(run)
    if p:
        return 10.0 ** p["log_ra"], 10.0 ** p["log_pr"]
    return float(np.ravel(f_tasks["Ra"])[-1]), float(np.ravel(f_tasks["Pr"])[-1])


def _uniform(t: np.ndarray, x: np.ndarray) -> tuple[np.ndarray, float]:
    """Resample onto a uniform grid at the median sample spacing (writes are
    nearly uniform: they land on the first step after each analysis interval)."""
    dt = float(np.median(np.diff(t)))
    tu = np.arange(t[0], t[-1] + 0.5 * dt, dt)
    return np.interp(tu, t, x), dt


def autocorr_time(t: np.ndarray, x: np.ndarray, c: float = 5.0) -> tuple[float, float]:
    """Integrated autocorrelation time and the standard error of the time mean.

    With rho(s) the autocorrelation function, tau_int = ∫₀^∞ rho(s) ds, and for a
    record of length T, Var(mean) ≈ 2 tau_int Var(x) / T. The sum is truncated at the
    smallest window W with W ≥ c tau_int(W) (Sokal 1997), which balances the bias of
    stopping early against the noise of summing the tail of rho.

    Parameters
    ----------
    t, x : numpy.ndarray
        Sample times and values (near-uniform spacing).
    c : float, optional
        Window constant (5 is the usual choice).

    Returns
    -------
    tau_int : float
        Integrated autocorrelation time, in simulation time units.
    se : float
        Standard error of the time mean of ``x``.
    """
    xu, dt = _uniform(t, x)
    n = xu.size
    y = xu - xu.mean()
    var = float(y @ y / n)
    if n < 4 or var == 0:
        return float("nan"), float("nan")
    f = np.fft.rfft(y, 2 * n)
    acf = np.fft.irfft(f * np.conj(f))[:n]
    acf /= acf[0]
    taus = 0.5 + np.cumsum(acf[1:])                    # tau(W) in samples, W = 1, 2, ...
    ok = np.arange(1, n) >= c * taus
    tau = float(taus[np.argmax(ok)] if ok.any() else taus[-1])
    tau = max(tau, 0.5)
    return tau * dt, float(np.sqrt(2 * tau * var / n))


def mser(t: np.ndarray, x: np.ndarray, batch: int = 5, max_frac: float = 0.5) -> float:
    """Truncation point for the initial transient by MSER-5 (White 1997).

    The series is averaged in batches of ``batch`` samples (after resampling to at
    most ~5000 batches); for each candidate start d the statistic
    Σ_{i≥d} (x_i - x̄_d)² / (n - d)² is computed, and its minimizer over the first
    ``max_frac`` of the record is returned. It balances the bias of keeping transient
    data against the variance of discarding steady data.

    Parameters
    ----------
    t, x : numpy.ndarray
        Sample times and values over the whole run.
    batch : int, optional
        Samples per batch.
    max_frac : float, optional
        Largest fraction of the record that may be discarded.

    Returns
    -------
    float
        Suggested start time of the averaging window.
    """
    xu, dt = _uniform(t, x)
    stride = max(1, xu.size // (5000 * batch))
    xu, dt = xu[::stride], dt * stride
    nb = xu.size // batch
    if nb < 10:
        return float(t[0])
    xb = xu[: nb * batch].reshape(nb, batch).mean(axis=1)
    rev = xb[::-1]
    k = np.arange(1, nb + 1)
    s1 = np.cumsum(rev)[::-1]
    s2 = np.cumsum(rev ** 2)[::-1]
    n_rem = k[::-1]
    stat = (s2 - s1 ** 2 / n_rem) / n_rem ** 2
    dmax = max(1, int(max_frac * nb))
    d = int(np.argmin(stat[:dmax]))
    return float(t[0] + d * batch * dt)


def wall_gradient_weights(z: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Linear functionals giving dT/dz at z = 0 and z = 1 from values at the nodes.

    The profile is represented by its Chebyshev interpolant through the nodes,
    T(ζ) = Σ c_k T_k(ζ) with ζ = 2z - 1, so c = V⁻¹ T for the Vandermonde matrix V.
    Since T_k'(±1) = (±1)^(k+1) k², dT/dz(wall) = 2 d·c = (2 V⁻ᵀ d)·T.

    Parameters
    ----------
    z : numpy.ndarray
        Node heights in [0, 1] (any order).

    Returns
    -------
    tuple of numpy.ndarray
        Weights ``(w0, w1)`` with dT/dz(0) = w0 @ T and dT/dz(1) = w1 @ T.
    """
    zeta = 2 * np.asarray(z, float) - 1
    n = zeta.size
    V = np.polynomial.chebyshev.chebvander(zeta, n - 1)
    k = np.arange(n, dtype=float)
    d0 = (-1.0) ** (k + 1) * k ** 2
    d1 = k ** 2
    return 2 * np.linalg.solve(V.T, d0), 2 * np.linalg.solve(V.T, d1)


def plate_fluxes(run: Path, t0: float, t1: float) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
    """Time series of the plate Nusselt numbers from ``horizontal_analysis``.

    Parameters
    ----------
    run : pathlib.Path
        Run directory.
    t0, t1 : float
        Time window.

    Returns
    -------
    tuple of numpy.ndarray or None
        ``(t, Nu_bot, Nu_top)``, or None if the profiles are not available.
    """
    hdir = run / "horizontal_analysis"
    if not hdir.exists():
        return None
    try:
        t, d = io.load_series(hdir, ["avg_T"], t0, t1)
        z = io.task_grid(hdir, "avg_T")
    except (FileNotFoundError, KeyError, IndexError):
        return None
    if z is None or t.size < 3:
        return None
    w0, w1 = wall_gradient_weights(z)
    T = d["avg_T"].reshape(t.size, -1)
    return t, -(T @ w0), -(T @ w1)


def _stats(t: np.ndarray, y: np.ndarray, value: float, scale: float) -> dict:
    """Convergence statistics of one Nu (or Re) series ``y`` with mean ``value``;
    relative quantities are relative to ``scale`` (Nu - 1 for Nu)."""
    secs, blocks = section_bounds(t, N_SECTIONS), section_bounds(t, N_BLOCKS)
    sec = np.array([tmean(t, y, a, b) for a, b in secs])
    blk = np.array([tmean(t, y, a, b) for a, b in blocks])
    se_b = float(blk.std(ddof=1) / np.sqrt(blk.size))
    tau, se_t = autocorr_time(t, y)
    T = t[-1] - t[0]
    from scipy.stats import chi2
    se_ref = se_t if np.isfinite(se_t) and se_t > 0 else se_b
    chi = float(np.sum((sec - sec.mean()) ** 2) / (N_SECTIONS * se_ref ** 2)) if se_ref > 0 else np.inf
    tc = np.array([0.5 * (t[a] + t[b]) for a, b in blocks])
    sb = se_ref * np.sqrt(N_BLOCKS)
    A = np.vstack([tc - tc.mean(), np.ones_like(tc)]).T / sb
    coef, *_ = np.linalg.lstsq(A, blk / sb, rcond=None)
    slope_se = float(np.sqrt(np.linalg.inv(A.T @ A)[0, 0]))
    return dict(value=value, sections=sec.tolist(),
                spread_pct=float(100 * np.ptp(sec) / scale),
                se=se_b, se_tau=float(se_t), tau_int=float(tau),
                n_eff=float(T / (2 * tau)) if tau > 0 else float("nan"),
                rel_se=float(se_ref / scale),
                p_sections=float(chi2.sf(chi, N_SECTIONS - 1)),
                drift_pct=float(100 * coef[0] * T / scale),
                drift_z=float(coef[0] / slope_se) if slope_se > 0 else float("nan"))


def _passes(s: dict) -> bool:
    """Whether one series' statistics meet :data:`CRITERIA` (plates checked separately)."""
    c = CRITERIA
    return bool(s["rel_se"] <= c["rel_se"] and s["n_eff"] >= c["n_eff"]
                and s["p_sections"] >= c["p_sections"] and abs(s["drift_z"]) < c["drift_z"])


def compute(run: Path, t0: float | str = 0.0, t1: float = np.inf) -> tuple[dict, dict]:
    """Nu (three identities and the plates), Re and convergence diagnostics.

    Parameters
    ----------
    run : pathlib.Path
        Run directory containing ``analysis/`` (and optionally ``horizontal_analysis/``).
    t0 : float or "auto", optional
        Start of the averaging window; ``"auto"`` uses the MSER-5 estimate (the later of
        the two checked Nu series).
    t1 : float, optional
        End of the averaging window.

    Returns
    -------
    res : dict
        Results written to ``outputs/nusselt.json``. Per identity (keys ``avg_wT``,
        ``avg_grad_T_sq`` and, if written, ``avg_vorticity_sq``): ``Nu``, ``sections``,
        ``spread_pct``, ``se``, ``se_tau``, ``tau_int``, ``n_eff``, ``rel_se``,
        ``p_sections``, ``drift_pct``, ``drift_z``. Top level: ``Nu`` and ``Nu_se`` (the
        ⟨wT⟩ values), ``Re``, ``plates``, ``qoi_spread_pct`` (⟨wT⟩ against ``⟨|∇T|²⟩``),
        ``vort_diff_pct``, ``t0_mser``, ``converged``, ``converged_stat``,
        ``convergence_qoi``.
    series : dict
        Time series for :func:`plot`.
    """
    tasks = list(QOI) + ["avg_K", "Ra", "Pr"]
    t_all, d_all = io.load_series(run / "analysis", tasks,
                                  optional=("avg_K", "Ra", "Pr", "avg_vorticity_sq"))
    Ra, Pr = run_params(run, d_all)
    qois = [k for k in QOI if k in d_all]
    nu_all = {k: QOI[k][0](d_all[k].ravel(), Ra, Pr) for k in qois}
    t0_mser = max(mser(t_all, nu_all[k]) for k in CONVERGENCE_QOI)
    if isinstance(t0, str):
        if t0 != "auto":
            raise ValueError(f"t0 must be a number or 'auto', not {t0!r}")
        t0 = t0_mser
    w = (t_all >= t0) & (t_all <= t1)
    t = t_all[w]
    if t.size < 2 * N_BLOCKS:
        raise ValueError(f"{run}: only {t.size} samples in [{t0}, {t1}]")

    res = dict(Ra=Ra, Pr=Pr, t_start=float(t_all[0]), t_end=float(t_all[-1]), t0=float(t[0]), t1=float(t[-1]),
               T_avg=float(t[-1] - t[0]), n_samples=int(t.size), t0_mser=float(t0_mser), criteria=CRITERIA,
               convergence_qoi=list(CONVERGENCE_QOI))
    series = dict(t=t, t_all=t_all, nu_all=nu_all)
    for task in qois:
        f = QOI[task][0]
        x = d_all[task][w].ravel()
        y = f(x, Ra, Pr)
        nu = float(f(tmean(t, x), Ra, Pr))
        s = _stats(t, y, nu, nu - 1)
        s["Nu"] = s.pop("value")
        res[task] = s
        series[task] = (y, f(cumulative_trapezoid(x, t) / (t[1:] - t[0]), Ra, Pr))
    if "avg_K" in d_all:
        K = d_all["avg_K"][w].ravel()
        re_inst = np.sqrt(Ra / Pr * K)
        re = float(np.sqrt(Ra / Pr * tmean(t, K)))
        s = _stats(t, re_inst, re, re)
        s["Re"] = s.pop("value")
        res["Re"] = s
        series["avg_K"] = (K, np.sqrt(Ra / Pr * cumulative_trapezoid(K, t) / (t[1:] - t[0])))
        series["K_all"] = d_all["avg_K"].ravel()

    nus = np.array([res[k]["Nu"] for k in CONVERGENCE_QOI]) - 1
    res["Nu"] = res["avg_wT"]["Nu"]
    res["Nu_se"] = res["avg_wT"]["se"]
    res["qoi_spread_pct"] = float(100 * (nus.max() - nus.min()) / nus.mean())
    if "avg_vorticity_sq" in res:
        res["vort_diff_pct"] = float(100 * (res["avg_vorticity_sq"]["Nu"] - res["Nu"]) / (res["Nu"] - 1))
    res["converged"] = bool(all(res[k]["spread_pct"] < 1.0 for k in CONVERGENCE_QOI))

    ok_plates = True
    pf = plate_fluxes(run, t[0], t[-1])
    if pf is not None:
        tp, nb, nt = pf
        diff = nb - nt
        _, se_d = autocorr_time(tp, diff)
        imb = tmean(tp, diff)
        res["plates"] = dict(Nu_bot=float(tmean(tp, nb)), Nu_top=float(tmean(tp, nt)),
                             imbalance_pct=float(100 * imb / (res["Nu"] - 1)),
                             imbalance_z=float(imb / se_d) if se_d > 0 else float("nan"))
        ok_plates = abs(res["plates"]["imbalance_z"]) < CRITERIA["plate_z"]
        series["plates"] = (tp, nb, nt)
    res["converged_stat"] = bool(all(_passes(res[k]) for k in CONVERGENCE_QOI) and ok_plates)
    series["t0_scan"] = t0_scan(t_all, nu_all["avg_wT"], t1)
    return res, series


def t0_scan(t: np.ndarray, y: np.ndarray, t1: float = np.inf, n: int = 40) -> tuple[np.ndarray, ...]:
    """Mean of ``y`` over [t0, t1] as a function of the window start t0.

    If the transient is over by some t0, the curve is flat from there on (within the
    standard error, which grows as the window shrinks).

    Parameters
    ----------
    t, y : numpy.ndarray
        Times and values over the whole run.
    t1 : float, optional
        End of every window.
    n : int, optional
        Number of start times, spread over the first half of the record.

    Returns
    -------
    tuple of numpy.ndarray
        ``(t0, mean, se)``.
    """
    m = t <= t1
    t, y = t[m], y[m]
    starts = np.linspace(t[0], t[0] + 0.5 * (t[-1] - t[0]), n)
    _, se_full = autocorr_time(t, y)
    T = t[-1] - t[0]
    means, ses = [], []
    for s in starts:
        k = t >= s
        means.append(tmean(t[k], y[k]))
        ses.append(se_full * np.sqrt(T / (t[-1] - s)))
    return starts, np.array(means), np.array(ses)


def write_info(res: dict, path: Path) -> None:
    """Write the human-readable summary ``outputs/info.txt``.

    Parameters
    ----------
    res : dict
        Output of :func:`compute`.
    path : pathlib.Path
        File to write.
    """
    c = res.get("criteria", CRITERIA)
    L = [f"Ra = {res['Ra']:.4e}, Pr = {res['Pr']:.4e}",
         f"Simulation time: {res['t_start']:.4f} to {res['t_end']:.4f}",
         f"Averaging window: {res['t0']:.4f} to {res['t1']:.4f} ({res['T_avg']:.2f} time units, "
         f"{res['n_samples']} samples)",
         f"MSER-5 estimate of the end of the transient: t = {res['t0_mser']:.2f}",
         "-" * 80]
    for k in QOI:
        if k not in res:
            continue
        r = res[k]
        L += [f"QoI {k}:" + ("" if k in CONVERGENCE_QOI else " (reported only, not a convergence check)")]
        L += [f"  section {i + 1}/{N_SECTIONS}: Nu = {v:.4f}" for i, v in enumerate(r["sections"])]
        L += [f"  max % difference of Nu-1 over sections: {r['spread_pct']:.4f}%",
              f"  Nu (full window): {r['Nu']:.6f}",
              f"  standard error: {r['se']:.2g} (batch means), {r['se_tau']:.2g} (autocorrelation); "
              f"{100 * r['rel_se']:.3f}% of Nu-1",
              f"  tau_int = {r['tau_int']:.3g}, effective samples = {r['n_eff']:.0f}",
              f"  sections consistent: p = {r['p_sections']:.3f}; "
              f"trend over window {r['drift_pct']:+.3f}% of Nu-1 ({r['drift_z']:+.1f} sigma)", "-" * 80]
    if "plates" in res:
        pl = res["plates"]
        L += [f"Plate fluxes: Nu_bot = {pl['Nu_bot']:.6f}, Nu_top = {pl['Nu_top']:.6f}, "
              f"imbalance {pl['imbalance_pct']:+.3f}% of Nu-1 ({pl['imbalance_z']:+.1f} sigma)"]
    L += [f"Percent difference in Nusselt-1, <wT> vs <|grad T|^2>: {res['qoi_spread_pct']:.4f}%"]
    if "vort_diff_pct" in res:
        L += [f"Vorticity identity vs <wT> (resolution diagnostic): {res['vort_diff_pct']:+.4f}% of Nu-1"]
    L += [f"Final Nusselt number: {res['Nu']:.6f}"]
    if "Re" in res:
        L += [f"Final Reynolds number: {res['Re']['Re']:.6f}"]
    L += ["-" * 80,
          "Convergence checks use <wT> and <|grad T|^2> only.",
          f"Converged (sections within 1%): {res['converged']}",
          f"Converged (statistical rule: SE <= {100 * c['rel_se']:g}% of Nu-1, n_eff >= {c['n_eff']}, "
          f"p_sections >= {c['p_sections']}, |drift| < {c['drift_z']} sigma, "
          f"|plate imbalance| < {c['plate_z']} sigma): {res['converged_stat']}"]
    path.write_text("\n".join(L) + "\n")


def plot(res: dict, s: dict, path: Path, max_pts: int = 20000) -> None:
    """Diagnostic figure ``outputs/info.pdf``.

    One row per identity (instantaneous Nu with section means | running mean), then
    the kinetic energy and Re, the plate fluxes, and the window-start scan.

    Parameters
    ----------
    res, s : dict
        Outputs of :func:`compute`.
    path : pathlib.Path
        File to write.
    max_pts : int, optional
        Series are thinned to about this many points.
    """
    import matplotlib.pyplot as plt
    from rbc import style

    style.use("diag")
    rows = [k for k in list(QOI) + ["avg_K"] if k in s] + ["scan"] + (["plates"] if "plates" in s else [])
    fig, axes = plt.subplots(len(rows), 2, figsize=(14, 2.6 * len(rows)))
    t = s["t"]
    k = max(1, t.size // max_pts)
    for ax, key in zip(axes, rows):
        if key == "scan":
            t0s, m, se = s["t0_scan"]
            ax[0].plot(s["t_all"][::k], s["nu_all"]["avg_wT"][::k], lw=0.5)
            ax[0].axvline(res["t0"], color="0.4", ls="--", lw=0.8, label="window start")
            ax[0].axvline(res["t0_mser"], color="C3", ls=":", lw=1, label="MSER-5")
            ax[0].set_ylabel(r"$Nu$ via $\langle wT\rangle$, whole run")
            ax[0].legend(fontsize=8)
            ax[1].fill_between(t0s, m - 2 * se, m + 2 * se, alpha=0.25, lw=0)
            ax[1].plot(t0s, m, ".-")
            ax[1].axvline(res["t0_mser"], color="C3", ls=":", lw=1)
            ax[1].set_ylabel(r"mean over $[t_0, t_1]$")
            ax[1].set_xlabel(r"window start $t_0$")
            continue
        if key == "plates":
            tp, nb, nt = s["plates"]
            kk = max(1, tp.size // max_pts)
            ax[0].plot(tp[::kk], nb[::kk], lw=0.5, label="bottom")
            ax[0].plot(tp[::kk], nt[::kk], lw=0.5, label="top")
            ax[0].set_ylabel("plate $Nu$")
            ax[0].legend(fontsize=8)
            d = nb - nt
            ax[1].plot(tp[1:][::kk], cumulative_trapezoid(d, tp)[::kk] / (tp[1:][::kk] - tp[0]))
            ax[1].axhline(0, color="0.5", lw=0.8)
            pl = res["plates"]
            ax[1].set_ylabel(r"running mean of $Nu_{bot}-Nu_{top}$")
            ax[1].set_title(f"imbalance {pl['imbalance_pct']:+.3f}% ({pl['imbalance_z']:+.1f}$\\sigma$)")
            continue
        inst, cumu = s[key]
        if key == "avg_K":
            ax[0].plot(s["t_all"][::k], s["K_all"][::k], lw=0.6)
            ax[0].axvline(t[0], color="0.5", ls="--", lw=0.8)
            ax[0].set_ylabel(r"$\langle|u|^2\rangle$")
            ax[1].plot(t[1::k], cumu[::k])
            ax[1].axhline(res["Re"]["Re"], color="C1", lw=0.8)
            ax[1].set_ylabel(r"running mean $Re$")
            continue
        r = res[key]
        ax[0].plot(t[::k], inst[::k], lw=0.6)
        for (a, b), v in zip(section_bounds(t, N_SECTIONS), r["sections"]):
            ax[0].hlines(v, t[a], t[b], color="k", lw=1)
        ax[0].set_ylabel(f"$Nu$ via {QOI[key][1]}")
        if key not in CONVERGENCE_QOI:
            ax[0].set_title("reported only (not a convergence check)", fontsize=9)
        ax[1].plot(t[1::k], cumu[::k])
        ax[1].axhline(r["Nu"], color="C1", lw=0.8)
        ax[1].set_title(f"$Nu = {r['Nu']:.5g} \\pm {r['se_tau']:.2g}$, sections {r['spread_pct']:.2f}%, "
                        f"$p={r['p_sections']:.2f}$")
    for a in axes[-1]:
        if not a.get_xlabel():
            a.set_xlabel(r"$t$")
    fig.savefig(path)
    plt.close(fig)


def _t0_arg(v: str) -> float | str:
    """argparse type: a number or ``"auto"``."""
    return v if v == "auto" else float(v)


def get_parser() -> argparse.ArgumentParser:
    """Command-line options (see the module docstring).

    Returns
    -------
    argparse.ArgumentParser
    """
    p = argparse.ArgumentParser(prog="python -m rbc.analysis.nusselt",
                                description="Nu and Re for one run, with convergence diagnostics.")
    p.add_argument("run", type=Path, help="run directory")
    p.add_argument("--t0", "--time", type=_t0_arg, default=0.0, dest="t0",
                   help="start of the averaging window, or 'auto' for the MSER-5 estimate")
    p.add_argument("--t1", "--end_time", type=float, default=0.0, dest="t1",
                   help="end of the averaging window (0 = end of the data)")
    p.add_argument("--no-plot", action="store_true", help="skip outputs/info.pdf")
    return p


def main(argv=None) -> None:
    """Run :func:`compute` on one run directory and write its outputs.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments (default: ``sys.argv[1:]``).
    """
    a = get_parser().parse_args(argv)
    run = a.run.resolve()
    res, series = compute(run, a.t0, a.t1 if a.t1 > 0 else np.inf)
    out = run / "outputs"
    out.mkdir(exist_ok=True)
    (out / "nusselt.json").write_text(json.dumps(res, indent=1))
    write_info(res, out / "info.txt")
    if not a.no_plot:
        plot(res, series, out / "info.pdf")
    print(f"{run}: Nu = {res['Nu']:.6f} ± {res['Nu_se']:.2g}, converged = {res['converged']}, "
          f"converged_stat = {res['converged_stat']}")


if __name__ == "__main__":
    main()