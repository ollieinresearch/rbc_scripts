"""Resolution tests of one run: spectral decades per direction, Kolmogorov scale, BL nodes.

Writes ``outputs/resolution.txt``, ``outputs/resolution.pdf`` and
``outputs/resolution.json`` (read by :mod:`rbc.collect`). Usage::

    python -m rbc.analysis.resolution RUN_DIR [--t0 T] [--target 4] [--no-plot] [--no-band]

The averaging window starts where ``outputs/nusselt.json`` says unless ``--t0`` is given.
These are screening diagnostics: none of them measures the error in Nu. Only a
comparison of runs at different resolutions does (:mod:`rbc.analysis.variants`).

Spectral decades, per direction and field
    The x, y and z marginals of velocity and of the temperature perturbation are
    tested separately, so the direction that limits resolution is visible. For each:
    ``D`` of the time-averaged spectrum (:func:`rbc.spectra.core.decades`), its
    standard error by block bootstrap (the spectra in ``N_BLOCKS`` equal-duration blocks
    are resampled with replacement), the mean and 5th percentile of the decades of the
    single writes, and the pile-up. ``status`` is ``pass`` if D reaches the target,
    ``within 2 SE`` if D + 2 SE does, and ``below`` otherwise. The horizontal-mean
    profiles (the mean boundary layers, Chebyshev in z) are reported the same way.
Derivative-weighted decades ``Dk2``
    Decades of n² E(n). ⟨|∇T|²⟩ and ⟨|ω|²⟩ are sums of squared derivatives, so their
    truncation error follows k² E rather than E. In x and y this is exact (a derivative
    multiplies a Fourier coefficient by its wavenumber). In z it is a proxy: with
    ζ = cos θ, T_n = cos(nθ), so n² E(n) is the spectrum of sin θ ∂f/∂ζ.
Kolmogorov and thermal scales
    η = (Pr² / (Ra (Nu - 1)))^(1/4) from the global dissipation rate, in units of H,
    and the thermal scale η Pr^(-3/4) for Pr < 1 (Obukhov–Corrsin) or η Pr^(-1/2) for
    Pr > 1 (Batchelor). ``kmax_eta_<dir>`` = π η / Δ with Δ the largest grid spacing in
    that direction (L/N horizontally, the mid-plane spacing of the Gauss–Chebyshev grid
    vertically). Criteria: ≥ 1 is Δ ≤ π η (Grötzbach 1983); ≥ 1.5 is the stricter
    spectral-DNS criterion (Pope 2000, Yeung et al. 2018). η is a volume average, so the
    local scale in plumes and boundary layers is smaller.
Boundary-layer nodes
    Chebyshev nodes inside the thermal and viscous layers, from ``outputs/profiles.json``.
    This is a vertical-only diagnostic: the horizontal directions are Fourier, with
    uniform spacing, so a node count in a layer says nothing about them; the horizontal
    spectra and Δx / η do. Reported, not tested.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from rbc import cases as rc
from rbc.analysis.nusselt import section_bounds
from rbc.spectra import core as cs

DIRS = ("x", "y", "z")
SPEC_KEY = {"x": "E_kx", "y": "E_ky", "z": "E_nz"}
MODE_KEY = {"x": "kx", "y": "ky", "z": "nz"}
FIELD_NAMES = {"vel": "velocity", "temp": "temperature perturbation"}
#: Default target decades for ``status``.
TARGET = 4.0
#: Pile-up above this many decades is flagged.
PILEUP_FLAG = 0.1
#: Blocks and resamples of the bootstrap standard error of D.
N_BLOCKS, N_BOOT = 10, 400


def block_means(time: np.ndarray, E: np.ndarray, n: int = N_BLOCKS) -> np.ndarray:
    """Time-averaged spectra over ``n`` equal-duration blocks.

    Parameters
    ----------
    time : numpy.ndarray
        Write times, shape (nt,).
    E : numpy.ndarray
        Spectra, shape (nt, n_modes).
    n : int, optional
        Number of blocks.

    Returns
    -------
    numpy.ndarray
        Shape (n, n_modes).
    """
    n = max(2, min(n, len(time) // 2))
    return np.array([cs.time_average(time[a:b + 1], E[a:b + 1]) for a, b in section_bounds(time, n)])


def bootstrap_se(modes: np.ndarray, blocks: np.ndarray, weight=None, n_boot: int = N_BOOT,
                 seed: int = 0) -> float:
    """Standard error of D of the mean spectrum, by resampling blocks.

    Parameters
    ----------
    modes : numpy.ndarray
        Mode numbers.
    blocks : numpy.ndarray
        Block-mean spectra (:func:`block_means`).
    weight : numpy.ndarray, optional
        Multiplies every spectrum before D is taken (e.g. n²).
    n_boot : int, optional
        Resamples.
    seed : int, optional
        Random seed, so the reports are reproducible.

    Returns
    -------
    float
        Standard deviation of D over the resamples.
    """
    rng = np.random.default_rng(seed)
    w = 1.0 if weight is None else weight
    idx = rng.integers(0, len(blocks), size=(n_boot, len(blocks)))
    D = cs.decades(modes, blocks[idx].mean(axis=1) * w)["D"]
    return float(np.nanstd(D, ddof=1))


def status(D: float, se: float, target: float) -> str:
    """Screening verdict of one decades value.

    Parameters
    ----------
    D, se : float
        Decades and its standard error.
    target : float
        Target decades.

    Returns
    -------
    str
        ``"pass"``, ``"within 2 SE"`` or ``"below"``.
    """
    if D >= target:
        return "pass"
    return "within 2 SE" if D + 2 * se >= target else "below"


def spectral_tests(time, modes, data, target=TARGET) -> dict:
    """Decades of every direction and field, and of the mean profiles.

    Parameters
    ----------
    time, modes, data
        As returned by :func:`rbc.spectra.core.load_spectra`, already windowed.
    target : float, optional
        Target decades for ``status``.

    Returns
    -------
    dict
        ``res[field][dir]`` for dir in x, y, z and ``"meanz"``, each with ``n``,
        ``E`` (time mean), ``E_lo``, ``E_hi`` (5th and 95th percentiles over writes),
        ``D_t`` (decades of each write) and the scalars ``D``, ``D_se``, ``D_mean_t``,
        ``D_p05``, ``D_min_t``, ``pileup``, ``Dk2``, ``Dk2_se`` and ``status``.
    """
    out = {}
    for fld in cs.FIELDS:
        res = {}
        series = {d: (modes[MODE_KEY[d]], data[fld][SPEC_KEY[d]]) for d in DIRS}
        series["meanz"] = (modes["nz"], data[fld]["mean_E_nz"])
        for d, (n, E_t) in series.items():
            n = np.asarray(n, float)
            E = cs.time_average(time, E_t)
            dm = cs.decades(n, E)
            dt = cs.decades(n, E_t)["D"]
            blocks = block_means(time, E_t)
            w2 = n ** 2
            r = dict(n=n, E=E, E_lo=np.percentile(E_t, 5, axis=0), E_hi=np.percentile(E_t, 95, axis=0),
                     D_t=dt, D=float(dm["D"]), D_se=bootstrap_se(n, blocks), D_mean_t=float(np.nanmean(dt)),
                     D_p05=float(np.nanpercentile(dt, 5)), D_min_t=float(np.nanmin(dt)),
                     pileup=float(dm["pileup"]), Dk2=float(cs.decades(n, E * w2)["D"]),
                     Dk2_se=bootstrap_se(n, blocks, w2))
            r["status"] = status(r["D"], r["D_se"], target)
            res[d] = r
        out[fld] = res
    return out


def grid_spacings(Nx: int, Ny: int, Nz: int, Lx: float, Ly: float, Lz: float = 1.0) -> dict:
    """Largest grid spacing in each direction.

    Parameters
    ----------
    Nx, Ny, Nz : int
        Modes in each direction.
    Lx, Ly, Lz : float
        Box size.

    Returns
    -------
    dict
        ``x``, ``y``: L / N; ``z``: the largest gap of the Gauss–Chebyshev grid (at mid
        height, about π Lz / (2 Nz)); ``z_wall``: the first node's distance from a wall.
    """
    z = cs.chebyshev_grid(Nz, (0.0, Lz))
    return dict(x=Lx / Nx, y=Ly / Ny, z=float(np.max(np.diff(z))), z_wall=float(z[0]))


def scale_tests(Ra: float, Pr: float, Nu: float, spacing: dict) -> dict:
    """Kolmogorov and thermal dissipation scales against the grid.

    Parameters
    ----------
    Ra, Pr, Nu : float
        Parameters and Nusselt number.
    spacing : dict
        Output of :func:`grid_spacings`.

    Returns
    -------
    dict
        ``eta``, ``eta_T``, and ``kmax_eta_<dir>``, ``kmax_etaT_<dir>`` = π η / Δ for dir
        in x, y, z.
    """
    eta = (Pr ** 2 / (Ra * max(Nu - 1, 1e-300))) ** 0.25
    eta_T = eta * Pr ** (-0.75 if Pr < 1 else -0.5)
    out = dict(eta=eta, eta_T=eta_T)
    for d in DIRS:
        out[f"kmax_eta_{d}"] = np.pi * eta / spacing[d]
        out[f"kmax_etaT_{d}"] = np.pi * eta_T / spacing[d]
    return out


def _read(p: Path) -> dict:
    """A JSON file as a dict ({} if missing)."""
    return json.loads(p.read_text()) if p.exists() else {}


def compute(run: Path, t0: float | None = None, t1: float = 0.0, target: float = TARGET,
            source: str = "auto") -> tuple[dict, dict, tuple]:
    """Every resolution test of one run.

    Parameters
    ----------
    run : pathlib.Path
        Run directory with spectra and ``outputs/nusselt.json``.
    t0 : float, optional
        Start of the window (default: ``nusselt.json``'s, else 0).
    t1 : float, optional
        End of the window (0: end of the data).
    target : float, optional
        Target decades.
    source : {"auto", "sim", "snapshots"}, optional
        Which spectra files to read.

    Returns
    -------
    summary : dict
        Flat scalars for ``resolution.json`` (see :func:`flatten`).
    spec : dict
        Output of :func:`spectral_tests`.
    window : tuple
        ``(time, attrs)`` of the windowed spectra.
    """
    run = Path(run)
    nu = _read(run / "outputs" / "nusselt.json")
    prof = _read(run / "outputs" / "profiles.json")
    if t0 is None:
        t0 = nu.get("t0", 0.0)
    files = cs.spectra_files(run, source)
    time, modes, data, attrs = cs.load_spectra(files)
    from rbc.spectra.cumulative import select_window
    time, data = select_window(time, data, t0, t1)
    spec = spectral_tests(time, modes, data, target)

    p = rc.parse_run_path(run) or {}
    prm = rc.read_params(run)
    Ra = nu.get("Ra", 10.0 ** p.get("log_ra", np.nan))
    Pr = nu.get("Pr", 10.0 ** p.get("log_pr", np.nan))
    Lx = float(attrs.get("Lx", prm.get("LX", 2.0)))
    Ly = float(attrs.get("Ly", prm.get("LY", 2.0)))
    Lz = float(attrs.get("Lz", 1.0))
    Nx, Ny, Nz = (int(attrs.get(k, 0)) for k in ("Nx", "Ny", "Nz"))
    spacing = grid_spacings(Nx, Ny, Nz, Lx, Ly, Lz)
    scales = scale_tests(Ra, Pr, nu["Nu"], spacing) if "Nu" in nu else {}
    summary = flatten(spec, scales, spacing, prof, nu, target)
    summary.update(t0=float(time[0]), t1=float(time[-1]), n_spectra=int(time.size), Nx=Nx, Ny=Ny, Nz=Nz,
                   Lx=Lx, Ly=Ly, Lz=Lz, Ra=float(Ra), Pr=float(Pr), config=rc.run_config(run))
    return summary, spec, (time, attrs)


def flatten(spec: dict, scales: dict, spacing: dict, prof: dict, nu: dict, target: float) -> dict:
    """The scalars of every test, for ``resolution.json`` and :mod:`rbc.collect`.

    Parameters
    ----------
    spec, scales, spacing : dict
        Outputs of :func:`spectral_tests`, :func:`scale_tests`, :func:`grid_spacings`.
    prof, nu : dict
        ``profiles.json`` and ``nusselt.json`` of the run (may be empty).
    target : float
        Target decades.

    Returns
    -------
    dict
        Per field f in {vel, temp} and direction d in {x, y, z, meanz}: ``D_<f>_<d>``,
        ``Dse_<f>_<d>``, ``Dp05_<f>_<d>``, ``Dk2_<f>_<d>``, ``pileup_<f>_<d>``,
        ``status_<f>_<d>``. ``D_run`` (smallest D over fields and x, y, z), ``D_run_at``,
        ``D_run_se``, ``screen`` (status of D_run, ``"pile-up"`` if any direction is
        flagged), the scale tests, ``dx``, ``dy``, ``dz``, ``dz_wall``, the node counts
        and the identity mismatches.
    """
    s = dict(target=target)
    worst = None
    pile = []
    for fld, res in spec.items():
        for d, r in res.items():
            k = f"{fld}_{d}"
            s.update({f"D_{k}": r["D"], f"Dse_{k}": r["D_se"], f"Dp05_{k}": r["D_p05"],
                      f"Dk2_{k}": r["Dk2"], f"pileup_{k}": r["pileup"], f"status_{k}": r["status"]})
            if d in DIRS:
                if worst is None or r["D"] < worst[1]["D"]:
                    worst = (k, r)
                if r["pileup"] > PILEUP_FLAG:
                    pile.append(k)
    s.update(D_run=worst[1]["D"], D_run_at=worst[0], D_run_se=worst[1]["D_se"],
             screen="pile-up" if pile else worst[1]["status"], pileup_in=pile)
    s.update({k: float(v) for k, v in scales.items()})
    s.update({f"d{k}": v for k, v in spacing.items()})
    s.update({k: prof[k] for k in prof if k.startswith("n_delta") or k.startswith("delta")})
    for k in ("qoi_spread_pct", "vort_diff_pct", "Nu", "Nu_se"):
        if k in nu:
            s[k] = nu[k]
    return s


def write_txt(s: dict, path: Path) -> None:
    """Write ``outputs/resolution.txt``.

    Parameters
    ----------
    s : dict
        Summary from :func:`compute`.
    path : pathlib.Path
        File to write.
    """
    c = s["config"]
    L = [f"Ra = {s['Ra']:.4e}, Pr = {s['Pr']:.4e}; N_x = {s['Nx']}, N_y = {s['Ny']}, N_z = {s['Nz']}; "
         f"stepper {c['stepper']}, CFL safety {c['safety']}",
         f"Window: t in [{s['t0']:.3f}, {s['t1']:.3f}], {s['n_spectra']} spectra"]
    L += [f"WARNING: {m}" for m in c["mismatch"]]
    L += ["", "These are screening diagnostics. The error in Nu is measured only by comparing runs",
          "(rbc.analysis.variants, <Ra>/<Pr>/outputs/variants.txt).", "",
          f"1. Spectral decades of the time-averaged marginal spectra (target {s['target']:g})",
          "   D = log10(peak) - mean log10 E over the last 3 modes; SE by block bootstrap;",
          "   p05 = 5th percentile of single writes; Dk2 = decades of n^2 E (derivative-weighted).",
          f"   {'field':<6}{'dir':<7}{'D':>7}{'SE':>7}{'p05':>7}{'Dk2':>7}{'pile-up':>9}  status"]
    for fld in ("vel", "temp"):
        for d in DIRS + ("meanz",):
            k = f"{fld}_{d}"
            flag = "  <- pile-up" if s[f"pileup_{k}"] > PILEUP_FLAG else ""
            note = "   (mean profile, reported only)" if d == "meanz" else ""
            L.append(f"   {fld:<6}{d:<7}{s[f'D_{k}']:7.2f}{s[f'Dse_{k}']:7.2f}{s[f'Dp05_{k}']:7.2f}"
                     f"{s[f'Dk2_{k}']:7.2f}{s[f'pileup_{k}']:9.2f}  {s[f'status_{k}']}{flag}{note}")
    L.append(f"   Limiting: {s['D_run_at']} with D = {s['D_run']:.2f} ± {s['D_run_se']:.2f}; screen: {s['screen']}")
    if "eta" in s:
        L += ["", "2. Dissipation scales against the largest grid spacing (pi*eta/Delta; >= 1 Grotzbach, "
                  ">= 1.5 strict)",
              f"   eta = {s['eta']:.3e}, thermal scale = {s['eta_T']:.3e} (volume averages; smaller in plumes "
              "and boundary layers)",
              f"   Delta: x {s['dx']:.3e}, y {s['dy']:.3e}, z (mid-plane) {s['dz']:.3e}, first node "
              f"{s['dz_wall']:.3e}"]
        for d in DIRS:
            ke, kt = s[f"kmax_eta_{d}"], s[f"kmax_etaT_{d}"]
            verdict = "strict" if ke >= 1.5 else "Grotzbach" if ke >= 1 else "BELOW"
            L.append(f"   {d}: kmax*eta = {ke:6.2f} ({verdict}),  kmax*eta_T = {kt:6.2f}")
    nodes = {k: v for k, v in s.items() if k.startswith("n_delta")}
    if nodes:
        L += ["", "3. Chebyshev nodes inside the boundary layers (vertical only; reported, not tested)"]
        L += [f"   {k[2:]}: {v}" for k, v in nodes.items()]
    if "qoi_spread_pct" in s:
        L += ["", "4. Exact identities (reported; in a dealiased spectral code they hold for the truncated",
              "   system, so agreement is weak evidence of resolution)",
              f"   <wT> vs <|grad T|^2>: {s['qoi_spread_pct']:.3f}% of Nu-1"]
        if "vort_diff_pct" in s:
            L.append(f"   vorticity vs <wT>: {s['vort_diff_pct']:+.3f}% of Nu-1")
    path.write_text("\n".join(L) + "\n")


def plot(s: dict, spec: dict, window: tuple, path: Path, band: bool = True) -> None:
    """Two-page ``outputs/resolution.pdf``.

    Page 1: the time-averaged x, y and z spectra of each field, with n² E (scaled to the
    same peak) and, in z, the mean-profile spectrum. Page 2: decades with ±2 SE against
    the target, the decades of single writes over time, and π η / Δ per direction.

    Parameters
    ----------
    s : dict
        Summary from :func:`compute`.
    spec : dict
        Output of :func:`spectral_tests`.
    window : tuple
        ``(time, attrs)`` from :func:`compute`.
    path : pathlib.Path
        File to write.
    band : bool, optional
        Shade the 5-95% range of single writes.
    """
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from rbc import style

    style.use("diag")
    time = window[0]
    with PdfPages(path) as pdf:
        fig, axes = plt.subplots(2, 3, figsize=(16, 9.5), layout="constrained")
        for i, fld in enumerate(cs.FIELDS):
            for j, d in enumerate(DIRS):
                ax = axes[i, j]
                r = spec[fld][d]
                n, E = r["n"], r["E"]
                m = (n >= 1) & (E > 0)
                if band:
                    ax.fill_between(n[m], r["E_lo"][m], r["E_hi"][m], alpha=0.2, lw=0, label="5-95% of writes")
                ax.loglog(n[m], E[m], ".-", ms=3, lw=1.2, label="$E$")
                E2 = n ** 2 * E
                ax.loglog(n[m], E2[m] * E[m].max() / E2[m].max(), "--", lw=1, label=r"$n^2E$ (scaled)")
                arrays = [E[m]]
                if d == "z":
                    q = spec[fld]["meanz"]
                    mm = (q["n"] >= 1) & (q["E"] > 0)
                    ax.loglog(q["n"][mm], q["E"][mm], ":", lw=1, label="horizontal mean")
                    arrays.append(q["E"][mm])
                ax.text(0.97, 0.97, f"D = {r['D']:.2f} ± {r['D_se']:.2f}\nDk2 = {r['Dk2']:.2f}\n"
                        f"pile-up {r['pileup']:.2f}\n{r['status']}", transform=ax.transAxes, ha="right",
                        va="top", fontsize=9, bbox=dict(boxstyle="round", fc="white", alpha=0.8))
                cs.setup_yaxis(ax, *cs.auto_ylim(*arrays))
                ax.set_title(f"{FIELD_NAMES[fld]}, {d}" + (" (Chebyshev degree)" if d == "z" else ""))
                ax.set_xlabel("$n_z$" if d == "z" else f"$k_{d}$")
                ax.legend(fontsize=8, loc="lower left")
        fig.suptitle(f"Time-averaged marginal spectra, t in [{s['t0']:.1f}, {s['t1']:.1f}]")
        pdf.savefig(fig)
        plt.close(fig)

        fig, axes = plt.subplots(1, 3, figsize=(16, 4.8), layout="constrained")
        ax = axes[0]
        labels, k = [], 0
        for fld in cs.FIELDS:
            for d in DIRS + ("meanz",):
                r = spec[fld][d]
                ax.errorbar(k, r["D"], yerr=2 * r["D_se"], fmt="o", capsize=3,
                            color="C0" if fld == "vel" else "C2", mfc="none" if d == "meanz" else None)
                ax.plot(k, r["Dk2"], "x", color="0.4")
                labels.append(f"{fld}\n{d}")
                k += 1
        ax.axhline(s["target"], color="C3", ls="--", lw=1)
        ax.set_xticks(range(k), labels, fontsize=7)
        ax.set_ylabel("decades (±2 SE); x: Dk2")
        ax.set_title(f"screen: {s['screen']}")
        ax = axes[1]
        for fld, c in zip(cs.FIELDS, ("C0", "C2")):
            for d, ls in zip(DIRS, ("-", "--", ":")):
                ax.plot(time, spec[fld][d]["D_t"], ls=ls, color=c, lw=0.7, label=f"{fld} {d}")
        ax.axhline(s["target"], color="C3", ls="--", lw=1)
        ax.set_xlabel("$t$")
        ax.set_ylabel("decades of single writes")
        ax.legend(fontsize=7, ncol=2)
        ax = axes[2]
        if "eta" in s:
            x = np.arange(len(DIRS))
            ax.bar(x - 0.18, [s[f"kmax_eta_{d}"] for d in DIRS], 0.36, label=r"$\pi\eta/\Delta$")
            ax.bar(x + 0.18, [s[f"kmax_etaT_{d}"] for d in DIRS], 0.36, label=r"$\pi\eta_T/\Delta$")
            ax.axhline(1, color="C3", ls="--", lw=1, label="Grötzbach")
            ax.axhline(1.5, color="C3", ls=":", lw=1, label="strict")
            ax.set_xticks(x, ["x", "y", "z (mid)"])
            ax.set_ylabel(r"$k_{\max}\eta$")
            ax.legend(fontsize=8)
        fig.suptitle(f"N = ({s['Nx']}, {s['Ny']}, {s['Nz']}), {s['config']['stepper']}, "
                     f"CFL safety {s['config']['safety']}")
        pdf.savefig(fig)
        plt.close(fig)


def run(run_dir, t0=None, t1=0.0, target=TARGET, source="auto", plots=True, band=True) -> dict:
    """Compute and write the resolution tests of one run.

    Parameters
    ----------
    run_dir : str or pathlib.Path
        Run directory.
    t0, t1 : float, optional
        Averaging window (default start: ``nusselt.json``'s).
    target : float, optional
        Target decades.
    source : {"auto", "sim", "snapshots"}, optional
        Spectra files to read.
    plots : bool, optional
        Write ``resolution.pdf``.
    band : bool, optional
        Shade the 5-95% range in the PDF.

    Returns
    -------
    dict
        The summary written to ``resolution.json``.
    """
    run_dir = Path(run_dir).resolve()
    s, spec, window = compute(run_dir, t0, t1, target, source)
    out = run_dir / "outputs"
    out.mkdir(exist_ok=True)
    (out / "resolution.json").write_text(json.dumps(s, indent=1, default=float))
    write_txt(s, out / "resolution.txt")
    if plots:
        plot(s, spec, window, out / "resolution.pdf", band)
    return s


def get_parser() -> argparse.ArgumentParser:
    """Command-line options (see the module docstring).

    Returns
    -------
    argparse.ArgumentParser
    """
    p = argparse.ArgumentParser(prog="python -m rbc.analysis.resolution", description=__doc__.split("\n")[0])
    p.add_argument("run", type=Path, help="run directory")
    p.add_argument("--t0", "--time", type=float, default=None, dest="t0",
                   help="start of the averaging window (default: the one in outputs/nusselt.json)")
    p.add_argument("--t1", "--end_time", type=float, default=0.0, dest="t1", help="end of the window (0 = end)")
    p.add_argument("--target", type=float, default=TARGET, help="target decades for the status column")
    p.add_argument("--source", choices=("auto", "sim", "snapshots"), default="auto", help="spectra files to read")
    p.add_argument("--no-plot", action="store_true", help="skip outputs/resolution.pdf")
    p.add_argument("--no-band", action="store_true", help="do not shade the 5-95%% range of single writes")
    return p


def main(argv=None) -> None:
    """Command-line entry point.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments (default: ``sys.argv[1:]``).
    """
    a = get_parser().parse_args(argv)
    s = run(a.run, a.t0, a.t1, a.target, a.source, not a.no_plot, not a.no_band)
    print(f"{a.run}: limiting {s['D_run_at']} D = {s['D_run']:.2f} ± {s['D_run_se']:.2f} ({s['screen']})")


if __name__ == "__main__":
    main()
