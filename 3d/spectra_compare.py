"""
Compare time-averaged power spectra across multiple simulations.

Discovers every simulation directory under <basepath> that contains a
snapshots/ folder with *_spectra.npz files (produced by spectra.py), loads
the time-averaged spectra from each, and plots them all on shared axes so
you can see convergence with resolution or differences across Ra/Pr.

Folder structure expected (but any depth works):

    <basepath>/
    ├── ra1e7/
    │   └── pr3e-3/
    │       ├── 144_Gam2_2/
    │       │   └── snapshots/*_spectra.npz
    │       ├── 192_Gam2_2/
    │       └── 256_Gam2_2/
    └── ra1e8/
        └── pr1e-3/
            └── 256_Gam2_2/

Each simulation is labelled by its path relative to <basepath>, so you can
see at a glance which line belongs to which Ra/Pr/resolution combination.

Output:
    <basepath>/outputs/compare_velocity.pdf
    <basepath>/outputs/compare_temperature.pdf

Usage:
    spectra_compare.py <basepath>
"""

import re
import numpy as np
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D

from scipy.integrate import simpson

plt.rcParams.update({"font.size": 20})
plt.ioff()


# ---------------------------------------------------------------------------
# Theory slopes
# ---------------------------------------------------------------------------

_THEORY_VEL = {
    "3d": ([-5/3, -11/5],
           [r"$n^{-5/3}$ (Kolmogorov)", r"$n^{-11/5}$ (Bolgiano)"]),
    "kz": ([-3],
           [r"$n_z^{-3}$ (BL)"]),
    "kx": ([-5/3],
           [r"$n_x^{-5/3}$ (Kolmogorov)"]),
    "xy": ([-5/3, -3],
           [r"$n_\perp^{-5/3}$ (Kolmogorov)", r"$n_\perp^{-3}$ (enstrophy)"]),
}

_THEORY_TEMP = {
    "3d": ([-5/3, -7/5],
           [r"$n^{-5/3}$ (Obukhov-Corrsin)", r"$n^{-7/5}$ (Bolgiano)"]),
    "kz": ([-3],
           [r"$n_z^{-3}$ (BL)"]),
    "kx": ([-5/3, -7/5],
           [r"$n_x^{-5/3}$ (OC)", r"$n_x^{-7/5}$ (BO)"]),
    "xy": ([-5/3, -7/5],
           [r"$n_\perp^{-5/3}$ (OC)", r"$n_\perp^{-7/5}$ (BO)"]),
}

# Panel definitions: (k_key, E_key, title, xlabel, ylabel, theory_key)
# Shared between velocity and temperature (field prefix vel_/temp_ prepended below)
_PANELS = [
    ("k_3d", "E_3d",
     "3D isotropic",
     r"$n = \sqrt{n_x^2+n_y^2+n_z^2}$  [mode]",
     r"$E(n)$",
     "3d"),
    ("k_kz", "E_kz",
     "Vertical marginal",
     r"$n_z$  [mode, max $= N_z/2$]",
     r"$E(n_z)$",
     "kz"),
    ("k_kx", "E_kx",
     "Horizontal marginal ($x$)",
     r"$n_x$  [mode, max $= N_x/2$]",
     r"$E(n_x)$",
     "kx"),
    ("k_xy", "E_xy",
     "2D planar",
     r"$n_\perp = \sqrt{n_x^2+n_y^2}$  [mode]",
     r"$E(n_\perp)$",
     "xy"),
]


# ---------------------------------------------------------------------------
# Helpers shared with spectra_cumulative.py
# ---------------------------------------------------------------------------

def time_average(time, E_mat):
    """
    True time average via Simpson's rule: <E> = ∫ E(t) dt / (t_end - t_0).

    Handles unequal time spacing correctly; falls back to arithmetic mean
    when fewer than 2 time steps are available.
    """
    if len(time) < 2 or (time[-1] - time[0]) <= 0.0:
        return E_mat.mean(axis=0)
    return simpson(E_mat, x=time, axis=0) / (time[-1] - time[0])


def _trim_trailing_drop(k, E, threshold=0.25):
    """
    Remove the Nyquist mode if it shows the characteristic ~0.3-decade drop.

    For a real-valued DFT of even length N, the Nyquist coefficient at
    n = N/2 appears only ONCE (not as a ±pair), so after folding it
    accumulates half the power of its neighbours → a systematic drop of
    log10(2) ≈ 0.301 decades.  We detect this by comparing the last point
    against a slope extrapolation from the two preceding points; if the
    excess drop exceeds `threshold` (default 0.25 dec, safely below 0.30),
    the last point is dropped.
    """
    mask = (k > 0) & (E > 0)
    if mask.sum() < 3:
        return k, E
    k_v = k[mask]
    E_v = E[mask]

    dk_log   = np.log10(k_v[-2]) - np.log10(k_v[-3])
    dE_log   = np.log10(E_v[-2]) - np.log10(E_v[-3])
    slope    = dE_log / dk_log if dk_log != 0 else 0.0
    expected = np.log10(E_v[-2]) + slope * (np.log10(k_v[-1]) - np.log10(k_v[-2]))
    excess   = expected - np.log10(E_v[-1])

    if excess > threshold:
        valid_idx = np.where(mask)[0]
        keep = np.ones(len(k), dtype=bool)
        keep[valid_idx[-1]] = False
        return k[keep], E[keep]
    return k, E


def add_theory_slopes(ax, k_data, E_data, slopes, labels, anchor_frac=0.15):
    """
    Full-range theory slope reference lines, anchored near the inertial range.
    """
    mask = (k_data > 0) & (E_data > 0)
    if not mask.any() or not slopes:
        return
    k_v = k_data[mask]
    E_v = E_data[mask]

    log_k_min = np.log10(k_v[0])
    log_k_max = np.log10(k_v[-1])
    log_k_ref = log_k_min + anchor_frac * (log_k_max - log_k_min)
    k_ref     = 10.0 ** log_k_ref
    E_ref     = 10.0 ** np.interp(log_k_ref, np.log10(k_v), np.log10(E_v))

    k_line = np.logspace(log_k_min, log_k_max, 300)
    colors = plt.cm.Greys(np.linspace(0.4, 0.7, len(slopes)))
    for slope, label, color in zip(slopes, labels, colors):
        ax.loglog(k_line, E_ref * (k_line / k_ref) ** slope,
                  ":", color=color, alpha=0.85, linewidth=1.5, label=label,
                  zorder=1)


def setup_yaxis(ax, ymin, ymax):
    """2-3 labelled major ticks, minor ticks + dashed lines at every decade."""
    ax.set_ylim(10**ymin, 10**ymax)
    span  = ymax - ymin
    ntick = 3 if span > 2 else 2
    log_ticks = np.linspace(ymin, ymax, ntick)
    ax.yaxis.set_major_locator(ticker.FixedLocator(10.0 ** log_ticks))
    labels = [rf"$10^{{{f'{l:.1f}'.rstrip('0').rstrip('.')}}}$" for l in log_ticks]
    ax.yaxis.set_major_formatter(ticker.FixedFormatter(labels))
    int_dec = np.arange(int(np.floor(ymin)), int(np.ceil(ymax)) + 1)
    ax.yaxis.set_minor_locator(ticker.FixedLocator(10.0 ** int_dec))
    ax.yaxis.set_minor_formatter(ticker.NullFormatter())
    for exp in int_dec:
        ax.axhline(10.0**exp, linestyle="--", color="gray", alpha=0.4, linewidth=0.7)


# ---------------------------------------------------------------------------
# Simulation discovery and data loading
# ---------------------------------------------------------------------------

def find_simulation_dirs(basepath):
    """
    Recursively find all directories that contain at least one
    snapshots/*_spectra.npz file.  Returns a list of Path objects, sorted
    so that shorter (shallower) paths come first and alphabetically within
    each depth level.

    This means the list is naturally ordered:
      ra1e7/pr3e-3/144_Gam2_2
      ra1e7/pr3e-3/192_Gam2_2
      ra1e7/pr3e-3/256_Gam2_2
      ra1e8/pr1e-3/256_Gam2_2
    """
    basepath = Path(basepath)
    dirs = set()
    for npz in basepath.rglob("snapshots/*_spectra.npz"):
        dirs.add(npz.parent.parent)   # the simulation root (parent of snapshots/)
    return sorted(dirs, key=lambda p: (len(p.parts), str(p)))


def load_simulation(sim_dir):
    """
    Load and time-average all spectra from one simulation directory.

    Reads every snapshots/*_spectra.npz file (sorted by segment number),
    enforces non-decreasing time (discards accidentally re-run segments),
    and returns the time-averaged spectra as a flat dict plus metadata.

    Returns
    -------
    dict with keys:
      label     : str — path relative to whatever the caller chooses
      t_range   : (t0, t1)
      n_frames  : int
      vel_k_*   : 1D arrays of mode-number x-axes
      temp_k_*  : 1D arrays of mode-number x-axes
      vel_E_*   : 1D time-averaged energy arrays
      temp_E_*  : 1D time-averaged energy arrays
    or None if no valid data found.
    """
    snapshots_dir = Path(sim_dir) / "snapshots"
    fs = sorted(
        snapshots_dir.glob("*_spectra.npz"),
        key=lambda p: int(
            re.search(r"snapshots_s(\d+)_spectra\.npz$", p.name).group(1)
        ),
    )
    if not fs:
        return None

    E_KEYS = ("vel_E_3d", "vel_E_xy", "vel_E_kx", "vel_E_kz",
              "temp_E_3d", "temp_E_xy", "temp_E_kx", "temp_E_kz")
    K_KEYS = ("vel_k_3d", "vel_k_xy", "vel_k_kx", "vel_k_kz",
              "temp_k_3d", "temp_k_xy", "temp_k_kx", "temp_k_kz")

    first  = np.load(fs[0])
    time   = first["time"].copy()
    result = {k: first[k].copy() for k in K_KEYS}
    E_acc  = {k: first[k].copy() for k in E_KEYS}   # (nt, nk)
    first.close()

    for fp in fs[1:]:
        d = np.load(fp)
        t_new = d["time"]
        keep  = t_new > time[-1]
        if keep.any():
            time = np.append(time, t_new[keep])
            for k in E_KEYS:
                E_acc[k] = np.concatenate([E_acc[k], d[k][keep]], axis=0)
        d.close()

    if len(time) == 0:
        return None

    # Time-average each spectrum
    for k in E_KEYS:
        result[k] = time_average(time, E_acc[k])   # now 1D (nk,)

    result["t_range"]  = (float(time[0]), float(time[-1]))
    result["n_frames"] = len(time)
    return result


# ---------------------------------------------------------------------------
# Colour/style assignment
# ---------------------------------------------------------------------------

def _build_style_cycle(n):
    """
    Return n (color, linestyle) pairs.

    Use a perceptually uniform colormap (viridis) for color so that up to ~8
    simulations are distinguishable without ambiguity, then cycle through
    solid/dashed/dotted linestyles for extra separation if n > the colormap's
    comfortable range.
    """
    cmap  = plt.cm.Set1
    colors = [cmap(i / max(n - 1, 1)) for i in range(n)]
    styles = [""]
    return [(colors[i], styles[i % len(styles)]) for i in range(n)]


# ---------------------------------------------------------------------------
# Main plotting function
# ---------------------------------------------------------------------------

def make_comparison_figure(simulations, field, theory_dict, basepath):
    """
    Build one 2×2 comparison figure for `field` ('velocity' or 'temperature').

    Each simulation is one line per panel, labelled by its relative path.
    Theory slope references are drawn in grey so they don't compete visually
    with the data lines.

    Parameters
    ----------
    simulations : list of (label, data_dict) tuples, sorted by resolution
    field       : 'velocity' or 'temperature'
    theory_dict : _THEORY_VEL or _THEORY_TEMP
    basepath    : Path — used to construct the output filename
    """
    prefix = "vel_" if field == "velocity" else "temp_"
    fig, axes = plt.subplots(2, 2, figsize=(18, 14), layout="constrained")
    ax_flat   = [axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]]

    styles = _build_style_cycle(len(simulations))

    for panel_idx, (k_sfx, E_sfx, title, xlabel, ylabel, tkey) in enumerate(_PANELS):
        ax        = ax_flat[panel_idx]
        k_key     = prefix + k_sfx
        E_key     = prefix + E_sfx
        is_1d_marginal = tkey in ("kx", "kz")

        # --- Draw each simulation ---
        all_E = []   # collect for auto y-limits
        for (label, data), (color, lstyle) in zip(simulations, styles):
            k = data[k_key]
            E = data[E_key]
            if is_1d_marginal:
                k, E = _trim_trailing_drop(k, E)
            mask = (k > 0) & (E > 0)
            if not mask.any():
                continue
            ax.loglog(k[mask], E[mask],
                      '.-', color=color,
                      linewidth=1.8, label=label, zorder=2, markersize=6)
            all_E.append(E[mask])

        if not all_E:
            continue

        # --- Theory slopes, anchored to the highest-resolution simulation ---
        # Use the last simulation in the list (assumed finest resolution)
        k_ref = simulations[-1][1][k_key]
        E_ref = simulations[-1][1][E_key]
        if is_1d_marginal:
            k_ref, E_ref = _trim_trailing_drop(k_ref, E_ref)
        slopes, slabels = theory_dict[tkey]
        #add_theory_slopes(ax, k_ref, E_ref, slopes, slabels)

        # --- Auto y-limits ---
        E_all = np.concatenate(all_E)
        E_pos = E_all[E_all > 0]
        if E_pos.size:
            log_max = np.log10(E_pos.max())
            log_min = np.log10(E_pos.min())
            setup_yaxis(ax, log_min - 0.5, log_max + 0.5)

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)

    # --- Legend: one entry per simulation, placed outside the top-right panel ---
    # Theory slopes are added to the legend automatically by add_theory_slopes
    legend_handles = []
    for (label, _), (color, lstyle) in zip(simulations, styles):
        legend_handles.append(
            Line2D([0], [0], color=color, linestyle=lstyle,
                   linewidth=1.8, label=label)
        )
    
    # Add theory entries from the last panel's legend
    extra = [h for h in ax_flat[-1].get_legend_handles_labels()[0]]
    if extra:
        legend_handles += ax_flat[-1].get_legend_handles_labels()[0]
        # Remove the per-panel legends; a single figure legend replaces them
        for ax in ax_flat:
            leg = ax.get_legend()
            if leg:
                leg.remove()
    

    fig.legend(
        handles=legend_handles,
        fontsize=11,
        framealpha=0.85,
        #title=field.capitalize() + " spectra",
    )

    t_ranges = [f"{d['t_range'][0]:.2f}–{d['t_range'][1]:.2f}"
                for _, d in simulations]
    fig.suptitle(
        f"Time-averaged {field} spectra\nResolution comparison\n",
        #f"t ranges: {', '.join(t_ranges)}",
        fontsize=24,
    )

    out = Path(basepath) / "outputs"
    out.mkdir(exist_ok=True)
    savepath = out / f"compare_{field}.pdf"
    fig.savefig(savepath, dpi=300, transparent=True)
    plt.close(fig)
    print(f"  Saved: {savepath}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main(basepath):
    basepath = Path(basepath)
    print(f"Searching for simulations under: {basepath}")

    sim_dirs = find_simulation_dirs(basepath)
    if not sim_dirs:
        raise FileNotFoundError(
            f"No *_spectra.npz files found under {basepath}.\n"
            f"Make sure spectra.py has been run on each simulation first."
        )

    print(f"Found {len(sim_dirs)} simulation(s):")
    simulations = []
    for d in sim_dirs:
        label = str(d.relative_to(basepath))
        print(f"  Loading: {label} … ", end="", flush=True)
        data = load_simulation(d)
        if data is None:
            print("SKIPPED (no valid frames)")
            continue
        data["label"] = label
        simulations.append((label, data))
        print(f"{data['n_frames']} frames, "
              f"t = {data['t_range'][0]:.3f} – {data['t_range'][1]:.3f}")

    if not simulations:
        raise RuntimeError("No simulation data could be loaded.")

    print(f"\nPlotting comparison figures …")
    make_comparison_figure(simulations, "velocity",    _THEORY_VEL,  basepath)
    make_comparison_figure(simulations, "temperature", _THEORY_TEMP, basepath)
    print("Done.")


if __name__ == "__main__":
    from docopt import docopt
    args = docopt(__doc__)
    main(Path(args["<basepath>"]))