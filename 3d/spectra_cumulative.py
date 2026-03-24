"""
Computes the time-averaged power spectra from all snapshot .npz files
produced by spectra.py.

Reads files matching  <basepath>/snapshots/*_spectra.npz
Writes to             <basepath>/outputs/

Output files
------------
  cumulative_spectra_velocity.pdf    — 4-panel time-averaged velocity spectra
  cumulative_spectra_temperature.pdf — 4-panel time-averaged temperature spectra
  spec_range.txt                     — decades of dynamic range for each spectrum

Usage:
    spectra_cumulative.py <basepath>
"""

import numpy as np
from pathlib import Path
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from scipy.integrate import simpson

plt.rcParams.update({"font.size": 14})
plt.ioff()


# ---------------------------------------------------------------------------
# Theory slopes (same rationale as spectra.py — see comments there)
# ---------------------------------------------------------------------------

_THEORY_VEL = {
    "3d":  ([-5/3, -11/5],
            [r"$k^{-5/3}$ (Kolmogorov)", r"$k^{-11/5}$ (Bolgiano)"]),
    "kz":  ([-3],
            [r"$k_z^{-3}$ (BL-dominated)"]),
    "kx":  ([-5/3],
            [r"$k_x^{-5/3}$ (Kolmogorov)"]),
    "xy":  ([-5/3, -3],
            [r"$k_\perp^{-5/3}$ (Kolmogorov)", r"$k_\perp^{-3}$ (enstrophy)"]),
}

_THEORY_TEMP = {
    "3d":  ([-5/3, -7/5],
            [r"$k^{-5/3}$ (Obukhov-Corrsin)", r"$k^{-7/5}$ (Bolgiano)"]),
    "kz":  ([-3],
            [r"$k_z^{-3}$ (BL-dominated)"]),
    "kx":  ([-5/3, -7/5],
            [r"$k_x^{-5/3}$ (Obukhov-Corrsin)", r"$k_x^{-7/5}$ (Bolgiano)"]),
    "xy":  ([-5/3, -7/5],
            [r"$k_\perp^{-5/3}$ (Obukhov-Corrsin)",
             r"$k_\perp^{-7/5}$ (Bolgiano)"]),
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def add_theory_slopes(ax, k_data, E_data, slopes, labels, position_frac=0.45):
    """
    Draw reference power-law slope lines on a log-log axis.

    Each line is anchored at a reference point `position_frac` of the way
    through the data and extends ±0.7 decades in k.  This lets you visually
    judge whether any stretch of the spectrum is consistent with a given
    theoretical prediction.
    """
    mask = E_data > 0
    if not mask.any() or not slopes:
        return
    k_v, E_v = k_data[mask], E_data[mask]
    ref       = int(np.clip(len(k_v) * position_frac, 0, len(k_v) - 1))
    k0, E0    = k_v[ref], E_v[ref]

    colors = plt.cm.Set1(np.linspace(0.05, 0.75, len(slopes)))
    for slope, label, color in zip(slopes, labels, colors):
        k_line = np.logspace(np.log10(k0) - 0.7, np.log10(k0) + 0.7, 60)
        ax.loglog(k_line, E0 * (k_line / k0) ** slope,
                  "--", color=color, alpha=0.75, linewidth=1.8, label=label)
    ax.legend(fontsize=10, framealpha=0.5, loc="lower left")


def setup_yaxis(ax, ymin, ymax):
    """Log y-axis with decade guide lines and clean tick labels."""
    ax.set_ylim(10**ymin, 10**ymax)
    ntick     = max(3, int(ymax - ymin) + 1)
    log_ticks = np.linspace(ymin, ymax, ntick)
    ax.yaxis.set_major_locator(ticker.FixedLocator(10.0**log_ticks))
    ax.yaxis.set_major_formatter(
        ticker.FuncFormatter(lambda v, _: rf"$10^{{{np.log10(v):.1f}}}$")
    )
    for exp in range(int(np.floor(ymin)), int(np.ceil(ymax)) + 1):
        ax.axhline(10.0**exp, linestyle="--", color="gray",
                   alpha=0.35, linewidth=0.8)


def time_average(time, E_mat):
    """
    True time average via Simpson's rule: <E(k)> = ∫ E(k,t) dt / (t_end - t_0).

    time  : shape (nt,)        — simulation times for each snapshot
    E_mat : shape (nt, nk)     — spectral power at each time
    returns: shape (nk,)       — time-averaged spectrum

    Why Simpson and not a simple mean?
    A simple arithmetic mean weights every snapshot equally regardless of
    the time gap between writes.  If your output frequency varies (or if
    you've concatenated segments with different Δt), Simpson's rule correctly
    accounts for the unequal spacing so that the average is a proper
    time-integral rather than a snapshot-count average.
    """
    dt = time[-1] - time[0]
    if dt <= 0.0 or len(time) < 2:
        return E_mat.mean(axis=0)
    return simpson(E_mat, x=time, axis=0) / dt


def count_decades(E_avg):
    """
    Count the number of decades (orders of magnitude) spanned by a spectrum.

    This is the primary diagnostic metric: how many decades of power does
    your simulation resolve?  A larger number means better scale separation
    and a cleaner inertial range.  You want at least 2-3 decades of clean
    power-law behaviour to draw meaningful conclusions.
    """
    mask = E_avg > 0
    if not mask.any():
        return 0.0
    return float(np.log10(E_avg[mask].max()) - np.log10(E_avg[mask].min()))


def plot_panel(ax, k, E_avg, title, xlabel, ylabel, slopes, slope_labels):
    """Render one averaged spectrum panel with theory overlays."""
    mask = (k > 0) & (E_avg > 0)          # skip DC mode on log scale
    ax.loglog(k[mask], E_avg[mask], ".-", linewidth=1.5, markersize=4)
    add_theory_slopes(ax, k[mask], E_avg[mask], slopes, slope_labels)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    # Annotate the number of decades directly on the plot
    decades = count_decades(E_avg)
    ax.text(0.04, 0.06, f"{decades:.2f} decades",
            transform=ax.transAxes, fontsize=11,
            bbox=dict(boxstyle="round", fc="white", alpha=0.75))


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_npz_files(snapshots_dir):
    """
    Load and concatenate all *_spectra.npz files in sorted order.

    Enforces non-decreasing time: if a later file starts before the last
    accepted time (e.g. an accidentally re-run segment), only frames with
    t > t_last are kept.  This prevents the Simpson time-average from being
    corrupted by overlapping or duplicated time windows.

    Returns
    -------
    time   : 1D array (nt,)
    data   : dict mapping key → array
               k-axis keys are 1D (nk,)
               E-value keys are 2D (nt, nk)
    """
    fs = sorted(
        Path(snapshots_dir).glob("*_spectra.npz"),
        key=lambda p: int(
            re.search(r"snapshots_s(\d+)_spectra\.npz$", p.name).group(1)
        ),
    )
    if not fs:
        raise FileNotFoundError(f"No *_spectra.npz files found in {snapshots_dir}")

    # Keys whose arrays should be concatenated along axis=0 (time dimension).
    # k-axis keys are constant and taken only from the first file.
    E_KEYS = ("vel_E_3d", "vel_E_xy", "vel_E_kx", "vel_E_kz",
              "temp_E_3d", "temp_E_xy", "temp_E_kx", "temp_E_kz")
    K_KEYS = ("vel_k_3d", "vel_k_xy", "vel_k_kx", "vel_k_kz",
              "temp_k_3d", "temp_k_xy", "temp_k_kx", "temp_k_kz")

    first    = np.load(fs[0])
    time     = first["time"].copy()
    data     = {k: first[k].copy() for k in K_KEYS}   # 1D k-axes
    E_data   = {k: first[k].copy() for k in E_KEYS}   # 2D (nt, nk) to start
    first.close()

    for fp in fs[1:]:
        d        = np.load(fp)
        new_time = d["time"]

        # Keep only frames strictly after the last accepted time
        keep = new_time > time[-1]
        if not keep.any():
            d.close()
            continue

        time = np.append(time, new_time[keep])
        for k in E_KEYS:
            E_data[k] = np.concatenate([E_data[k], d[k][keep]], axis=0)
        d.close()

    data.update(E_data)
    return time, data


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(basepath):
    basepath  = Path(basepath)
    snapshots = basepath / "snapshots"
    output    = basepath / "outputs"
    output.mkdir(exist_ok=True)

    print("Loading snapshot files …")
    time, data = load_npz_files(snapshots)
    print(f"  {len(time)} time frames, t ∈ [{time[0]:.4f}, {time[-1]:.4f}]")

    # -----------------------------------------------------------------------
    # Define all eight panels (4 velocity + 4 temperature)
    # -----------------------------------------------------------------------
    #
    # Each entry:
    #   (k_key, E_key, title, xlabel, ylabel, theory_key, field_type)
    #
    panels = [
        # --- Velocity ---
        ("vel_k_3d", "vel_E_3d",
         "Velocity - 3D isotropic",
         r"$k$ [rad/L]", r"$E(k)$",
         "3d", "velocity"),

        ("vel_k_kz", "vel_E_kz",
         "Velocity - vertical marginal",
         r"$k_z$ [rad/L]", r"$E(k_z)$",
         "kz", "velocity"),

        ("vel_k_kx", "vel_E_kx",
         "Velocity - horizontal marginal ($x$)",
         r"$k_x$ [rad/L]", r"$E(k_x)$",
         "kx", "velocity"),

        ("vel_k_xy", "vel_E_xy",
         "Velocity - 2D planar",
         r"$k_\perp = \sqrt{k_x^2 + k_y^2}$ [rad/L]", r"$E(k_\perp)$",
         "xy", "velocity"),

        # --- Temperature ---
        ("temp_k_3d", "temp_E_3d",
         "Temperature - 3D isotropic",
         r"$k$ [rad/L]", r"$E_T(k)$",
         "3d", "temperature"),

        ("temp_k_kz", "temp_E_kz",
         "Temperature - vertical marginal",
         r"$k_z$ [rad/L]", r"$E_T(k_z)$",
         "kz", "temperature"),

        ("temp_k_kx", "temp_E_kx",
         "Temperature - horizontal marginal ($x$)",
         r"$k_x$ [rad/L]", r"$E_T(k_x)$",
         "kx", "temperature"),

        ("temp_k_xy", "temp_E_xy",
         "Temperature - 2D planar",
         r"$k_\perp = \sqrt{k_x^2 + k_y^2}$ [rad/L]", r"$E_T(k_\perp)$",
         "xy", "temperature"),
    ]

    # -----------------------------------------------------------------------
    # Two separate 2×2 figures: one for velocity, one for temperature
    # -----------------------------------------------------------------------
    fig_v, axes_v = plt.subplots(2, 2, figsize=(16, 14), layout="constrained")
    fig_t, axes_t = plt.subplots(2, 2, figsize=(16, 14), layout="constrained")

    # Panel layout order: 3D | kz  /  kx | xy   (matches spectra.py)
    axes_v_flat = [axes_v[0, 0], axes_v[0, 1], axes_v[1, 0], axes_v[1, 1]]
    axes_t_flat = [axes_t[0, 0], axes_t[0, 1], axes_t[1, 0], axes_t[1, 1]]

    with open(output / "spec_range.txt", "w") as report:
        report.write(f"Time-averaged spectra: t ∈ [{time[0]:.4f}, {time[-1]:.4f}]\n")
        report.write("=" * 60 + "\n\n")

        for idx, (k_key, E_key, title, xlabel, ylabel, tkey, ftype) in enumerate(panels):
            k_vals = data[k_key]           # 1D (nk,)
            E_mat  = data[E_key]           # 2D (nt, nk)
            E_avg  = time_average(time, E_mat)
            decades = count_decades(E_avg)

            theory  = _THEORY_VEL if ftype == "velocity" else _THEORY_TEMP
            slopes, slope_labels = theory[tkey]

            is_temp = ftype == "temperature"
            panel_idx = idx % 4
            ax = axes_t_flat[panel_idx] if is_temp else axes_v_flat[panel_idx]

            plot_panel(ax, k_vals, E_avg, title, xlabel, ylabel,
                       slopes, slope_labels)

            report.write(f"{title}\n")
            report.write(f"  Decades of dynamic range: {decades:.3f}\n")

            # Optional: auto-set y-limits from the data
            # (comment these two lines out if you prefer to set limits manually)
            mask = (k_vals > 0) & (E_avg > 0)
            if mask.any():
                log_max = np.log10(E_avg[mask].max())
                log_min = np.log10(E_avg[mask].min())
                pad = 0.5
                setup_yaxis(ax, log_min - pad, log_max + pad)

            report.write("\n")

    t_range = f"t ∈ [{time[0]:.3f}, {time[-1]:.3f}]"
    fig_v.suptitle(f"Time-averaged Velocity Spectra\n{t_range}")
    fig_t.suptitle(f"Time-averaged Temperature Spectra (fluctuations T')\n{t_range}")

    vel_path  = output / "cumulative_spectra_velocity.pdf"
    temp_path = output / "cumulative_spectra_temperature.pdf"
    fig_v.savefig(vel_path,  transparent=True, dpi=300)
    fig_t.savefig(temp_path, transparent=True, dpi=300)
    plt.close("all")

    print(f"Saved: {vel_path}")
    print(f"Saved: {temp_path}")
    print(f"Range report: {output / 'spec_range.txt'}")


if __name__ == "__main__":
    from docopt import docopt
    args = docopt(__doc__)
    main(Path(args["<basepath>"]))