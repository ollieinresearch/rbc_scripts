"""
Script to calculate power spectra, plot them at each time, and save them for cumulative time averaging.

Physical setup
--------------
Domain : Lx = Ly = 2 (periodic), Lz = 1 (Chebyshev / wall-bounded)
Bases  : Fourier in x, y  |  Chebyshev type-2 (Gauss-Lobatto) in z

Usage:
    spectra.py <files>... [--vmins=<vmins>] [--vmaxs=<vmaxs>]
                          [--tmins=<tmins>] [--tmaxs=<tmaxs>]

Options:
    --vmins=<vmins>   Comma-separated log10 y-axis minima for velocity (4 values).
    --vmaxs=<vmaxs>   Comma-separated log10 y-axis maxima for velocity (4 values).
    --tmins=<tmins>   Comma-separated log10 y-axis minima for temperature (4 values).
    --tmaxs=<tmaxs>   Comma-separated log10 y-axis maxima for temperature (4 values).

Panel order (both velocity and temperature figures):
    [0] 3D isotropic E(k)
    [1] Vertical marginal E(kz)
    [2] Horizontal marginal E(kx)
    [3] 2D planar E(k_perp)
"""

import h5py as h5
import numpy as np
from numpy.polynomial.chebyshev import chebpts2
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from scipy.interpolate import interp1d

# ---------------------------------------------------------------------------
# Module-level constants — adjust if your domain dimensions differ
# ---------------------------------------------------------------------------
LX, LY, LZ = 2.0, 2.0, 1.0

s = 30
plt.rcParams.update({"font.size": 0.75 * s})
plt.ioff()


# ===========================================================================
#  Theoretical slope reference lines
# ===========================================================================

# Velocity spectra scaling laws:
# ┌─────────────────────────────────────────────────────────────────────────┐
# │  Kolmogorov (K41): E_v(k) ~ k^{-5/3}                                    │
# │    The classical inertial-range prediction for isotropic turbulence.    │
# │    Valid when the kinetic energy cascade rate ε dominates buoyancy.     │
# │                                                                         │
# │  Bolgiano-Obukhov (BO59): E_v(k) ~ k^{-11/5}                            │
# │    Applies when buoyancy forces are dominant at a given scale (the      │
# │    "Bolgiano scale" L_B = (N^3/ε_T)^{1/2} is in the inertial range).    │
# │    In most laboratory/numerical RBC, L_B ≈ L (whole domain), so BO      │
# │    scaling is rarely cleanly observed for velocity.                     │
# └─────────────────────────────────────────────────────────────────────────┘
#
# Temperature spectra scaling laws:
# ┌─────────────────────────────────────────────────────────────────────────┐
# │  Obukhov-Corrsin (K41 passive scalar): E_T(k) ~ k^{-5/3}                │
# │    Valid when T is advected passively and the scalar dissipation rate   │
# │    χ is the relevant quantity. k^{-5/3} for the scalar in K41 regime.   │
# │                                                                         │
# │  Bolgiano-Obukhov: E_T(k) ~ k^{-7/5}                                    │
# │    The counterpart BO slope for the temperature spectrum. Note it is    │
# │    shallower than the K41 -5/3, which seems counter-intuitive but       │
# │    follows from the BO dimensional argument.                            │
# └─────────────────────────────────────────────────────────────────────────┘
#
# Vertical marginal slope:
# ┌─────────────────────────────────────────────────────────────────────────┐
# │  E(kz) ~ kz^{-3}                                                        │
# │    Steep rolloff expected from thin boundary-layer structure. Because   │
# │    the BL thickness δ ~ Ra^{-1/3}, energy at scales < δ decays rapidly. │
# │    This is not a universal law but a useful empirical reference.        │
# └─────────────────────────────────────────────────────────────────────────┘
#
# 2D planar slope (horizontal):
# ┌─────────────────────────────────────────────────────────────────────────┐
# │  E(k_perp) ~ k_perp^{-5/3}  (Kolmogorov)                                │
# │  E(k_perp) ~ k_perp^{-3}    (2D enstrophy cascade / quasi-2D regime)    │
# │    At large Ra the flow can show a forward enstrophy cascade for        │
# │    scales larger than the convective cell size.                         │
# └─────────────────────────────────────────────────────────────────────────┘

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


# ===========================================================================
#  Core analysis class
# ===========================================================================

class RBCSpectraAnalyzer:
    """
    Computes power spectra (3D isotropic, 2D planar, 1D marginals) for
    velocity and temperature fields from a Dedalus RBC simulation.

    Grid assumptions
    ----------------
    x, y : uniform Fourier grids on [0, Lx) and [0, Ly)
    z    : Chebyshev type-2 (Gauss-Lobatto interior) on [-Lz/2, Lz/2]

    Wavenumber convention
    ---------------------
    All wavenumbers are ANGULAR (rad / unit length): k = 2π x (cycles / L).
    This makes kx_max = kz_max = π x N/L for a matched grid, so comparisons
    across directions are on an equal footing.
    """

    def __init__(self, Lx=2.0, Ly=2.0, Lz=1.0):
        self.Lx = Lx
        self.Ly = Ly
        self.Lz = Lz
        self._grid = None   # cached grid dict; rebuilt when resolution changes

    # ------------------------------------------------------------------
    # Grid construction
    # ------------------------------------------------------------------

    def _build_grid(self, nx, ny, nz):
        """
        Build (and cache) all grid-space and wavenumber arrays.

        Returns a dict with:
          z_cheb, z_uniform, sort_idx     — z grids
          kx, ky, kz                      — 1D angular wavenumber arrays
          Knorm                           — flattened 3D |k| array (nx*ny*nz,)
          Kperp                           — flattened 2D k_perp array (nx*ny,)
          bins_3d, k_centers_3d          — isotropic histogram bins / centres
          bins_xy, k_centers_xy          — planar histogram bins / centres
        """
        if self._grid is not None and self._grid["shape"] == (nx, ny, nz):
            return self._grid

        # --- z grid (Chebyshev) ---
        # chebpts2(n) returns cos(πj/(n-1)) for j=0…n-1 → DESCENDING [1, …, -1].
        # Scale to [-Lz/2, Lz/2] and sort ascending for interpolation.
        z_raw    = chebpts2(nz) * (self.Lz / 2.0)   # descending
        sort_idx = np.argsort(z_raw)                  # ascending permutation
        z_cheb   = z_raw[sort_idx]                    # now ascending
        z_uniform = np.linspace(-self.Lz / 2.0, self.Lz / 2.0, nz, endpoint=False)

        # --- Angular wavenumbers ---
        # k = 2π × fftfreq / (dx) = 2π × fftfreq × (N/L)
        # Equivalently: np.fft.fftfreq(N, d=L/N) gives cycles/unit,
        # multiplying by 2π converts to angular (rad/unit).
        kx = 2 * np.pi * np.fft.fftfreq(nx, d=self.Lx / nx)  # (nx,)
        ky = 2 * np.pi * np.fft.fftfreq(ny, d=self.Ly / ny)  # (ny,)
        kz = 2 * np.pi * np.fft.fftfreq(nz, d=self.Lz / nz)  # (nz,)

        # --- 3D isotropic binning ---
        KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing="ij")
        Knorm = np.sqrt(KX**2 + KY**2 + KZ**2).ravel()  # (nx*ny*nz,)

        # Bin width = smallest angular wavenumber step among all directions.
        # δk_x = 2π/Lx = π,  δk_y = π,  δk_z = 2π/Lz = 2π  → min = π.
        # Using this width ensures no mode is ever split across two bins.
        dk = 2.0 * np.pi / max(self.Lx, self.Ly, self.Lz)   # = π here
        Kmax       = Knorm.max()
        bins_3d    = np.arange(0.5 * dk, Kmax + dk, dk)
        k_centers_3d = 0.5 * (bins_3d[:-1] + bins_3d[1:])

        # --- 2D planar binning ---
        KX2D, KY2D = np.meshgrid(kx, ky, indexing="ij")
        Kperp       = np.sqrt(KX2D**2 + KY2D**2).ravel()  # (nx*ny,)
        Kperp_max   = Kperp.max()
        bins_xy     = np.arange(0.5 * dk, Kperp_max + dk, dk)
        k_centers_xy = 0.5 * (bins_xy[:-1] + bins_xy[1:])

        self._grid = dict(
            shape=(nx, ny, nz),
            z_cheb=z_cheb, z_uniform=z_uniform, sort_idx=sort_idx,
            kx=kx, ky=ky, kz=kz,
            Knorm=Knorm,
            Kperp=Kperp,
            bins_3d=bins_3d, k_centers_3d=k_centers_3d,
            bins_xy=bins_xy, k_centers_xy=k_centers_xy,
            dk=dk,
        )
        return self._grid

    # ------------------------------------------------------------------
    # z-interpolation: Chebyshev → uniform
    # ------------------------------------------------------------------

    def _regrid_z(self, field, g):
        """
        Interpolate a (nx, ny, nz) field from the (ascending) Chebyshev z-grid
        to the uniform z-grid using a cubic spline along the z-axis.

        Why cubic?  Near the walls the Chebyshev points are dense and the
        boundary layer has high curvature; linear interpolation introduces
        visible spectral pollution at high kz.

        Memory note: interp1d with axis=-1 is vectorized over all (nx, ny)
        pairs simultaneously without building an O(N^3) coordinate array.
        """
        # Reorder data z-axis to match the ascending z_cheb
        field_asc = field[:, :, g["sort_idx"]]

        f_interp = interp1d(
            g["z_cheb"], field_asc,
            axis=-1, kind="cubic",
            bounds_error=False,
            fill_value="extrapolate",
        )
        return f_interp(g["z_uniform"])

    # ------------------------------------------------------------------
    # Core spectral machinery
    # ------------------------------------------------------------------

    @staticmethod
    def _one_sided_marginal(E_3d, k_1d, sum_axes):
        """
        Compute a one-sided 1D marginal power spectrum.

        Steps
        -----
        1. Sum E_3d over `sum_axes` (the two directions being marginalised).
        2. Map every wavenumber index to |k|.
        3. Fold ±k together by SUMMING their contributions.
           For a real-valued field, |X(-k)|² = |X(+k)|² exactly, so this
           doubles the power at every non-zero k compared with looking at
           positive frequencies only.  The DC component (k=0) and the Nyquist
           component each appear once and are not doubled.
        4. Return (k_unique_positive, E_folded), both ascending in k.

        Integrating E_folded x dk over k ≥ 0 gives the total energy in
        that marginal direction (up to the field-normalisation convention).
        """
        E_raw     = E_3d.sum(axis=sum_axes)              # 1D, indexed by k_1d
        k_abs     = np.abs(k_1d)
        k_unique, inv = np.unique(k_abs, return_inverse=True)
        E_folded  = np.bincount(inv, weights=E_raw, minlength=len(k_unique))
        return k_unique, E_folded

    def _bin_isotropic(self, E_3d, g):
        """Bin a 3D power array into the isotropic 1D spectrum."""
        pk, _   = np.histogram(g["Knorm"], bins=g["bins_3d"], weights=E_3d.ravel())
        nm, _   = np.histogram(g["Knorm"], bins=g["bins_3d"])
        mask    = nm > 0
        E_iso   = np.where(mask, pk / np.where(mask, nm, 1.0), 0.0)
        return g["k_centers_3d"], E_iso

    def _bin_planar(self, E_3d, g):
        """Sum over kz, then bin the (nx×ny) 2D power by k_perp."""
        E_z_sum   = E_3d.sum(axis=2).ravel()
        pk_xy, _  = np.histogram(g["Kperp"], bins=g["bins_xy"], weights=E_z_sum)
        nm_xy, _  = np.histogram(g["Kperp"], bins=g["bins_xy"])
        mask      = nm_xy > 0
        E_xy      = np.where(mask, pk_xy / np.where(mask, nm_xy, 1.0), 0.0)
        return g["k_centers_xy"], E_xy

    # ------------------------------------------------------------------
    # Public: velocity spectra
    # ------------------------------------------------------------------

    def compute_velocity_spectra(self, u, v, w, g):
        """
        Compute kinetic-energy spectra: E = |U|² + |V|² + |W|².

        All three components are interpolated, transformed, and their power
        arrays accumulated before any binning.  This is more efficient than
        binning each component separately (3 FFTs, 1 binning pass).

        Returns dict with keys: k_3d, E_3d, k_xy, E_xy, k_kx, E_kx, k_kz, E_kz,
                                 parseval_grid, parseval_spec.
        """
        parseval_grid = 0.0
        E_total       = None

        for comp in (u, v, w):
            comp_u = self._regrid_z(comp, g)
            parseval_grid += float(np.sum(comp_u**2))
            F   = np.fft.fftn(comp_u, norm="ortho")
            E_c = np.abs(F) ** 2
            del comp_u, F
            E_total = E_c if E_total is None else E_total + E_c
            del E_c

        parseval_spec = float(np.sum(E_total))

        k_3d, E_iso      = self._bin_isotropic(E_total, g)
        k_xy, E_xy       = self._bin_planar(E_total, g)
        k_kx, E_kx       = self._one_sided_marginal(E_total, g["kx"], (1, 2))
        k_kz, E_kz       = self._one_sided_marginal(E_total, g["kz"], (0, 1))
        del E_total

        return dict(
            k_3d=k_3d, E_3d=E_iso,
            k_xy=k_xy, E_xy=E_xy,
            k_kx=k_kx, E_kx=E_kx,
            k_kz=k_kz, E_kz=E_kz,
            parseval_grid=parseval_grid,
            parseval_spec=parseval_spec,
        )

    # ------------------------------------------------------------------
    # Public: scalar (temperature) spectra
    # ------------------------------------------------------------------

    def compute_scalar_spectra(self, field, g):
        """
        Compute power spectra for a scalar field (temperature).

        IMPORTANT — mean subtraction
        ─────────────────────────────
        Before taking the FFT we subtract the horizontal (x,y) mean at each
        z-level: T'(x,y,z) = T(x,y,z) - <T>(z).

        Why?  In RBC the temperature has a mean stratification T̄(z) that
        varies smoothly in z and carries almost all its energy in the lowest
        k_z modes.  If you FFT without subtracting T̄(z), the vertical
        spectrum is completely dominated by the background gradient rather
        than the turbulent fluctuations, making comparisons across directions
        meaningless and the inertial-range slopes unreadable.

        Parseval is checked on the fluctuation field T', not the full T.

        Returns same dict keys as compute_velocity_spectra.
        """
        field_u = self._regrid_z(field, g)

        # Subtract horizontal mean at each z to isolate turbulent fluctuations
        T_mean  = field_u.mean(axis=(0, 1), keepdims=True)   # shape (1, 1, nz)
        field_u = field_u - T_mean

        # Parseval on the fluctuation field (not the full T)
        parseval_grid = float(np.sum(field_u**2))

        F   = np.fft.fftn(field_u, norm="ortho")
        E   = np.abs(F) ** 2
        del field_u, F

        parseval_spec = float(np.sum(E))

        k_3d, E_iso  = self._bin_isotropic(E, g)
        k_xy, E_xy   = self._bin_planar(E, g)
        k_kx, E_kx   = self._one_sided_marginal(E, g["kx"], (1, 2))
        k_kz, E_kz   = self._one_sided_marginal(E, g["kz"], (0, 1))
        del E

        return dict(
            k_3d=k_3d, E_3d=E_iso,
            k_xy=k_xy, E_xy=E_xy,
            k_kx=k_kx, E_kx=E_kx,
            k_kz=k_kz, E_kz=E_kz,
            parseval_grid=parseval_grid,
            parseval_spec=parseval_spec,
        )

    # ------------------------------------------------------------------
    # Plotting helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _add_theory_slopes(ax, k_data, E_data, slopes, labels,
                           position_frac=0.45):
        """
        Draw reference power-law lines anchored at position_frac of the way
        through the plotted k range.  Each line spans ±0.7 decades in k.
        """
        mask = E_data > 0
        if not mask.any() or len(slopes) == 0:
            return
        k_v = k_data[mask]
        E_v = E_data[mask]
        ref  = int(np.clip(len(k_v) * position_frac, 0, len(k_v) - 1))
        k0, E0 = k_v[ref], E_v[ref]

        colors = plt.cm.Set1(np.linspace(0.05, 0.75, len(slopes)))
        for slope, label, color in zip(slopes, labels, colors):
            k_line = np.logspace(np.log10(k0) - 0.7, np.log10(k0) + 0.7, 60)
            ax.loglog(k_line, E0 * (k_line / k0) ** slope,
                      "--", color=color, alpha=0.75, linewidth=1.8, label=label)
        ax.legend(fontsize=max(8, plt.rcParams["font.size"] * 0.45),
                  framealpha=0.5, loc="lower left")

    @staticmethod
    def _setup_yaxis(ax, ymin, ymax):
        """Consistent log y-axis with grey decade guide lines."""
        ax.set_ylim(10**ymin, 10**ymax)
        ntick     = max(3, int(ymax - ymin) + 1)
        log_ticks = np.linspace(ymin, ymax, ntick)
        ax.yaxis.set_major_locator(ticker.FixedLocator(10.0 ** log_ticks))
        ax.yaxis.set_major_formatter(
            ticker.FuncFormatter(lambda v, _: rf"$10^{{{np.log10(v):.1f}}}$")
        )
        for exp in range(int(np.floor(ymin)), int(np.ceil(ymax)) + 1):
            ax.axhline(10.0**exp, linestyle="--", color="gray",
                       alpha=0.35, linewidth=0.8)

    def plot_spectra(self, spec, title, savepath, mins, maxs, field_type="velocity"):
        """
        Four-panel figure:  3D isotropic | vertical marginal
                            horizontal marginal | 2D planar

        Parameters
        ----------
        spec      : dict returned by compute_velocity_spectra / compute_scalar_spectra
        title     : figure suptitle string
        savepath  : Path or str
        mins/maxs : list of 4 log10 y-axis limits [3D, kz, kx, xy]
        field_type: "velocity" or "temperature"  — selects theory slopes
        """
        theory = _THEORY_VEL if field_type == "velocity" else _THEORY_TEMP

        fig, axes = plt.subplots(2, 2, figsize=(16, 14), layout="constrained")

        panels = [
            # (ax, k_key, E_key, xlabel, ylabel, title_str, theory_key)
            (axes[0, 0], "k_3d", "E_3d",
             r"$k$ [rad/L]", r"$E(k)$",
             "3D isotropic spectrum", "3d"),

            (axes[0, 1], "k_kz", "E_kz",
             r"$k_z$ [rad/L]", r"$E(k_z)$",
             "Vertical marginal", "kz"),

            (axes[1, 0], "k_kx", "E_kx",
             r"$k_x$ [rad/L]", r"$E(k_x)$",
             "Horizontal marginal ($x$)", "kx"),

            (axes[1, 1], "k_xy", "E_xy",
             r"$k_\perp = \sqrt{k_x^2+k_y^2}$ [rad/L]",
             r"$E(k_\perp)$",
             "2D planar spectrum", "xy"),
        ]

        for idx, (ax, kk, Ek, xl, yl, ttl, tkey) in enumerate(panels):
            k = spec[kk]
            E = spec[Ek]
            # Skip DC / zero mode on log scale
            mask = (k > 0) & (E > 0)
            ax.loglog(k[mask], E[mask], ".-", linewidth=1.2, markersize=4)
            slopes, labels = theory[tkey]
            self._add_theory_slopes(ax, k[mask], E[mask], slopes, labels)
            ax.set_xlabel(xl)
            ax.set_ylabel(yl)
            ax.set_title(ttl)
            self._setup_yaxis(ax, mins[idx], maxs[idx])

        parseval_ok = np.isclose(
            spec["parseval_grid"], spec["parseval_spec"], rtol=1e-2
        )
        fig.suptitle(
            f"{title}\n"
            f"Parseval: {parseval_ok}"
            f"  (grid {spec['parseval_grid']:.4e},"
            f"  spec {spec['parseval_spec']:.4e})"
        )
        fig.savefig(savepath, dpi=600, transparent=True)
        plt.close(fig)


# ===========================================================================
#  Entry point called by post.visit_writes
# ===========================================================================

def main(file, start, count, vmins, vmaxs, tmins, tmaxs):
    fp       = Path(file)
    basepath = Path(fp.parents[1])

    (basepath / "res_check_3d").mkdir(exist_ok=True)
    (basepath / "res_check_temp").mkdir(exist_ok=True)

    analyzer   = RBCSpectraAnalyzer(Lx=LX, Ly=LY, Lz=LZ)
    prev_shape = (0, 0, 0)

    with h5.File(fp, "r") as f:
        times = []

        # Storage: k-axes are identical for every time step at fixed resolution.
        # We keep lists of 1D E-arrays and convert at the end.
        vel  = {k: [] for k in ("E_3d", "E_xy", "E_kx", "E_kz")}
        temp = {k: [] for k in ("E_3d", "E_xy", "E_kx", "E_kz")}
        k_axes_vel  = {}
        k_axes_temp = {}

        nt = f["scales/sim_time"].shape[0]

        for ti in range(nt):
            time  = float(f["scales/sim_time"][ti])
            write = int(f["scales/write_number"][ti])

            # Skip the initial condition (t=0); spectra are trivial there
            if np.isclose(time, 0.0):
                continue

            # --- Load fields ---
            raw_temp = np.array(f["tasks"]["temperature"][ti])

            # Velocity may be stored in several different ways depending on
            # how the simulation was set up.
            try:
                u = np.array(f["tasks"]["u"][ti, 0])
                v = np.array(f["tasks"]["u"][ti, 1])
                w = np.array(f["tasks"]["u"][ti, 2])
            except Exception:
                try:
                    u = np.array(f["tasks"]["x"][ti])
                    v = np.array(f["tasks"]["y"][ti])
                    w = np.array(f["tasks"]["w"][ti])
                except Exception:
                    u = np.array(f["tasks"]["u"][ti])
                    v = np.array(f["tasks"]["v"][ti])
                    w = np.array(f["tasks"]["w"][ti])

            nx, ny, nz = u.shape

            # Rebuild grid only when resolution changes
            if (nx, ny, nz) != prev_shape:
                g          = analyzer._build_grid(nx, ny, nz)
                prev_shape = (nx, ny, nz)

            # --- Velocity spectra ---
            vspec = analyzer.compute_velocity_spectra(u, v, w, g)
            del u, v, w

            analyzer.plot_spectra(
                vspec,
                title=f"Velocity Power Spectra  |  t = {time:.4f}",
                savepath=basepath / "res_check_3d" / f"write_{write:06}.png",
                mins=vmins, maxs=vmaxs,
                field_type="velocity",
            )

            # --- Temperature spectra ---
            tspec = analyzer.compute_scalar_spectra(raw_temp, g)
            del raw_temp

            analyzer.plot_spectra(
                tspec,
                title=f"Temperature Power Spectra  |  t = {time:.4f}",
                savepath=basepath / "res_check_temp" / f"write_{write:06}.png",
                mins=tmins, maxs=tmaxs,
                field_type="temperature",
            )

            # Save k-axes once (they are the same for all t at fixed resolution)
            if not k_axes_vel:
                k_axes_vel  = {kk: vspec[kk] for kk in ("k_3d", "k_xy", "k_kx", "k_kz")}
                k_axes_temp = {kk: tspec[kk] for kk in ("k_3d", "k_xy", "k_kx", "k_kz")}

            for kk in ("E_3d", "E_xy", "E_kx", "E_kz"):
                vel[kk].append(vspec[kk])
                temp[kk].append(tspec[kk])
            times.append(time)

        # Stack time series: each list element is a 1D array (nk,)
        # → result is 2D (nt, nk)
        np.savez(
            fp.parent / f"{fp.stem}_spectra.npz",
            time=np.array(times),
            # k-axes (1D, same for all time steps)
            vel_k_3d=k_axes_vel["k_3d"],
            vel_k_xy=k_axes_vel["k_xy"],
            vel_k_kx=k_axes_vel["k_kx"],
            vel_k_kz=k_axes_vel["k_kz"],
            temp_k_3d=k_axes_temp["k_3d"],
            temp_k_xy=k_axes_temp["k_xy"],
            temp_k_kx=k_axes_temp["k_kx"],
            temp_k_kz=k_axes_temp["k_kz"],
            # E(t, k) arrays (2D: nt × nk)
            vel_E_3d=np.array(vel["E_3d"]),
            vel_E_xy=np.array(vel["E_xy"]),
            vel_E_kx=np.array(vel["E_kx"]),
            vel_E_kz=np.array(vel["E_kz"]),
            temp_E_3d=np.array(temp["E_3d"]),
            temp_E_xy=np.array(temp["E_xy"]),
            temp_E_kx=np.array(temp["E_kx"]),
            temp_E_kz=np.array(temp["E_kz"]),
        )


if __name__ == "__main__":
    from docopt import docopt
    from dedalus.tools import logging, post

    args  = docopt(__doc__)
    vmins = list(map(float, args["--vmins"].split(",")))
    vmaxs = list(map(float, args["--vmaxs"].split(",")))
    tmins = list(map(float, args["--tmins"].split(",")))
    tmaxs = list(map(float, args["--tmaxs"].split(",")))

    post.visit_writes(
        args["<files>"],
        main,
        vmins=vmins, vmaxs=vmaxs,
        tmins=tmins, tmaxs=tmaxs,
    )