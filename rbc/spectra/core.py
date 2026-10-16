"""Shared spectral machinery: Dedalus coefficient conventions, spectra, decades.

Used by ``sim/rayleigh_benard.py`` (in-simulation spectra, MPI-distributed) and by
:mod:`rbc.spectra.snapshots`, :mod:`rbc.spectra.cumulative` and
:mod:`rbc.spectra.compare`. The conventions below were checked against the Dedalus v3
source (``core/basis.py``, ``core/transforms.py``, ``tools/jacobi.py``).

**RealFourier** (x and y), size N (even), bounds (0, L):

* grid x_j = L j / N, j = 0..N-1;
* coefficients are unit-amplitude and interleaved::

      f(x) = a_0 + sum_{k=1}^{N/2-1} [ a_k cos(2 pi k x/L) + b_k (-sin(2 pi k x/L)) ]
      storage: [a_0, b_0 (= 0), a_1, b_1, a_2, b_2, ...]

  so element index e is mode number k = e // 2;
* the Nyquist mode k = N/2 is not stored (no Nyquist artefact to trim);
* mean square: <f²>_x = a_0² + (1/2) sum_k (a_k² + b_k²).

**ChebyshevT** (z), size N, bounds (0, Lz):

* grid = Gauss–Chebyshev roots of T_N, ascending: zeta_j = -cos(pi (j + 1/2) / N)
  (numpy's ``chebpts1``), z_j = Lz (zeta_j + 1) / 2; not the Gauss–Lobatto points;
* coefficients in the orthonormal Jacobi normalization: f(zeta) = sum_n c_n p_n(zeta)
  with p_0 = T_0/sqrt(pi), p_n = sqrt(2/pi) T_n, so c_0 = sqrt(pi) a_0 and
  c_n = sqrt(pi/2) a_n for the unit-amplitude a_n, and
  sum_n c_n² = ∫ f² (1 - zeta²)^(-1/2) dzeta.

**Energy per coefficient**, used everywhere::

    E[e_x, e_y, n] = w(k_x) w(k_y) c² / pi,    w(0) = 1, w(k > 0) = 1/2

In every direction this is the mean square contributed by that mode (a cosine of
amplitude a gives a²/2, a constant a_0 gives a_0²), so the three directions are on the
same footing. Summing over all coefficients gives exactly the mean of f² over the
Gauss grid (Gauss–Chebyshev quadrature with N nodes is exact for f²): the Parseval
identity of these spectra. That is a Chebyshev-weighted volume average, leaning
slightly towards the walls, which is the right measure for a resolution diagnostic,
since per-mode coefficient decay controls truncation error. Dedalus divides by N in
the forward transform, so these coefficients are already resolution independent; no
further normalization is needed.

**Spectra produced** (all are sums of E over a bin, never per-bin averages):

``E_kx[k]``, ``E_ky[k]``
    marginals: sum over the other two indices with abs(k_x) = k (k = 0..Nx/2-1).
``E_nz[n]``
    sum over k_x, k_y at Chebyshev degree n (n = 0..Nz-1).
``E_perp[s]``
    horizontal shells, s = rint(sqrt(k_x² + k_y²)).
``E_3d[s]``
    index-space shells, s = rint(sqrt(k_x² + k_y² + n²)). This mixes Fourier
    wavenumbers with Chebyshev degrees, so it is a resolution-space diagnostic, not a
    physical isotropic spectrum.
``mean_E_nz[n]``
    Chebyshev spectrum of the horizontal-mean profile.

The horizontal mean (k_x = k_y = 0) is removed from every fluctuation spectrum and
reported separately as ``mean_E_nz``. For temperature the code uses b = T - (1 - z), so
``mean_E_nz`` is the spectrum of the mean boundary-layer deviation from conduction.
With Chebyshev coefficients there is no Gibbs problem in z, so the mean thermal
boundary layer, usually the hardest feature to resolve in z, is included as a check.

Shell sums versus per-bin averages: the standard turbulence spectrum E(k) is a shell
sum (sum E(k) = total energy), which is what Kolmogorov's k^(-5/3) refers to. Dividing
each shell by its mode count removes the shell-area factor (k² in 3D, k in 2D) and
changes every slope by -2 or -1. The mode counts are saved (``count_perp``,
``count_3d``) so the per-mode average can still be recovered.
"""

import re
from pathlib import Path

import numpy as np
import scipy.fft
from scipy.integrate import trapezoid   # np.trapezoid needs numpy>=2

SQRT_PI = np.sqrt(np.pi)
FIELDS = ("vel", "temp")
SPEC_KEYS = ("E_kx", "E_ky", "E_nz", "E_perp", "E_3d", "mean_E_nz")
MODE_KEYS = ("kx", "ky", "nz", "n_perp", "n_3d", "count_perp", "count_3d")


# ===========================================================================
#  Grid -> coefficient transforms (serial; reproduce Dedalus exactly)
# ===========================================================================

def _axslice(ndim, axis, sl):
    """Index tuple selecting ``sl`` along ``axis``."""
    idx = [slice(None)] * ndim
    idx[axis] = sl
    return tuple(idx)


def realfourier_forward(data, axis):
    """Grid values to Dedalus RealFourier coefficients along one axis.

    Parameters
    ----------
    data : numpy.ndarray
        Values on the uniform grid; the axis length must be even.
    axis : int
        Axis to transform.

    Returns
    -------
    numpy.ndarray
        Interleaved (cos, -sin) unit-amplitude coefficients, Nyquist dropped.
    """
    N = data.shape[axis]
    if N % 2:
        raise ValueError("realfourier_forward assumes an even grid size.")
    F = np.fft.rfft(data, axis=axis)
    kmax = (N - 1) // 2                         # Nyquist (N/2) dropped
    F = F[_axslice(data.ndim, axis, slice(0, kmax + 1))]
    cos = F.real * (2.0 / N)
    msin = F.imag * (2.0 / N)
    cos[_axslice(data.ndim, axis, 0)] *= 0.5    # a_0 = Re F_0 / N
    msin[_axslice(data.ndim, axis, 0)] = 0.0    # invalid element, kept as zero
    shape = list(data.shape)
    shape[axis] = 2 * (kmax + 1)
    out = np.empty(shape, dtype=np.float64)
    out[_axslice(data.ndim, axis, slice(0, None, 2))] = cos
    out[_axslice(data.ndim, axis, slice(1, None, 2))] = msin
    return out


def chebyshev_forward(data, axis):
    """Grid values to orthonormal ChebyshevT coefficients along one axis.

    Same steps as Dedalus' ``FastChebyshevTransform``.

    Parameters
    ----------
    data : numpy.ndarray
        Values on Dedalus' ascending Gauss–Chebyshev grid.
    axis : int
        Axis to transform.

    Returns
    -------
    numpy.ndarray
        Coefficients c_n in the orthonormal normalization.
    """
    N = data.shape[axis]
    y = scipy.fft.dct(data, type=2, axis=axis)          # 2 sum x_j cos(pi k(2j+1)/2N)
    scale = np.full(N, np.sqrt(np.pi / 2) / N)
    scale[0] = SQRT_PI / (2 * N)
    scale[1::2] *= -1.0                                 # ascending grid -> odd sign flip
    shape = [1] * data.ndim
    shape[axis] = N
    return y * scale.reshape(shape)


def chebyshev_grid(N, bounds=(0.0, 1.0)):
    """Dedalus ChebyshevT grid (scales = 1).

    Parameters
    ----------
    N : int
        Number of modes.
    bounds : tuple of float, optional
        Interval.

    Returns
    -------
    numpy.ndarray
        Ascending Gauss–Chebyshev nodes mapped to ``bounds``.
    """
    zeta = -np.cos(np.pi * (np.arange(N) + 0.5) / N)
    a, b = bounds
    return a + (b - a) * (zeta + 1.0) / 2.0


def grid_to_coeff(field, scales=(1, 1, 1)):
    """3D grid data to Dedalus coefficients.

    Parameters
    ----------
    field : numpy.ndarray
        Values on the (nx, ny, nz) grid.
    scales : tuple of float, optional
        Grid scales the data was written at; for a dealiased grid (scales > 1) the
        padding modes are dropped.

    Returns
    -------
    numpy.ndarray
        (Nx, Ny, Nz) coefficients.
    """
    c = realfourier_forward(field, 0)
    c = realfourier_forward(c, 1)
    c = chebyshev_forward(c, 2)
    if any(s != 1 for s in scales):
        Nx, Ny, Nz = (int(round(n / s)) for n, s in zip(field.shape, scales))
        Nx, Ny = Nx - Nx % 2, Ny - Ny % 2
        c = c[:Nx, :Ny, :Nz]
    return c


# ===========================================================================
#  Coefficients -> spectra
# ===========================================================================

class SpectralReducer:
    """Reduce coefficient arrays to 1D spectra.

    Works on whatever block of coefficients the caller holds: in the simulation, the
    rank-local block with element indices from ``dist.coeff_layout.local_elements``,
    partial sums being combined with one MPI Reduce; in post-processing, the global
    array with ``ex = arange(Nx)`` and so on.

    Parameters
    ----------
    ex, ey, ez : array_like of int
        Element indices of the block along x, y, z.
    Nx, Ny, Nz : int
        Global coefficient sizes.
    comm : mpi4py.MPI.Comm, optional
        Communicator (None for serial use).
    """

    def __init__(self, ex, ey, ez, Nx, Ny, Nz, comm=None):
        self.comm = comm
        self.Nx, self.Ny, self.Nz = int(Nx), int(Ny), int(Nz)
        self.kx = np.asarray(ex, dtype=np.int64) // 2
        self.ky = np.asarray(ey, dtype=np.int64) // 2
        self.nz = np.asarray(ez, dtype=np.int64)

        self.Kx, self.Ky = self.Nx // 2, self.Ny // 2        # k = 0 .. N/2-1
        wx = np.where(self.kx == 0, 1.0, 0.5)
        wy = np.where(self.ky == 0, 1.0, 0.5)
        self.w_h = (wx[:, None] * wy[None, :]) / np.pi        # (lx, ly)
        self.mean_mask = (self.kx[:, None] == 0) & (self.ky[None, :] == 0)

        self.kh2 = (self.kx[:, None] ** 2 + self.ky[None, :] ** 2).astype(np.float64)
        self.perp_idx = np.rint(np.sqrt(self.kh2)).astype(np.int64).ravel()

        kxm, kym, nzm = self.Kx - 1, self.Ky - 1, self.Nz - 1
        self.n_perp_len = int(np.rint(np.hypot(kxm, kym))) + 1
        self.n_3d_len = int(np.rint(np.sqrt(kxm**2 + kym**2 + nzm**2))) + 1
        # Largest shell radius that is fully inside the coefficient box;
        # shells beyond it are only partly filled and drop off artificially.
        self.perp_complete = int(min(kxm, kym))
        self.n3d_complete = int(min(kxm, kym, nzm))
        self._modes = self._global_modes()
        self._buf_len = (self.Kx + self.Ky + 2 * self.Nz
                         + self.n_perp_len + self.n_3d_len)

    # -- static mode tables (computed from global sizes, no MPI needed) -----
    def _global_modes(self):
        """Mode-number axes and shell mode counts (global, no MPI needed)."""
        kx = np.arange(self.Kx)
        ky = np.arange(self.Ky)
        dof_x = np.where(kx == 0, 1, 2)              # cos only for k=0
        dof_y = np.where(ky == 0, 1, 2)
        s2 = np.rint(np.hypot(kx[:, None], ky[None, :])).astype(np.int64)
        count_perp = np.bincount(s2.ravel(), weights=(dof_x[:, None] * dof_y[None, :]).ravel(),
                                 minlength=self.n_perp_len)
        count_3d = np.zeros(self.n_3d_len)
        kh2 = (kx[:, None] ** 2 + ky[None, :] ** 2).astype(float)
        dof_h = (dof_x[:, None] * dof_y[None, :]).ravel()
        for n in range(self.Nz):
            s3 = np.rint(np.sqrt(kh2 + n * n)).astype(np.int64).ravel()
            count_3d += np.bincount(s3, weights=dof_h, minlength=self.n_3d_len)
        return dict(
            kx=kx.astype(float), ky=ky.astype(float), nz=np.arange(self.Nz, dtype=float),
            n_perp=np.arange(self.n_perp_len, dtype=float),
            n_3d=np.arange(self.n_3d_len, dtype=float),
            count_perp=count_perp, count_3d=count_3d,
        )

    @property
    def modes(self):
        """dict: mode-number axes (``kx``, ``ky``, ``nz``, ``n_perp``, ``n_3d``) and shell
        mode counts (``count_perp``, ``count_3d``)."""
        return self._modes

    # -- core -------------------------------------------------------------
    def power(self, components):
        """Energy per coefficient, summed over vector components.

        Parameters
        ----------
        components : list of numpy.ndarray
            Coefficient blocks (one per component).

        Returns
        -------
        numpy.ndarray
            E = sum_i w c_i² / pi on this block.
        """
        P = None
        for c in components:
            c = np.asarray(c)
            P = c * c if P is None else P + c * c
        P *= self.w_h[:, :, None]
        return P

    def reduce(self, P):
        """Bin a block of energies into the global 1D spectra.

        Parameters
        ----------
        P : numpy.ndarray
            Energies of the local block (modified in place: the mean column is zeroed).

        Returns
        -------
        dict or None
            ``E_kx``, ``E_ky``, ``E_nz``, ``mean_E_nz``, ``E_perp``, ``E_3d`` on rank 0
            (summed over ranks with one Reduce); None on other ranks.
        """
        # 1) Horizontal-mean column: Chebyshev spectrum of the mean profile.
        mean_part = P[self.mean_mask]                       # (n_mean, lz)
        mean_E_nz = np.bincount(self.nz, weights=mean_part.sum(axis=0),
                                minlength=self.Nz) if mean_part.size else np.zeros(self.Nz)
        P[self.mean_mask] = 0.0                             # fluctuations only from here

        # 2) Marginals: collapse the other two axes first, then bin (cheap).
        E_kx = np.bincount(self.kx, weights=P.sum(axis=(1, 2)), minlength=self.Kx)
        E_ky = np.bincount(self.ky, weights=P.sum(axis=(0, 2)), minlength=self.Ky)
        E_nz = np.bincount(self.nz, weights=P.sum(axis=(0, 1)), minlength=self.Nz)

        # 3) Horizontal shells.
        E_perp = np.bincount(self.perp_idx, weights=P.sum(axis=2).ravel(),
                             minlength=self.n_perp_len)

        # 4) 3D index-space shells, one z-slab at a time (no 3D int array).
        E_3d = np.zeros(self.n_3d_len)
        for j, n in enumerate(self.nz):
            s3 = np.rint(np.sqrt(self.kh2 + float(n * n))).astype(np.int64).ravel()
            E_3d += np.bincount(s3, weights=P[:, :, j].ravel(), minlength=self.n_3d_len)

        local = np.concatenate([E_kx, E_ky, E_nz, mean_E_nz, E_perp, E_3d])
        if self.comm is not None and self.comm.size > 1:
            total = np.zeros_like(local) if self.comm.rank == 0 else None
            self.comm.Reduce(local, total, root=0)
            if self.comm.rank != 0:
                return None
        else:
            total = local
        return self._split(total)

    def _split(self, buf):
        """Cut the reduced buffer back into named spectra."""
        out, i = {}, 0
        for key, n in (("E_kx", self.Kx), ("E_ky", self.Ky), ("E_nz", self.Nz),
                       ("mean_E_nz", self.Nz), ("E_perp", self.n_perp_len),
                       ("E_3d", self.n_3d_len)):
            out[key] = buf[i:i + n].copy()
            i += n
        return out


# ===========================================================================
#  HDF5 output (one file per run segment / snapshot set; rank 0 only)
# ===========================================================================

def create_spectra_file(path, reducer, attrs):
    """Create an empty spectra file with resizable time-series datasets.

    Parameters
    ----------
    path : pathlib.Path
        File to create.
    reducer : SpectralReducer
        Supplies the mode tables and spectrum lengths.
    attrs : dict
        File attributes (Ra, Pr, box size, Nx, Ny, Nz, source).
    """
    import h5py
    with h5py.File(path, "w") as f:
        for k, v in attrs.items():
            f.attrs[k] = v
        f.attrs["perp_complete"] = reducer.perp_complete
        f.attrs["n3d_complete"] = reducer.n3d_complete
        f.attrs["energy_definition"] = (
            "E = w(kx) w(ky) c^2 / pi per Dedalus coefficient; w(0)=1, w(k>0)=1/2; "
            "fluctuation spectra exclude kx=ky=0; spectra are bin SUMS")
        g = f.create_group("modes")
        for k, v in reducer.modes.items():
            g.create_dataset(k, data=v)
        s = f.create_group("scales")
        s.create_dataset("sim_time", shape=(0,), maxshape=(None,), dtype="f8")
        s.create_dataset("iteration", shape=(0,), maxshape=(None,), dtype="i8")
        s.create_dataset("grid_mean_square_vel", shape=(0,), maxshape=(None,), dtype="f8")
        s.create_dataset("grid_mean_square_temp", shape=(0,), maxshape=(None,), dtype="f8")
        lengths = dict(E_kx=reducer.Kx, E_ky=reducer.Ky, E_nz=reducer.Nz,
                       mean_E_nz=reducer.Nz, E_perp=reducer.n_perp_len,
                       E_3d=reducer.n_3d_len)
        for fld in FIELDS:
            g = f.create_group(fld)
            for key, n in lengths.items():
                g.create_dataset(key, shape=(0, n), maxshape=(None, n), dtype="f8",
                                 chunks=(64, n))


def append_spectra(path, sim_time, iteration, vel, temp, gms_vel=np.nan, gms_temp=np.nan):
    """Append one write to a spectra file.

    The file is opened and closed every time, so a killed job never leaves a corrupt
    or unflushed file.

    Parameters
    ----------
    path : pathlib.Path
        Spectra file.
    sim_time : float
        Simulation time.
    iteration : int
        Solver iteration.
    vel, temp : dict
        Spectra from :meth:`SpectralReducer.reduce`.
    gms_vel, gms_temp : float, optional
        Grid mean squares, for the Parseval check (snapshot-derived files only).
    """
    import h5py
    with h5py.File(path, "a") as f:
        nt = f["scales/sim_time"].shape[0]
        for name, val in (("sim_time", sim_time), ("iteration", iteration),
                          ("grid_mean_square_vel", gms_vel),
                          ("grid_mean_square_temp", gms_temp)):
            d = f["scales"][name]
            d.resize((nt + 1,))
            d[nt] = val
        for fld, spec in (("vel", vel), ("temp", temp)):
            for key in SPEC_KEYS:
                d = f[fld][key]
                d.resize((nt + 1, d.shape[1]))
                d[nt] = spec[key]


def next_set_path(folder, prefix):
    """Path of the next set file.

    Parameters
    ----------
    folder : str or pathlib.Path
        Output folder.
    prefix : str
        Set name prefix.

    Returns
    -------
    pathlib.Path
        ``folder/prefix_s{N+1}.h5`` with N the largest existing set number.
    """
    folder = Path(folder)
    nums = [int(m.group(1)) for p in folder.glob(f"{prefix}_s*.h5")
            if (m := re.search(rf"{prefix}_s(\d+)\.h5$", p.name))]
    return folder / f"{prefix}_s{max(nums, default=0) + 1}.h5"


# ===========================================================================
#  In-simulation driver (used by rbc_script.py)
# ===========================================================================

class InSimSpectra:
    """Spectra computed during the simulation from the solver's coefficients.

    Every ``sim_dt`` of simulation time, appends the spectra of u and b to
    ``outdir/spectra_s<N>.h5`` (a new set per job). Cost per write: one pass over the
    local coefficients, a few ``np.bincount`` calls and one MPI Reduce of a vector of
    length ~5N, far cheaper than a timestep; output is a few kB per write instead of
    a full 3D field.

    Parameters
    ----------
    dist : dedalus.core.distributor.Distributor
        The simulation's distributor.
    u, b : dedalus.core.field.Field
        Velocity and temperature perturbation fields.
    sim_dt : float
        Simulation time between writes.
    outdir : str or pathlib.Path
        Output folder.
    attrs : dict
        Attributes stored in the file (Ra, Pr, box size).
    comm : mpi4py.MPI.Comm
        Communicator.
    """

    def __init__(self, dist, u, b, sim_dt, outdir, attrs, comm):
        self.u, self.b, self.sim_dt, self.comm = u, b, float(sim_dt), comm
        layout = dist.coeff_layout
        ex, ey, ez = layout.local_elements(b.domain, scales=1)
        gshape = layout.global_shape(b.domain, scales=1)
        self.red = SpectralReducer(ex, ey, ez, *gshape, comm=comm)
        self.path = None
        if comm.rank == 0:
            outdir = Path(outdir)
            outdir.mkdir(parents=True, exist_ok=True)
            self.path = next_set_path(outdir, "spectra")
            create_spectra_file(self.path, self.red, dict(attrs, source="simulation",
                                                          Nx=gshape[0], Ny=gshape[1], Nz=gshape[2]))
        self.last_div = None

    def maybe_write(self, solver):
        """Write spectra if a new ``sim_dt`` interval has started.

        Call once per iteration on every rank. Every rank evaluates the same condition
        (``sim_time`` is identical across ranks), so the collective Reduce is safe.

        Parameters
        ----------
        solver : dedalus.core.solvers.InitialValueSolver
            The running solver.
        """
        div = int(np.floor(solver.sim_time / self.sim_dt))
        if self.last_div is None:          # first call: do not re-write on restart
            self.last_div = div
            return
        if div <= self.last_div:
            return
        self.last_div = div
        self.write(solver)

    def write(self, solver):
        """Compute and append one write (collective).

        Parameters
        ----------
        solver : dedalus.core.solvers.InitialValueSolver
            The running solver.
        """
        uc = self.u["c"]                                   # (3, lx, ly, lz) coefficient view
        vel = self.red.reduce(self.red.power([uc[0], uc[1], uc[2]]))
        temp = self.red.reduce(self.red.power([self.b["c"]]))
        if self.comm.rank == 0:
            append_spectra(self.path, solver.sim_time, solver.iteration, vel, temp)


# ===========================================================================
#  Loading spectra files (post-processing)
# ===========================================================================

def _set_number(p):
    """Set number from a file name (0 if none)."""
    m = re.search(r"_s(\d+)(?:_spectra)?\.h5$", p.name)
    return int(m.group(1)) if m else 0


def spectra_files(basepath, source="auto"):
    """Spectra files of a run.

    Parameters
    ----------
    basepath : str or pathlib.Path
        Run directory.
    source : {"auto", "sim", "snapshots"}, optional
        In-simulation files (``spectra/spectra_s*.h5``), snapshot-derived files
        (``snapshots/snapshots_s*_spectra.h5``), or in-simulation if any exist.

    Returns
    -------
    list of pathlib.Path
        Files in set order.
    """
    basepath = Path(basepath)
    sim = sorted((basepath / "spectra").glob("spectra_s*.h5"), key=_set_number)
    snap = sorted((basepath / "snapshots").glob("snapshots_s*_spectra.h5"), key=_set_number)
    if source == "sim":
        return sim
    if source == "snapshots":
        return snap
    return sim if sim else snap


def load_spectra(files):
    """Concatenate spectra files.

    Only strictly increasing times are kept, which drops frames repeated by a restart
    from an earlier checkpoint.

    Parameters
    ----------
    files : list of pathlib.Path
        Spectra files in set order.

    Returns
    -------
    time : numpy.ndarray
        Shape (nt,).
    modes : dict
        Mode tables of the first file.
    data : dict
        ``data[field][key]`` with shape (nt, n), field in :data:`FIELDS`.
    attrs : dict
        Attributes of the first file.

    Raises
    ------
    FileNotFoundError
        If ``files`` is empty.
    """
    import h5py
    if not files:
        raise FileNotFoundError("No spectra files found.")
    time = np.empty(0)
    data = {fld: {k: [] for k in SPEC_KEYS} for fld in FIELDS}
    modes, attrs = None, None
    for fp in files:
        with h5py.File(fp, "r") as f:
            if modes is None:
                modes = {k: f["modes"][k][:] for k in f["modes"]}
                attrs = dict(f.attrs)
            t = f["scales/sim_time"][:]
            keep = t > (time[-1] if time.size else -np.inf)
            # also guard against non-monotone times inside one file
            keep &= np.concatenate([[True], np.diff(t) > 0]) if t.size else keep
            if not keep.any():
                continue
            time = np.concatenate([time, t[keep]])
            for fld in FIELDS:
                for k in SPEC_KEYS:
                    data[fld][k].append(f[fld][k][:][keep])
    for fld in FIELDS:
        for k in SPEC_KEYS:
            data[fld][k] = np.concatenate(data[fld][k], axis=0)
    return time, modes, data, attrs


# ===========================================================================
#  Panels, time averaging and the "decades" metrics
# ===========================================================================

def panel_series(modes, spec, attrs):
    """The four plotted spectra of one field.

    Horizontal is the mean of the x and y marginals (identical statistics in a square
    box; averaging halves the noise). Shells are cut at the last complete shell.

    Parameters
    ----------
    modes : dict
        Mode tables.
    spec : dict
        Spectra of one field, each shaped (..., n).
    attrs : dict
        File attributes (``perp_complete``, ``n3d_complete``).

    Returns
    -------
    dict
        ``"3d"``, ``"nz"``, ``"h"``, ``"perp"`` -> ``(n, E)``.
    """
    pc = int(attrs.get("perp_complete", len(modes["n_perp"]) - 1))
    c3 = int(attrs.get("n3d_complete", len(modes["n_3d"]) - 1))
    nh = min(len(modes["kx"]), len(modes["ky"]))
    E_h = 0.5 * (spec["E_kx"][..., :nh] + spec["E_ky"][..., :nh])
    return {
        "3d":   (modes["n_3d"][:c3 + 1], spec["E_3d"][..., :c3 + 1]),
        "nz":   (modes["nz"], spec["E_nz"]),
        "h":    (modes["kx"][:nh], E_h),
        "perp": (modes["n_perp"][:pc + 1], spec["E_perp"][..., :pc + 1]),
    }


def time_average(time, E):
    """Trapezoidal time average.

    Trapezoid rather than Simpson: the samples are irregular (writes land on the first
    step after each interval) and the signal is not smooth on the sampling scale, so
    Simpson buys no accuracy and its alternating 4/3, 2/3 weights add noise.

    Parameters
    ----------
    time : numpy.ndarray
        Shape (nt,).
    E : numpy.ndarray
        Shape (nt, n).

    Returns
    -------
    numpy.ndarray
        Shape (n,).
    """
    if len(time) < 2 or time[-1] <= time[0]:
        return E.mean(axis=0)
    return trapezoid(E, time, axis=0) / (time[-1] - time[0])


def decades(n, E, n_min=1, tail_pts=3):
    """Dynamic range of a spectrum, in decades.

    Works on the last axis, so ``E`` can be one spectrum (n,) or a time series (nt, n).

    * ``D`` = peak - tail (recommended): peak = max of log10 E over n >= ``n_min``,
      tail = mean of log10 E over the last ``tail_pts`` modes.
    * ``D_maxmin`` = peak - min(log10 E) (the old definition).
    * ``pileup`` = tail - (smallest ``tail_pts``-mode running mean of log10 E): 0 if the
      spectrum is lowest at its end; above ~0.1 energy is piling up at the grid scale.

    Why not max - min: with a pile-up the minimum sits before the cutoff, so max - min
    counts the dip and reports more decades for a worse-resolved run. Why not the last
    point alone: one mode is noisy. The mean of log E over the last few modes is robust
    to that and still penalizes an uptick at the end.

    Parameters
    ----------
    n : array_like
        Mode numbers.
    E : array_like
        Energies, modes along the last axis.
    n_min : int, optional
        Smallest mode considered for the peak (1 skips the mean).
    tail_pts : int, optional
        Modes averaged for the tail.

    Returns
    -------
    dict
        ``D``, ``D_maxmin`` and ``pileup`` (arrays if ``E`` is 2D).
    """
    n = np.asarray(n)
    E = np.asarray(E, dtype=float)
    sel = n >= n_min
    Es = np.where(E[..., sel] > 0, E[..., sel], np.nan)
    with np.errstate(divide="ignore", invalid="ignore"):
        L = np.log10(Es)
        peak = np.nanmax(L, axis=-1)
        tail = np.nanmean(L[..., -tail_pts:], axis=-1)
        k = np.ones(tail_pts) / tail_pts
        run = np.apply_along_axis(lambda v: np.convolve(v, k, mode="valid"), -1, L)
        lo_smooth = np.nanmin(run, axis=-1)
        lo = np.nanmin(L, axis=-1)
    return dict(D=peak - tail, D_maxmin=peak - lo, pileup=tail - lo_smooth)


# ===========================================================================
#  Plot helpers
# ===========================================================================

def setup_yaxis(ax, ymin, ymax):
    """Log y axis with two or three labelled ticks and a dashed line at every decade.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes.
    ymin, ymax : float
        Limits in log10 units.
    """
    import matplotlib.ticker as ticker
    ax.set_ylim(10**ymin, 10**ymax)
    ntick = 3 if (ymax - ymin) > 2 else 2
    log_ticks = np.linspace(ymin, ymax, ntick)
    ax.yaxis.set_major_locator(ticker.FixedLocator(10.0**log_ticks))
    labels = [rf"$10^{{{f'{l:.1f}'.rstrip('0').rstrip('.')}}}$" for l in log_ticks]
    ax.yaxis.set_major_formatter(ticker.FixedFormatter(labels))
    dec = np.arange(int(np.floor(ymin)), int(np.ceil(ymax)) + 1)
    ax.yaxis.set_minor_locator(ticker.FixedLocator(10.0**dec))
    ax.yaxis.set_minor_formatter(ticker.NullFormatter())
    for e in dec:
        ax.axhline(10.0**e, linestyle="--", color="gray", alpha=0.5, linewidth=0.8)


def auto_ylim(*arrays, pad=0.5):
    """Log10 y limits covering every positive finite value.

    Parameters
    ----------
    *arrays : array_like
        Data to cover.
    pad : float, optional
        Padding in decades.

    Returns
    -------
    tuple of float
        ``(ymin, ymax)`` in log10 units.
    """
    vals = np.concatenate([np.ravel(a) for a in arrays])
    vals = vals[np.isfinite(vals) & (vals > 0)]
    if vals.size == 0:
        return -1.0, 1.0
    return np.log10(vals.min()) - pad, np.log10(vals.max()) + pad


PANEL_META = {
    #        position, title, xlabel
    "3d":   ((0, 0), "3D index-space shells",
             r"$s = \mathrm{rint}\sqrt{k_x^2+k_y^2+n_z^2}$"),
    "nz":   ((0, 1), "Vertical (Chebyshev degree)", r"$n_z$  [max $= N_z-1$]"),
    "h":    ((1, 0), r"Horizontal marginal ($\frac{1}{2}(E_x+E_y)$)",
             r"$k_{x,y}$  [mode number, max $= N_x/2-1$]"),
    "perp": ((1, 1), "Horizontal shells",
             r"$k_\perp = \mathrm{rint}\sqrt{k_x^2+k_y^2}$"),
}

# Plain-text panel names for text reports.
REPORT_NAMES = {"3d": "3D index-space shells", "nz": "vertical (Chebyshev degree)",
                "h": "horizontal marginal (x,y avg)", "perp": "horizontal shells"}