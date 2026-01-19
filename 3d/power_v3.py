"""
Script to perform analysis tasks...

Usage:
    power_v3.py <files>... [--mins=<mins>] [--maxs=<maxs>]

Options:
    --mins=<mins>   Minimum log10 limit. String of csv floats.
    --maxs=<maxs>   Maximum log10 limit. String of csv floats.
"""



import h5py as h5
import numpy as np
from numpy.polynomial.chebyshev import chebpts2
from pathlib import Path
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator

#TODO: get this working for temp as well!
# Questions:
# Should I be checking that the energy is the same as the energy computed using the analysis files?
# Any extra multiplicative factors that I forgot? does it really matter if things are scaled?

def main(file, start, count, mins, maxs):
    fp = Path(file)
    basepath = Path(fp.parents[1])
    with h5.File(fp, 'r') as f:
        # Collect the state file data (start/count is for parallelization)
        time = np.array(f["scales/sim_time"][start:start+count])
        writes = np.array(f["scales/write_number"][start:start+count])
        u = np.array(f['tasks']['u'][start:start+count, 0])
        v = np.array(f['tasks']['u'][start:start+count, 1])
        w = np.array(f['tasks']['u'][start:start+count, 2])
        nt, nx, ny, nz = u.shape       
        
    # Discard the first timestep
    c=0
    if np.isclose(time[0],0):
        time = time[1:]
        u = u[1:]
        v = v[1:]
        w = w[1:]
        nt -= 1
        c = 1

    # Physical grid for interpolation
    x = np.linspace(0, 2, nx, endpoint=False)
    y = np.linspace(0, 2, ny, endpoint=False)
    z_cheb = chebpts2(nz) * 0.5
    z_uniform = np.linspace(-1/2, 1/2, nz, endpoint=False)

    # Points for interpolation - see this google group post:
    # https://groups.google.com/d/msgid/dedalus-users/7040a52e-83fc-49fa-90e5-8a8820ddd7ea%40googlegroups.com
    X, Y, Z = np.meshgrid(x, y, z_uniform, indexing='ij')
    pts = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=-1)

    # Prepare k-space grids
    kx = np.fft.fftfreq(nx, d=2/nx)
    ky = np.fft.fftfreq(ny, d=2/ny)
    kz = np.fft.fftfreq(nz, d=1/nz)

    # Get all possible combinations of wavenumbers
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing='ij')
    Knorm_full = np.sqrt(KX**2 + KY**2 + KZ**2).flatten()
    Kmax = int(np.ceil(Knorm_full.max()))
    bins = np.arange(0.5, Kmax + 1, 1.0)
    k_centers = 0.5 * (bins[:-1] + bins[1:])
    eps = 1e-12
    
    # Main loop: interpolate, FFT, bin, plot
    for ti in range(nt):
        # Interpolators for current time slice
        interp_u = RegularGridInterpolator((x, y, z_cheb), u[ti], method='linear', bounds_error=False, fill_value=0)
        interp_v = RegularGridInterpolator((x, y, z_cheb), v[ti], method='linear', bounds_error=False, fill_value=0)
        interp_w = RegularGridInterpolator((x, y, z_cheb), w[ti], method='linear', bounds_error=False, fill_value=0)

        # Interpolate onto uniform z grid
        u_uniform = interp_u(pts).reshape((nx, ny, nz))
        v_uniform = interp_v(pts).reshape((nx, ny, nz))
        w_uniform = interp_w(pts).reshape((nx, ny, nz))

        # Compute 3D FFT and energy. Supposedly the ortho norm negates the
        # need for renormalization. No dividing by resolution :)
        U = np.fft.fftn(u_uniform, norm='ortho')
        V = np.fft.fftn(v_uniform, norm='ortho')
        W = np.fft.fftn(w_uniform, norm='ortho')



        # Bin power spectrum into magnitude of 3d wavevector - Full 3d spectra
        E3d = np.abs(U)**2 + np.abs(V)**2 + np.abs(W)**2
        E_k_flat = E3d.flatten()
        pk, _ = np.histogram(Knorm_full, bins=bins, weights=E_k_flat)
        nm, _ = np.histogram(Knorm_full, bins=bins)
        E_k = pk / (nm + eps)


        # Bin power spectrum into z wavenumber - vertical marginal
        E_kz = np.sum(E3d, axis=(0,1))
        kz_plot = np.abs(kz)

        # Bin power into x wavenumber - horizontal marginal
        E_kx = np.sum(E3d, axis=(1,2))
        kx_plot = np.abs(kx)

        # Bin power into xy wavevector - horizontal marginal
        E_xy = np.sum(E3d, axis=2)
        KX2D, KY2D = np.meshgrid(kx, ky, indexing='ij')
        Kperp = np.sqrt(KX2D**2 + KY2D**2)

        E_xy_flat = E_xy.ravel()
        Kxy_flat = Kperp.ravel()
        bins_xy = np.arange(0.5, Kperp.max() + 1.0, 1.0)

        pk_xy, _ = np.histogram(Kxy_flat, bins=bins_xy, weights=E_xy_flat)
        nm_xy, _ = np.histogram(Kxy_flat, bins=bins_xy)

        E_kxy = pk_xy / (nm_xy + eps)
        k_xy_centers = 0.5 * (bins_xy[:-1] + bins_xy[1:])

        pars_grid = np.sum(u_uniform**2 + v_uniform**2 + w_uniform**2)
        pars_3d = np.sum(E3d)
        pars_kz = np.sum(E_kz)
        pars_kx = np.sum(E_kx)
        pars_kxy = np.sum(E_xy_flat)                

        fig, axes = plt.subplots(2, 2, figsize=(22, 18))

        # --- Dealiasing cutoffs (2/3 rule) ---
        kx_cut = (2/3) * np.max(np.abs(kx))
        ky_cut = (2/3) * np.max(np.abs(ky))
        kz_cut = (2/3) * np.max(np.abs(kz))

        # For isotropic 3D and planar spectra, use the smallest relevant cutoff
        k3d_cut = min(kx_cut, ky_cut, kz_cut)
        kxy_cut = min(kx_cut, ky_cut)

        # --- Top-left: full 3D isotropic spectrum ---
        ax = axes[0, 0]
        ax.loglog(k_centers, E_k, '.-')
        ax.axvline(k3d_cut, color='k', linestyle='--', alpha=0.7, label='2/3 cutoff')
        ax.set_xlabel(r'$k$')
        ax.set_ylabel(r'$E(k)$')
        ax.set_title('3D isotropic spectrum')
        ax.set_ylim((10**float(mins[0]), 10**float(maxs[0])))
        ax.legend(frameon=False)

        # --- Top-right: vertical spectrum (kz) ---
        ax = axes[0, 1]
        ax.loglog(np.abs(kz), E_kz, '.-')
        ax.axvline(kz_cut, color='k', linestyle='--', alpha=0.7)
        ax.set_xlabel(r'$k_z$')
        ax.set_ylabel(r'$E(k_z)$')
        ax.set_title('Vertical spectrum')
        ax.set_ylim((10**float(mins[1]), 10**float(maxs[1])))

        # --- Bottom-left: horizontal spectrum (kx) ---
        ax = axes[1, 0]
        ax.loglog(np.abs(kx), E_kx, '.-')
        ax.axvline(kx_cut, color='k', linestyle='--', alpha=0.7)
        ax.set_xlabel(r'$k_x$')
        ax.set_ylabel(r'$E(k_x)$')
        ax.set_title('Horizontal spectrum (x)')
        ax.set_ylim((10**float(mins[2]), 10**float(maxs[2])))

        # --- Bottom-right: 2D planar spectrum (xy) ---
        ax = axes[1, 1]
        ax.loglog(k_xy_centers, E_kxy, '.-')
        ax.axvline(kxy_cut, color='k', linestyle='--', alpha=0.7)
        ax.set_xlabel(r'$k = \sqrt{k_x^2 + k_y^2}$')
        ax.set_ylabel(r'$E_{xy}(k_{xy})$')
        ax.set_title('Horizontal planar spectrum (xy)')
        ax.set_ylim((10**float(mins[3]), 10**float(maxs[3])))


        # --- Global title with energy check ---
        fig.suptitle(
            f't: {time[ti]:.4f}, grid: {pars_grid:.5e}, 3d: {pars_3d:.5e}, z: {pars_kz:.5e}, x: {pars_kx:.5e}, xy: {pars_kxy:.5e}',
            fontsize=18
        )

        fig.tight_layout(rect=[0, 0, 1, 0.95])
        savename=f"{str(basepath)}/res_check_3d/write_{writes[ti]+c:06}.png"
        fig.savefig(savename)
        fig.clear()



if __name__ == "__main__":

    from docopt import docopt
    from dedalus.tools import logging
    from dedalus.tools import post

    args = docopt(__doc__)
    mins = list(map(float, args["--mins"].split(",")))
    maxs = list(map(float, args["--maxs"].split(",")))



    post.visit_writes(
        args["<files>"],
        main,
        mins=mins,
        maxs=maxs
    )

