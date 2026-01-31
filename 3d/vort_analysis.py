"""
Checks the magnitude of the vorticity of snapshot files for negative values and volume averages.

Usage:
    vort_analysis.py <files>... [--mins=<mins>] [--maxs=<maxs>]

Options:
    --mins=<mins>   Minimum log10 limit. String of csv floats.
    --maxs=<maxs>   Maximum log10 limit. String of csv floats.
"""



import h5py as h5
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt




def main(file, start, count):
    fp = Path(file)
    basepath = Path(fp.parents[1])
    with h5.File(fp, 'r') as f:
        # Collect the file data (start/count is for parallelization)

        time = np.array(f["scales/sim_time"][start:start+count])
        writes = np.array(f["scales/write_number"][start:start+count])
        omega = np.array(f['tasks']['vorticity'][start:start+count])
        nt, nx, ny, nz = omega.shape       
    
    # Discard the first timestep
    c=0
    if np.isclose(time[0],0):
        if nt == 1:
            return
        time = time[1:]
        omega = omega[1:]
        nt -= 1
        c = 1

    
    # Main loop: interpolate, FFT, bin, plot
    for ti in range(nt):
        # Interpolators for current time slice
        interp_omega = RegularGridInterpolator((x, y, z_cheb), omega[ti], method='linear', bounds_error=False, fill_value=0)

        # Interpolate onto uniform z grid
        omega_uniform = interp_u(pts).reshape((nx, ny, nz))


        # Compute 3D FFT and energy. Supposedly the ortho norm negates the
        # need for renormalization. No dividing by resolution :)
        Omega = np.fft.fftn(omega_uniform, norm='ortho')




        # Bin power spectrum into magnitude of 3d wavevector - Full 3d spectra
        E3d = np.abs(Omega)**2
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

        pars_grid = np.sum(omega_uniform**2)
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
        savename=f"{str(basepath)}/res_check_vort/write_{writes[ti]+c:06}.png"
        fig.savefig(savename)
        fig.clear()



if __name__ == "__main__":

    from docopt import docopt
    from dedalus.tools import logging
    from dedalus.tools import post

    args = docopt(__doc__)



    post.visit_writes(
        args["<files>"],
        main
    )

