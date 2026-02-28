"""
Script to perform analysis tasks...

Usage:
    power_temp.py <files>... [--mins=<mins>] [--maxs=<maxs>]

Options:
    --mins=<mins>   Minimum log10 limit. String of csv floats.
    --maxs=<maxs>   Maximum log10 limit. String of csv floats.
"""



import h5py as h5
import numpy as np
from numpy.polynomial.chebyshev import chebpts2
from pathlib import Path
import matplotlib
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator
import matplotlib.ticker as ticker
matplotlib.use("Agg")
s = 30

plt.rcParams.update({"font.size": 0.75*s})
plt.ioff()


# Questions:
# Should I be checking that the energy is the same as the energy computed using the analysis files?

def main(file, start, count, mins, maxs):
    fp = Path(file)
    basepath = Path(fp.parents[1])
        
    prev_nx, prev_ny, prev_nz = (0,0,0)


    with h5.File(fp, 'r') as f:
        for ti in range(start, start+count):
            # Collect the file data (start/count is for parallelization)
            time = f["scales/sim_time"][ti]
            write = f["scales/write_number"][ti]
            temp = np.array(f['tasks']['temperature'][ti])

            nx, ny, nz = temp.shape       
            
            # Discard the first timestep
            c=0
            if np.isclose(time,0):
                continue

            # Only build the x,y,z and kx,ky,kz arrays if the resolution has changed; will never happen within a file, but it's easiest this way rather than finding some way to take this out of the loop
            if prev_nx != nx or prev_ny != ny or prev_nz != nz:
                prev_nx, prev_ny, prev_nz = nx, ny, nz
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
                eps = 2 * np.finfo(np.float64).eps


            # Main loop: interpolate, FFT, bin, plot
            # Interpolators for current time slice
            interp_temp = RegularGridInterpolator((x, y, z_cheb), temp, method='linear', bounds_error=False, fill_value=0)

            # Interpolate onto uniform z grid
            temp_uniform = interp_temp(pts).reshape((nx, ny, nz))
            pars_grid = np.sum(temp_uniform**2)

            # Compute 3D FFT and energy. Supposedly the ortho norm negates the
            # need for renormalization. No dividing by resolution :)
            Temp = np.fft.fftn(temp_uniform, norm='ortho')
            E_3d = np.abs(Temp)**2

            del interp_temp, temp_uniform, Temp



            # Bin power spectrum into magnitude of 3d wavevector - Full 3d spectra
            E_k_flat = E_3d.flatten()
            pk, _ = np.histogram(Knorm_full, bins=bins, weights=E_k_flat)
            nm, _ = np.histogram(Knorm_full, bins=bins)
            E_k = pk / (nm + eps)


            # Bin power spectrum into z wavenumber - vertical marginal
            E_kz = np.sum(E_3d, axis=(0,1))
            kz_plot = np.abs(kz)

            # Bin power into x wavenumber - horizontal marginal
            E_kx = np.sum(E_3d, axis=(1,2))
            kx_plot = np.abs(kx)

            # Bin power into xy wavevector - horizontal marginal
            E_xy = np.sum(E_3d, axis=2)
            KX2D, KY2D = np.meshgrid(kx, ky, indexing='ij')
            Kperp = np.sqrt(KX2D**2 + KY2D**2)

            E_xy_flat = E_xy.ravel()
            Kxy_flat = Kperp.ravel()
            bins_xy = np.arange(0.5, Kperp.max() + 1.0, 1.0)

            pk_xy, _ = np.histogram(Kxy_flat, bins=bins_xy, weights=E_xy_flat)
            nm_xy, _ = np.histogram(Kxy_flat, bins=bins_xy)

            E_kxy = pk_xy / (nm_xy + eps)
            k_xy_centers = 0.5 * (bins_xy[:-1] + bins_xy[1:])

            pars_3d = np.sum(E_3d)
            #pars_kz = np.sum(E_kz)
            #pars_kx = np.sum(E_kx)
            #pars_kxy = np.sum(E_xy_flat)                


            # --- Dealiasing cutoffs (2/3 rule) --- IDK IF THIS IS USEFUL
            #kx_cut = (2/3) * np.max(np.abs(kx))
            #ky_cut = (2/3) * np.max(np.abs(ky))
            #kz_cut = (2/3) * np.max(np.abs(kz))

            # For isotropic 3D and planar spectra, use the smallest relevant cutoff
            #k3d_cut = min(kx_cut, ky_cut, kz_cut)
            #kxy_cut = min(kx_cut, ky_cut)

            fig, axes = plt.subplots(2, 2, figsize=(16, 14), layout='constrained')

            # --- Top-left: full 3D isotropic spectrum ---
            ax = axes[0, 0]
            ax.loglog(k_centers, E_k, '.-')
            #ax.axvline(k3d_cut, color='k', linestyle='--', alpha=0.7, label='2/3 cutoff')
            ax.set_xlabel(r'$k$')
            ax.set_ylabel(r'$E(k)$')
            ax.set_title(r'3D isotropic spectrum')
            ax.set_ylim((10**float(mins[0]), 10**float(maxs[0])))
            ntick=4
            if maxs[0]-mins[0] < 4:
                ntick = int(maxs[0]-mins[0]) + 1
            log_ticks = np.linspace(np.log10(10**float(mins[0])), np.log10(10**float(maxs[0])), ntick)
            yticks = 10**log_ticks
            ax.yaxis.set_major_locator(ticker.FixedLocator(yticks))
            labels = [rf"$10^{{{int(np.floor(l))}}}$"
                    for l in log_ticks]

            ax.yaxis.set_major_formatter(ticker.FixedFormatter(labels))
            ax.yaxis.set_minor_locator(ticker.NullLocator())


            # --- Top-right: vertical spectrum (kz) ---
            ax = axes[0, 1]
            ax.loglog(np.abs(kz), E_kz, '.-')
            #ax.axvline(kz_cut, color='k', linestyle='--', alpha=0.7)
            ax.set_xlabel(r'$k_z$')
            ax.set_ylabel(r'$E(k_z)$')
            ax.set_title(r'Vertical spectrum ($z$)')
            ax.set_ylim((10**float(mins[1]), 10**float(maxs[1])))
            ntick=4
            if maxs[1]-mins[1] < 4:
                ntick = int(maxs[1]-mins[1]) + 1
            log_ticks = np.linspace(np.log10(10**float(mins[1])), np.log10(10**float(maxs[1])), ntick)            
            yticks = 10**log_ticks
            ax.yaxis.set_major_locator(ticker.FixedLocator(yticks))
            labels = [rf"$10^{{{int(np.floor(l))}}}$"
                    for l in log_ticks]

            ax.yaxis.set_major_formatter(ticker.FixedFormatter(labels))
            ax.yaxis.set_minor_locator(ticker.NullLocator())

            # --- Bottom-left: horizontal spectrum (kx) ---
            ax = axes[1, 0]
            ax.loglog(np.abs(kx), E_kx, '.-')
            #ax.axvline(kx_cut, color='k', linestyle='--', alpha=0.7)
            ax.set_xlabel(r'$k_x$')
            ax.set_ylabel(r'$E(k_x)$')
            ax.set_title(r'Horizontal spectrum ($x$)')
            ax.set_ylim((10**float(mins[2]), 10**float(maxs[2])))
            ntick=4
            if maxs[2]-mins[2] < 4:
                ntick = int(maxs[2]-mins[2]) + 1
            log_ticks = np.linspace(np.log10(10**float(mins[2])), np.log10(10**float(maxs[2])), ntick)              
            yticks = 10**log_ticks
            ax.yaxis.set_major_locator(ticker.FixedLocator(yticks))
            labels = [rf"$10^{{{int(np.floor(l))}}}$"
                    for l in log_ticks]

            ax.yaxis.set_major_formatter(ticker.FixedFormatter(labels))
            ax.yaxis.set_minor_locator(ticker.NullLocator())

            # --- Bottom-right: 2D planar spectrum (xy) ---
            ax = axes[1, 1]
            ax.loglog(k_xy_centers, E_kxy, '.-')
            #ax.axvline(kxy_cut, color='k', linestyle='--', alpha=0.7)
            ax.set_xlabel(r'$k = \sqrt{k_x^2 + k_y^2}$')
            ax.set_ylabel(r'$E_{xy}(k_{xy})$')
            ax.set_title(r'Horizontal planar spectrum ($xy$)')
            ax.set_ylim((10**float(mins[3]), 10**float(maxs[3])))
            ntick=4
            if maxs[3]-mins[3] < 4:
                ntick = int(maxs[3]-mins[3]) + 1
            log_ticks = np.linspace(np.log10(10**float(mins[3])), np.log10(10**float(maxs[3])), ntick)              
            yticks = 10**log_ticks
            ax.yaxis.set_major_locator(ticker.FixedLocator(yticks))
            labels = [rf"$10^{{{int(np.floor(l))}}}$"
                    for l in log_ticks]

            ax.yaxis.set_major_formatter(ticker.FixedFormatter(labels))
            ax.yaxis.set_minor_locator(ticker.NullLocator())


            # --- Global title with energy check ---
            fig.suptitle(
                f'Temperature Power Spectra \n Time: {time:.4f} \n Energy Check via Parseval: {np.isclose(pars_grid, pars_3d)}'#, grid: {pars_grid:.5e}, 3d: {pars_3d:.5e}, z: {pars_kz:.5e}, x: {pars_kx:.5e}, xy: {pars_kxy:.5e}',
            )

            savename=f"{str(basepath)}/res_check_temp/write_{write:06}.png"
            fig.savefig(savename)
            plt.close()



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

