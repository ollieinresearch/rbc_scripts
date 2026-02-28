"""
Script to perform analysis tasks...

Usage:
    power_combined.py <files>... [--vmins=<vmins>] [--vmaxs=<vmaxs>] [--tmins=<tmins>] [--tmaxs=<tmaxs>]

Options:
    --vmins=<vmins>   Minimum log10 limit. String of csv floats.
    --vmaxs=<vmaxs>   Maximum log10 limit. String of csv floats.
    --tmins=<tmins>   Minimum log10 limit. String of csv floats.
    --tmaxs=<tmaxs>   Maximum log10 limit. String of csv floats.
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

#TODO: get this working for temp as well!
# Questions:
# Should I be checking that the energy is the same as the energy computed using the analysis files?
# Any extra multiplicative factors that I forgot? does it really matter if things are scaled?

def main(file, start, count, vmins, vmaxs, tmins, tmaxs):
    fp = Path(file)
    basepath = Path(fp.parents[1])

    prev_nx, prev_ny, prev_nz = (0,0,0)

    with h5.File(fp, 'r') as f:
        times=[]
        full_spectra = [[], []]
        hor_spectra = [[], []]
        vert_spectra = [[], []]
        nt = f["scales/sim_time"].shape[0]

        for ti in range(nt):
            # Collect the state file data (start/count is for parallelization)
            time = f["scales/sim_time"][ti]
            write = f["scales/write_number"][ti]
            temp = np.array(f['tasks']['temperature'][ti])


            try:
                u = np.array(f['tasks']['u'][ti, 0])
                v = np.array(f['tasks']['u'][ti, 1])
                w = np.array(f['tasks']['u'][ti, 2])
                nx, ny, nz = u.shape

            except:
                try:
                    u = np.array(f['tasks']['x'][ti])
                    v = np.array(f['tasks']['y'][ti])
                    w = np.array(f['tasks']['w'][ti])
                    nx, ny, nz = u.shape

                except:
                    u = np.array(f['tasks']['u'][ti])
                    v = np.array(f['tasks']['v'][ti])
                    w = np.array(f['tasks']['w'][ti])
                    nx, ny, nz = u.shape

            
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

            # Create interpolators for current time slice,
            # Interpolate onto uniform z grid,
            # Compute 3D FFT and energy. Supposedly the ortho norm negates the
            # need for renormalization. No dividing by resolution :)

            interp_u = RegularGridInterpolator((x, y, z_cheb), u, method='linear', bounds_error=False, fill_value=0)
            u_uniform = interp_u(pts).reshape((nx, ny, nz))
            pars_grid = np.sum(u_uniform**2)
            U = np.fft.fftn(u_uniform, norm='ortho')
            E_3d = np.abs(U)**2
            del interp_u, u_uniform, U

            interp_v = RegularGridInterpolator((x, y, z_cheb), v, method='linear', bounds_error=False, fill_value=0)
            v_uniform = interp_v(pts).reshape((nx, ny, nz))
            pars_grid += np.sum(v_uniform**2)
            V = np.fft.fftn(v_uniform, norm='ortho')
            E_3d += np.abs(V)**2
            del interp_v, v_uniform, V

            interp_w = RegularGridInterpolator((x, y, z_cheb), w, method='linear', bounds_error=False, fill_value=0)
            w_uniform = interp_w(pts).reshape((nx, ny, nz))
            pars_grid += np.sum(w_uniform**2)
            W = np.fft.fftn(w_uniform, norm='ortho')
            E_3d += np.abs(W)**2
            del interp_w, w_uniform, W
            

            # Bin power spectrum into magnitude of 3d wavevector - Full 3d spectra
            E_k_flat = E_3d.flatten()
            pk, _ = np.histogram(Knorm_full, bins=bins, weights=E_k_flat)
            nm, _ = np.histogram(Knorm_full, bins=bins)
            # safe division
            with np.errstate(invalid='ignore', divide='ignore'):
                E_k = pk / nm
            mask = nm > 0
            k_centers_plot = k_centers[mask]
            E_k_plot = E_k[mask]

            E_kz_raw = np.sum(E_3d, axis=(0,1))   # one value per kz index
            kz_vals = np.abs(kz)                  # positive-frequency values (may repeat)
            # Aggregate identical |kz| values (e.g., +k and -k if present)
            kz_unique, inv = np.unique(kz_vals, return_inverse=True)
            E_kz_modes = np.bincount(inv, weights=E_kz_raw)
            counts_kz = np.bincount(inv)
            E_kz_modes_avg = E_kz_modes / (counts_kz + eps)

            # same for kx:
            E_kx_raw = np.sum(E_3d, axis=(1,2))
            kx_vals = np.abs(kx)
            kx_unique, invx = np.unique(kx_vals, return_inverse=True)
            E_kx_modes = np.bincount(invx, weights=E_kx_raw)
            counts_kx = np.bincount(invx)
            E_kx_modes_avg = E_kx_modes / (counts_kx + eps)


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

            fig, axes = plt.subplots(2, 2, figsize=(16, 14), layout='constrained')
            maxs = vmaxs
            mins = vmins
            # --- Dealiasing cutoffs (2/3 rule) ---  IDK IF THIS IS USEFUL
            #kx_cut = (2/3) * np.max(np.abs(kx))
            #ky_cut = (2/3) * np.max(np.abs(ky))
            #kz_cut = (2/3) * np.max(np.abs(kz))

            # For isotropic 3D and planar spectra, use the smallest relevant cutoff
            #k3d_cut = min(kx_cut, ky_cut, kz_cut)
            #kxy_cut = min(kx_cut, ky_cut)

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
            idx_z = np.argsort(kz_unique)
            ax.loglog(kz_unique[idx_z], E_kz_modes_avg[idx_z], '.-')
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
            idx_x = np.argsort(kx_unique)
            ax.loglog(kx_unique[idx_x], E_kx_modes_avg[idx_x], '.-')
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
            ax.set_title(rf'Horizontal planar spectrum ($xy$)')
            ax.set_ylim((10**float(mins[3]), 10**float(maxs[3])))
            ntick = 4
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
                f'Velocity Power Spectra \n Time: {time:.4f} \n Energy Check via Parseval: {np.isclose(pars_grid, pars_3d)}'#, grid: {pars_grid:.5e}, 3d: {pars_3d:.5e}, z: {pars_kz:.5e}, x: {pars_kx:.5e}, xy: {pars_kxy:.5e}',
            )

            savename=f"{str(basepath)}/res_check_3d/write_{write:06}.png"
            fig.savefig(savename)
            plt.close(fig)

            full_spectra[0].append(np.array([k_centers, E_k]))
            hor_spectra[0].append(np.array([kx_unique[idx_x], E_kx_modes_avg[idx_x]]))
            vert_spectra[0].append(np.array([kz_unique[idx_z], E_kz_modes_avg[idx_z]]))



            #######################################################
            # temp
            #######################################################

            # Main loop: interpolate, FFT, bin, plot
            # Interpolators for current time slice
            interp_temp = RegularGridInterpolator((x, y, z_cheb), temp, method='linear', bounds_error=False, fill_value=0)

            # Interpolate onto uniform z grid
            temp_uniform = interp_temp(pts).reshape((nx, ny, nz))
            # pars_grid = np.sum(temp_uniform**2)

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


            E_kz_raw = np.sum(E_3d, axis=(0,1))   # one value per kz index
            kz_vals = np.abs(kz)                  # positive-frequency values (may repeat)
            # Aggregate identical |kz| values (e.g., +k and -k if present)
            kz_unique, inv = np.unique(kz_vals, return_inverse=True)
            E_kz_modes = np.bincount(inv, weights=E_kz_raw)
            counts_kz = np.bincount(inv)
            E_kz_modes_avg = E_kz_modes / (counts_kz + eps)

            # same for kx:
            E_kx_raw = np.sum(E_3d, axis=(1,2))
            kx_vals = np.abs(kx)
            kx_unique, invx = np.unique(kx_vals, return_inverse=True)
            E_kx_modes = np.bincount(invx, weights=E_kx_raw)
            counts_kx = np.bincount(invx)
            E_kx_modes_avg = E_kx_modes / (counts_kx + eps)

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
            maxs = tmaxs
            mins = tmins
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
            idx_z = np.argsort(kz_unique)
            ax.loglog(kz_unique[idx_z], E_kz_modes_avg[idx_z], '.-')
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
            idx_x = np.argsort(kx_unique)
            ax.loglog(kx_unique[idx_x], E_kx_modes_avg[idx_x], '.-')
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
            plt.close(fig)

            full_spectra[1].append(np.array([k_centers, E_k]))
            hor_spectra[1].append(np.array([kx_unique[idx_x], E_kx_modes_avg[idx_x]]))
            vert_spectra[1].append(np.array([kz_unique[idx_z], E_kz_modes_avg[idx_z]]))

            times.append(time)

        # We have now saved the temp and velocity spectra for each time for time avg; write to storage
        full_spectra = np.array(full_spectra)
        hor_spectra = np.array(hor_spectra)
        vert_spectra = np.array(vert_spectra)
        np.savez(fp.parent / f"{fp.stem}_spectra.npz", time=times,full_spectra=full_spectra, hor_spectra=hor_spectra, vert_spectra=vert_spectra)





if __name__ == "__main__":

    from docopt import docopt
    from dedalus.tools import logging
    from dedalus.tools import post

    args = docopt(__doc__)
    vmins = list(map(float, args["--vmins"].split(",")))
    vmaxs = list(map(float, args["--vmaxs"].split(",")))
    tmins = list(map(float, args["--tmins"].split(",")))
    tmaxs = list(map(float, args["--tmaxs"].split(",")))



    post.visit_writes(
        args["<files>"],
        main,
        vmins=vmins,
        vmaxs=vmaxs,
        tmins=tmins,
        tmaxs=tmaxs
    )

