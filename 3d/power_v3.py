"""
Script to perform analysis tasks on data obtained from running a
Rayleigh-Benard convection simulation. Accepts an analysis file as a command
line argument and takes time averages from the specified time. All output is
written to the (optionally) specified file. Outputs plots of the Nusselt number,
various profiles and information about the simulation to a text document.

Usage:
    power.py <files>... [--ymin=<ymin>] [--ymax=<ymax>]
    power.py <files>...

Options:
    --ymin=<ymin>  Min for y axis, in powers of 10 [default: -12.0]
    --ymax=<ymax>  Max for y axis, in powers of 10 [default: 0.0]
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

def main(file, start, count, ymin=-12.0, ymax=0.0):
    fp = Path(file)
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
            E3d = np.abs(U)**2 + np.abs(V)**2 + np.abs(W)**2

            # Bin power spectrum into magnitude of wavevector
            Eflat = E3d.flatten()
            pk, _ = np.histogram(Knorm_full, bins=bins, weights=Eflat)
            nm, _ = np.histogram(Knorm_full, bins=bins)
            E_k = pk / (nm + eps)

            # Plot
            plt.figure()
            plt.loglog(k_centers, E_k, '.-')
            plt.xlabel('k')
            plt.ylabel('E(k)')
            plt.title(f'E(real)={np.sum(u_uniform**2 + v_uniform**2 + w_uniform**2):.5e}, E(spec)={np.sum(E3d):.5e}, t={time[ti]:.4f}')
            plt.ylim((10**ymin, 10**ymax))
            plt.tight_layout()
            plt.savefig(f"res_check_3d/write_{writes[ti]+c:06}.png")
            plt.close()



if __name__ == "__main__":

    from docopt import docopt
    from dedalus.tools import logging
    from dedalus.tools import post

    args = docopt(__doc__)
    ymin=float(args["--ymin"])
    ymax=float(args["--ymax"])
    post.visit_writes(
        args["<files>"],
        main,
        ymin=ymin,
        ymax=ymax
        
    )

