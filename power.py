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

#TODO: add the y limits as arguments

import h5py as h5
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt


def main(file, start, count, ylims):
    fp = Path(file)
    with h5.File(fp, 'r') as f:

        time = np.array(f["scales/sim_time"][start:start+count])
        writes = np.array(f["scales/write_number"][start:start+count])
        u = np.array(f['tasks']['u'][start:start+count])
        w = np.array(f['tasks']['w'][start:start+count])
        shp = u.shape
        dim = len(shp) - 1
        mid_z = 0
        if dim == 3:
            v = np.array(f['tasks']['v'][start:start+count])

            nt, nx, ny, nz = shp
            horz_res = nx
            res = horz_res**2

        else:
            nt, nx, nz = shp
            horz_res = nx
            vert_res = nz
            res = horz_res * vert_res
            mid_z = vert_res / 2            
        
        c=0
        if np.isclose(time[0],0):
            time = time[1:]
            u = u[1:]
            w = w[1:]
            if dim == 3:
                v = v[1:]
            count -= 1
            c = 1

        


        # Setup frequency bins
        k1 = np.fft.fftfreq(horz_res, 1/horz_res)[:, None]
        k2 = np.fft.fftfreq(horz_res, 1/horz_res)[None, :] if dim == 3 else np.fft.fftfreq(vert_res, 1/vert_res)[None, :]

        # Wavenumbers
        k = np.sqrt((k1**2 + k2**2))
        kmax = int(np.ceil(np.max(k)))
        bins = np.arange(1, kmax+1, 2)
        kcen = bins[:-1] + np.diff(bins)/2

        for i in range(count):
            # Use renormalized numpy FFT to compute power spectrum AT time=t

            #TODO: Do i want w? - yes should i look at multiple slices of xy, or should
            # i also look at xz or yz? - no Do i need the 2*pi*k? -yes  or the other factors
            # from the tutorial? - no are the resolutions in
            # the right place for different size resolutions? are there any other
            # factors that i need to include?
            # Unrelated, how do i change the scale for snapshots? for 2d, can i run
            # with nodes = vertical resolution (1/2 horiz res)? or is this the hard max
            # of 128? - do some tests

            # For 2d, there is only one 2d slice to use
            slices = [u[i], w[i]] if dim ==2 else [e[i, :, :, mid_z] for e in [u,v,w]]
            ffts = np.power(np.abs(np.fft.fft2(slices) / res), 2)
            
            E_k2 = np.sum(ffts, axis=0)
            E_k1 = E_k2 * 2 * np.pi * k
            # Build histogram over modes, weighted by energy
            # FROM DAVID using like a finite difference of the energy/T(k)^2
            pow_samples, _ = np.histogram(k, bins=bins, weights=E_k1)
            hist_samples, _ = np.histogram(k, bins=bins)
            spectrum = pow_samples / hist_samples

            # Plot histogram
            plt.figure()
            plt.loglog(kcen, spectrum, '.-')
            plt.xlabel("k")
            plt.ylabel("E(k)")
            plt.title(f't={time[i]:.4f}')
            plt.ylim((10**ylims[0], 10**ylims[1]))
            plt.tight_layout()
            plt.savefig(f"res_check/write_{writes[i]+c:06}.png")
            plt.close()


if __name__ == "__main__":

    from docopt import docopt
    from dedalus.tools import logging
    from dedalus.tools import post

    args = docopt(__doc__)
    ylims = (args["--ymin"], args["--ymax"])

    post.visit_writes(
        args["<files>"],
        main,
        ylims
    )