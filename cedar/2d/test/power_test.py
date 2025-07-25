"""
Script to perform analysis tasks on data obtained from running a 2D
Rayleigh-Benard convection simulation. Accepts an analysis file as a command
line argument and takes time averages from the specified time. All output is
written to the (optionally) specified file. Outputs plots of the Nusselt number,
various profiles and information about the simulation to a text document.

Usage:
    power_test.py [--res=<res>] [--file=<file>]
    power_test.py 

Options:
    --res=<res>   resolution [default: 240]
    --file=<file>  state file for checkpoint data [default: /state/state_s1.h5]
"""

import h5py as h5
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt


def main(z_res, file):
    x_res = 2 * z_res
    res = x_res * z_res
    with h5.File(file, 'r') as f:
        # Get time from file, then take data only from those times beyond where
        # averaging begins
        time = np.array(f["scales/sim_time"])
        if np.isclose(time[0],0):
            time = time[1:]
            u = f['tasks']['u'][1:, :, :]
            w = f['tasks']['w'][1:, :, :]
            
        else:
            u = f['tasks']['u'][:, :, :]
            w = f['tasks']['w'][:, :, :]

        # Setup frequency bins
        kx = np.fft.fftfreq(x_res, 1/x_res)[:, None]
        kz = np.fft.fftfreq(z_res, 1/z_res)[None, :]
        k = (kx**2 + kz**2)**0.5
        kmax = int(np.ceil(np.max(k)))
        bins = np.arange(1, kmax+1, 2)
        kcen = bins[:-1] + np.diff(bins)/2

        for i, t in enumerate(time):
            # Use renormalized numpy FFT to compute power spectrum AT time=t

            #TODO: Do i want w? - yes should i look at multiple slices of xy, or should
            # i also look at xz or yz? - no Do i need the 2*pi*k? or the other factors
            # from the tutorial? - no are the resolutions in
            # the right place for different size resolutions? are there any other
            # factors that i need to include?
            # Unrelated, how do i change the scale for snapshots? for 2d, can i run
            # with nodes = vertical resolution (1/2 horiz res)? or is this the hard max
            # of 128? - do some tests

            u_fft = np.fft.fft2(u[i]) / res
            w_fft = np.fft.fft2(w[i]) / res
            E_k2 = (np.abs(u_fft)**2 + np.abs(w_fft)**2)
            E_k1 = E_k2 * 2 * np.pi * k
            # Build histogram over modes, weighted by energy
            pow_samples, _ = np.histogram(k, bins=bins, weights=E_k1)
            hist_samples, _ = np.histogram(k, bins=bins)
            spectrum = pow_samples / hist_samples

            # Plot histogram
            plt.figure()
            plt.loglog(kcen, spectrum, '.-')
            plt.xlabel("k")
            plt.ylabel("E(k)")
            plt.title(f't={t:.4f}')
            plt.ylim((1e-12, 1e-1))
            plt.tight_layout()
            plt.savefig(f"res_check/write_{i:06}.png")
            plt.close()


if __name__ == "__main__":

    from docopt import docopt
    from dedalus.tools import logging
    from dedalus.tools import post
    from dedalus.tools.parallel import Sync

    args = docopt(__doc__)

    z_res = int(args["--res"])
    file = Path(args["--file"])

    main(z_res, file)