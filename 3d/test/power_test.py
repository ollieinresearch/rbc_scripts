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
    --res=<res>   resolution [default: 120]
    --file=<file>  state file for checkpoint data [default: /state/state_s1.h5]
"""

import h5py as h5
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt


def main(res: int, file: Path):
    z_res = res
    horz_res = 2 * res
    half_res_z = int(z_res/2)

    with h5.File(file, 'r') as f:
        # Get time from file, then take data only from those times beyond where
        # averaging begins
        time = np.array(f["scales/sim_time"])
        if time[0] == 0:
            time = time[1:]
            u_slice = f['tasks']['u'][1:, :, :, half_res_z]
            v_slice = f['tasks']['v'][1:, :, :, half_res_z] 
            w_slice = f['tasks']['w'][1:, :, :, half_res_z] 
        else:
            u_slice = f['tasks']['u'][:, :, :, half_res_z]
            v_slice = f['tasks']['v'][:, :, :, half_res_z]
            w_slice = f['tasks']['w'][:, :, :, half_res_z]

        # Setup frequency bins
        kx = np.fft.fftfreq(horz_res, 1/horz_res)[:, None]
        ky = np.fft.fftfreq(horz_res, 1/horz_res)[None, :]
        k = (kx**2 + ky**2)**0.5
        kmax = int(np.ceil(np.max(k)))
        bins = np.arange(1, kmax+1, 2)
        kcen = bins[:-1] + np.diff(bins)/2

        for i, t in enumerate(time):
            # Use renormalized numpy FFT to compute power spectrum AT time=t
            #TODO: Do i need w? should i look at multiple slices of xy, or should
            # i also look at xz or yz? Do i need the 2*pi*k?
            u_fft = np.fft.fft2(u_slice[i]) / horz_res**2
            v_fft = np.fft.fft2(v_slice[i]) / horz_res**2
            w_fft = np.fft.fft2(w_slice[i]) / horz_res**2
            E_k2 = (np.abs(u_fft)**2 + np.abs(v_fft)**2 + np.abs(w_fft)**2)
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
            plt.ylim((1e-10, 1e0))
            plt.tight_layout()
            plt.savefig(f"res_check/write_{i:06}.png")
            plt.close()


if __name__ == "__main__":

    from docopt import docopt
    from dedalus.tools import logging
    from dedalus.tools import post
    from dedalus.tools.parallel import Sync

    args = docopt(__doc__)

    res = int(args["--res"])
    file = Path(args["--file"])

    main(res, file)
