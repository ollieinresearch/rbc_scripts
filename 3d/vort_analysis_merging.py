"""
Script to perform analysis tasks...

Usage:
    power_vort.py <files>... [--mins=<mins>] [--maxs=<maxs>]

Options:
    --mins=<mins>   Minimum log10 limit. String of csv floats.
    --maxs=<maxs>   Maximum log10 limit. String of csv floats.
"""



import h5py as h5
import numpy as np
from pathlib import Path
from scipy.integrate import simpson, cumulative_trapezoid

#TODO: get this working for temp as well!
# Questions:
# Should I be checking that the energy is the same as the energy computed using the analysis files?
# Any extra multiplicative factors that I forgot? does it really matter if things are scaled?

def main(file, start, count):
    fp = Path(file)
    basepath = Path(fp.parents[1])
    with h5.File(fp, 'r') as f:
        time = np.array(f["scales/sim_time"][:])
        writes = np.array(f["scales/write_number"][:])
        omega = np.array(f['tasks']['vorticity'][:])
        scales = f["scales"]

        x_key = next(k for k in scales.keys() if k.startswith("x_"))
        y_key = next(k for k in scales.keys() if k.startswith("y_"))
        z_key = next(k for k in scales.keys() if k.startswith("z_"))

        x = np.array(scales[x_key])
        y = np.array(scales[y_key])
        z = np.array(scales[z_key])       
    

    dset_omega_full = simpson(
        simpson(
            simpson(
                omega, x, axis=1
            ), y, axis=1
        ), z, axis=1
    )


    omeg_mask = omega < 0
    omega = np.where(omeg_mask, omega, 0)
    num_neg = np.sum(omeg_mask, axis=(1,2,3))
    num_tot = x.shape[0] * y.shape[0] * z.shape[0]
    num_neg = num_neg / num_tot
    
    dset_omega = simpson(
        simpson(
            simpson(
                omega, x, axis=1
            ), y, axis=1
        ), z, axis=1
    )

    np.savez(fp.parent / f"{fp.stem}_avgd.npz", time=time, num=num_neg, omega=dset_omega, omega_full=dset_omega_full)




if __name__ == "__main__":

    from docopt import docopt
    from dedalus.tools import logging
    from dedalus.tools import post

    args = docopt(__doc__)


    post.visit_writes(
        args["<files>"],
        main
    )

