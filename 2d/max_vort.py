'''
Script to perform preliminary analysis tasks on data obtained from running a convection
simulation. Writes total simulation time and Rayleigh and Prandtl numbers to a text file,
and plots the instantaneous Nusselt number and the kinetic energy to help determine when
time averaging should begin. Both a plot over the entire integration time and a zoomed
plot are made.

Usage:
    prelim.py <files>...
'''

import h5py # pyright: ignore
import numpy as np
from pathlib import Path

def main(filenames):
    n = len(filenames)
    if n == 0:
        return
    output = Path(filenames[0]).parent / "max_vort.txt"

    maxes = np.ones(n)
    for i, fp in enumerate(filenames):
        with h5py.File(fp, mode='r') as file:
            
            dset = np.array(file['tasks']['vorticity'][:]).ravel()
            maxes[i] = np.max(np.abs(dset))

    mega_max = np.max(maxes)

    info = open(output, 'w')
    info.write(f"{mega_max}")
    info.close()



if __name__ == "__main__":

    from docopt import docopt
    

    args = docopt(__doc__)

    main(args['<files>'])