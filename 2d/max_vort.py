'''
Finds the maximum vorticity in a simulation to correctly set the movie colour 
scale.

Usage:
    prelim.py <files>...
'''

import h5py # pyright: ignore
import numpy as np
from pathlib import Path

def main(filenames):
    
    if n == 0:
        return
    output = Path(filenames[0]).parent / "max_vort.txt"

    # To store maximums
    maxes = np.ones(len(filenames))

    # Check the maximum of each file; the max of this will be used for the limit
    for i, fp in enumerate(filenames):
        with h5py.File(fp, mode='r') as file:
            
            dset = np.array(file['tasks']['vorticity'][:]).ravel()
            maxes[i] = np.max(np.abs(dset))

    mega_max = np.max(maxes)

    # Write the max to a text file
    txt = open(output, 'w')
    txt.write(f"{mega_max}")
    txt.close()



if __name__ == "__main__":

    from docopt import docopt
    

    args = docopt(__doc__)

    main(args['<files>'])