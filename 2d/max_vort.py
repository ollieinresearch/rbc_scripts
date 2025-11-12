'''
Finds the maximum vorticity in a simulation to correctly set the movie colour 
scale.

Usage:
    max_vort.py <files>...
'''

import h5py # pyright: ignore
import numpy as np
from pathlib import Path

def main(filename, start, count):
    
    basepath = Path(filename).parent

    # Check the maximum of each chunk; the max of the maxes will be the limit
    with h5py.File(filename, mode='r') as file:
        
        dset = np.array(file['tasks']['vorticity'][start:start+count])
        mv = np.max(np.abs(dset))

    # Write the max to a text file
    with open(basepath / f"max_vort_{start}.txt", 'w') as txt:
        txt.write(f"{mv}")



if __name__ == "__main__":

    from docopt import docopt
    from dedalus.tools import post


    args = docopt(__doc__)

    post.visit_writes(args['<files>'], main)


