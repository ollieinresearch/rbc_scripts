'''
Tiny script to find the index of the restart file we want to load from.
Change the value of load_time to get different spots in an existing simulation.
The printed value ic_index should be passed as index into the RB script.
Load time of 0 implies the most recent index should be used.

Usage:
    initial_condition.py <time> --file=<file>

Options:
    --file=<file> Path to the restart file.
'''
import h5py # pyright: ignore
import numpy as np
from docopt import docopt

# Collect arguments
args = docopt(__doc__)
load_time = float(args['<time>'])
fp = str(args['--file'])

# Open file
with h5py.File(fp, mode='r') as file:
    # Load the time values
    time = np.array(file['scales/sim_time'])

    # Find the closest time to the specified load_time
    if load_time == 0:
        ic_index = len(time)-1
    else:
        ic_index = np.argmin(np.abs(time-load_time))

    # Print out for the shell script to collect
    print(ic_index)
    print(time[ic_index])
