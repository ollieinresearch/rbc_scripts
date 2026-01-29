'''
Finds the maximum vorticity in a simulation to correctly set the movie colour 
scale.

Usage:
    max_vort_2.py <base>...
'''



import numpy as np
from pathlib import Path
from docopt import docopt

args = docopt(__doc__)

basepath = Path(args['<base>'][0])

txts = list(basepath.glob("*.txt"))
maxes = np.ones(len(txts))
for i, f in enumerate(txts):
    with f.open() as txt:
        val = float(txt.read())
        maxes[i] = val
    f.unlink()

mega_max = np.max(maxes)

final = basepath / "max_vort.txt"
with open(final, 'w') as txt:
    txt.write(f"{mega_max}")