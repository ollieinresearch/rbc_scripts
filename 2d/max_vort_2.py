import numpy as np
from pathlib import Path

basepath = Path(args['<files>'][0]).parent

txts = list(basepath.glob("*.txt"))

maxes = np.ones(len(txts))
for i, f in enumerate(txts):
    with f.open() as txt:
        val = float(txt.read())
        print(val)
        maxes[i] = val
    f.unlink()

print(maxes)
mega_max = np.max(maxes)

final = basepath / "max_vort.txt"
with open(final, 'w') as txt:
    txt.write(f"{mega_max}")