"""
Checks the magnitude of the vorticity of snapshot files for negative values and volume averages. REQUIRES SAME RESOLUTION BETWEEN SNAPSHOT FILES!

Usage:
    v_analysis.py <files>...

"""



import h5py as h5
import numpy as np
from pathlib import Path
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("Agg")
plt.rcParams.update({"font.size": 14})
plt.ioff()
from scipy.integrate import simpson, cumulative_trapezoid
import os
import re



def main(basepath):
    snapshots = basepath / "snapshots"
    output = basepath / 'outputs'
    output.mkdir(exist_ok='true')


    # We need the snapshot files in the right order; sometimes portions of simulations are accidentally ran multiple
    # times, and in that case we may only consider one of the outputs (otherwise time integration can get a bit wacky).
    # To only consider one set of outputs, we must ensure the times are nondecreasing, which requires the right
    # file order.
    fs = sorted(
        snapshots.glob("*_avgd.npz"),
        key=lambda p: int(re.search(r"snapshots_s(\d+)_avgd\.npz$", p.name).group(1))
    )

    
    fp = fs.pop(0)

    # Save the data from the first file; now we may just check whether or not the times are decreasing with each
    # file we load.
    with np.load(fp) as f:
        time = f['time']
        num = f['num']
        omega = f['omega']
        omega_full = f['omega_full']
  
    
    # Load all the files, ensuring the time at the end of the new file is more than the last one.
    for fi in fs:
        with np.load(fi) as fi:
            time = np.append(time, fi['time'])
            num = np.append(num, fi['num'])
            omega = np.append(omega, fi['omega'])
            omega_full = np.append(omega_full, fi['omega_full'])
    
    fig, axes = plt.subplots(
            nrows=2,
            ncols=1,
            figsize=(8, 3 * 3),
            layout="constrained",
        )
    
    axes[0].plot(time, num)
    #axes[0].set_ylabel(r"Fraction of grid points with $\omega^2 < 0$")
    axes[0].set_xlabel(r"$t$")
    #axes[0].tick_params(axis='y', labelrotation=90)

    
    axes[1].plot(time, -1 * omega / omega_full)
    #axes[1].set_title(r"Normalized volume average of $\min\{\omega^2, 0\}$ over time")
    axes[1].set_xlabel(r"$t$")
    #axes[1].set_ylabel(r"$\frac{\int_\Omega \min\{\omega^2, 0\} dx}{\int_\Omega \omega^2 dx}$")
    #axes[1].tick_params(axis='y', labelrotation=90)
    fig.savefig(output / "vort_avg.pdf", dpi=400)
    plt.close(fig)

if __name__ == "__main__":
    from docopt import docopt

    # Collect arguments
    args = docopt(__doc__)

    file = Path(args["<files>"][0])

    main(file)

