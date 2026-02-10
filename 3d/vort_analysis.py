"""
Checks the magnitude of the vorticity of snapshot files for negative values and volume averages. REQUIRES SAME RESOLUTION BETWEEN SNAPSHOT FILES!

Usage:
    vort_analysis.py <files>...

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
        snapshots.glob("snapshots_s*.h5"),
        key=lambda p: int(re.search(r"snapshots_s(\d+)\.h5", p.name).group(1))
    )

    print(len(fs))
    
    fp = fs.pop(0)

    # Save the data from the first file; now we may just check whether or not the times are decreasing with each
    # file we load.
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

    # Load all the files, ensuring the time at the end of the new file is more than the last one.
    for fi in fs:
        with h5.File(fi, 'r') as f:
            start_time = time[-1]
            if f['scales']['sim_time'][-1] > start_time:
                start_ind = np.searchsorted(f['scales']['sim_time'][:], start_time)

                time = np.append(time, f['scales']['sim_time'][start_ind:], axis=0)
                writes = np.append(writes, f["scales/write_number"][start_ind:], axis=0)
                omega = np.append(omega, f['tasks']['vorticity'][start_ind:], axis=0)
                print(time.shape, writes.shape, omega.shape)

                
                # Collect the file data (start/count is for parallelization)

    print(time.shape, writes.shape, omega.shape)
    increased = time[1:] >= time[:-1]
    if not np.all(increased):
        print("not all increased")
        indxs = np.argsort(time)

        time = time[indxs] 
        writes = writes[indxs] 
        omega = omega[indxs] 


    
    fig, axes = plt.subplots(
            nrows=3,
            ncols=1,
            figsize=(12, 3 * 3),
            layout="constrained",
        )
    
    c=0
    omeg_mask = omega < 0
    num_neg = np.sum(omeg_mask, axis=(1,2,3))
    num_tot = x.shape[0] * y.shape[0] * z.shape[0]

    omega = np.where(omeg_mask, omega, 0)
    
    dset_omega = simpson(
        simpson(
            simpson(
                omega, x, axis=1
            ), y, axis=1
        ), z, axis=1
    )

    print(dset_omega.shape, time.shape)

    cumu_omega = cumulative_trapezoid(dset_omega, time, axis=0) / (time[1:]-time[0])

    axes[0].plot(time, num_neg/num_tot)
    axes[0].set_title(r"Percentage of grid points with a negative $\omega^2$")
    axes[0].set_xlabel(r"$t$")
    axes[1].plot(time, dset_omega)
    axes[1].set_title(r"Volume average of $\min\{\omega^2, 0\}$ over time")
    axes[1].set_xlabel(r"$t$")
    axes[1].set_ylabel(r"$\int_\Omega \min\{\omega^2, 0\} dx$")
    axes[2].plot(time[1:], cumu_omega)
    axes[2].set_title(r"$\langle\overline{\min\{\omega^2, 0\}}\rangle$")
    axes[2].set_xlabel(r"$t$")

    fig.savefig(output / "vort_avg.pdf", dpi=400)


if __name__ == "__main__":
    from docopt import docopt

    # Collect arguments
    args = docopt(__doc__)

    file = Path(args["<files>"][0])

    main(file)

