"""
Computes the time average of the power spectra. Requires same resolution between all snapshots files!!!

Usage:
    spectra_avg.py <files>...

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
        snapshots.glob("*_spectra.npz"),
        key=lambda p: int(re.search(r"snapshots_s(\d+)_spectra\.npz$", p.name).group(1))
    )

    
    fp = fs.pop(0)

    # Save the data from the first file; now we may just check whether or not the times are decreasing with each
    # file we load.
    with np.load(fp) as f:
        time = f['time']
        # each spectra in the file is [vert,temp]
        full_spectra = f['full_spectra']
        hor_spectra = f['hor_spectra']
        vert_spectra = f['vert_spectra']
        print(full_spectra.shape[1])
  
    
    # Load all the files, ensuring the time at the end of the new file is more than the last one.
    for fi in fs:
        with np.load(fi) as fi:
            time = np.append(time, fi['time'])
            full_spectra = np.append(full_spectra, fi['full_spectra'], 1)
            hor_spectra = np.append(hor_spectra, fi['hor_spectra'], 1)
            vert_spectra = np.append(vert_spectra, fi['vert_spectra'], 1)


    labels = [
        (full_spectra, r'3D isotropic spectrum', r'$E(k)$', r'$k$'),
        (hor_spectra, r'Horizontal spectrum', r'$E(k_x)$', r'$k_x$'),
        (vert_spectra, r'Vertical spectrum', r'$E(k_z)$', r'$k_z$'),
    ]



    full_time = time[-1] - time[0]
    fig, axes = plt.subplots(2, 2, figsize=(16, 14), layout='constrained')
    with open(output / "spec_range.txt", 'w') as f:
        for k, (spectra, title, y_ax, x_ax) in enumerate(labels):
            j = int(np.floor(k / 2))
            
            for i in range(2):
                quantity = 'velocity' if i == 0 else 'temperature'
                spec = spectra[i]
                x = spec[0,0,:]
                spec = spec[:,1,:]
                #time = time[:spec.shape[0]]
                cumu_spec = simpson(spec, time, axis=0) / full_time
                mask = cumu_spec > 1e-20
                cumu_spec = cumu_spec[mask]
                ma = np.max(np.log10(cumu_spec))
                mi = np.min(np.log10(cumu_spec))
                
                span = ma - mi

                ax = axes[j, k%2]
                ax.loglog(x[mask], cumu_spec, '.-')
                #ax.axvline(k3d_cut, color='k', linestyle='--', alpha=0.7, label='2/3 cutoff')
                ax.set_xlabel(x_ax)
                ax.set_ylabel(y_ax)
                ax.set_title(title)
                ax.hlines([10**ma,10**mi],xmin=np.min(x[mask]), xmax=np.max(x[mask]), color='k', linestyle='--')
                

                f.write(f"Range of spectra for {quantity} {title}: {span:.3f}")
        fig.savefig(output / "cumulative_spectra.pdf", transparent=True, dpi=600)
        plt.close(fig)


if __name__ == "__main__":
    from docopt import docopt

    # Collect arguments
    args = docopt(__doc__)

    file = Path(args["<files>"][0])

    main(file)

