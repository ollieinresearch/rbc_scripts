"""
Plot planes from snapshot files to produce frames for visualization. Accepts a
list of .h5 files as a command line argument. Adapted from the plot_slices.py
sample script found on the dedalus github page.

Note: if you have too many writes to plot, change the modulus on line 55 (set to
1 by default but comes in handy if you're working at low Pr and forgot to change
how often snapshots get saved. # TODO: add this as a command line option.


Usage:
    plot_slices.py <files>... [--output=<dir>] [--max_vort=<max_vort>] [--plot_freq=<plot_freq>]

Options:
    --output=<dir>           Output directory [default: ./frames]
    --max_vort=<max_vort>    Maximum achieved vorticity [default: 10000]
    --plot_freq=<plot_freq>  How often to plot a snapshot [default: 1]

"""

from matplotlib.scale import AsinhTransform
import h5py # pyright: ignore
from copy import deepcopy
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 14})
plt.ioff()
from dedalus.extras import plot_tools



asinh = AsinhTransform(linear_width=1.5)


def sym(xmesh, ymesh, data):
    """Function to transform scalar values for colorbar"""

    og_shape = data.shape
    new_dset = data.ravel()
    new_dset = asinh.transform(new_dset)
    data = np.reshape(new_dset, og_shape)
    
    return xmesh, ymesh, data

def main(filename, start, count, output, max_vort, plot_freq):
    """Save plot of specified tasks for given range of analysis writes."""
    # fields to plot. Change as needed.
    tasks = ['temp', 'vorticity']
    filename = Path(filename)
    
    
    # Plot settings
    scale = 2.5
    dpi = 400
    title_func = lambda sim_time: 't = {:.3f}'.format(sim_time)
    savename_func = lambda write: 'write_{:06}.png'.format(write)
    # Layout
    nrows, ncols = 2, 1
    # aspect ratio
    image = plot_tools.Box(2, 1)
    # padding per image
    pad = plot_tools.Frame(0.2, 0.2, 0.1, 0.1)
    margin = plot_tools.Frame(0.4, 0.2, 0.1, 0.1)

    # Create multifigure
    mfig = plot_tools.MultiFigure(nrows, ncols, image, pad, margin, scale)
    fig = mfig.figure


    # Plot writes
    with h5py.File(filename, mode='r') as file:
        for index in range(start, start+count):
            # dummy check; if you have too many writes to plot, change the modulus
            # to only plot certain writes; happens at low Pr
            if file['scales/write_number'][index]%plot_freq==0:
                for n, task in enumerate(tasks):
                    # Build subfigure axes
                    i, j = divmod(n, ncols)
                    axes = mfig.add_axes(i, j, [0, 0, 1, 1])
                    # Call 3D plotting helper, slicing in time
                    dset = file['tasks'][task]
                    if task == 'vorticity':
                        new_mv = 0.95*asinh.transform(max_vort)
                        plot_tools.plot_bot_3d(dset, 0, index, axes=axes, title=task, clim=(-new_mv, new_mv), even_scale=True, func=sym)
                    else:
                        plot_tools.plot_bot_3d(dset, 0, index, axes=axes, title=task)
                # Add time title
                title = title_func(file['scales/sim_time'][index])
                title_height = 1 - 0.5 * mfig.margin.top / mfig.fig.y
                fig.suptitle(title, x=0.48, y=title_height, ha='left')
                # Save figure
                savename = savename_func(file['scales/write_number'][index])
                savepath = output.joinpath(savename)
                fig.savefig(str(savepath), dpi=dpi)
                fig.clear()
        plt.close(fig)


if __name__ == "__main__":

    from pathlib import Path
    from docopt import docopt
    from dedalus.tools import logging
    from dedalus.tools import post
    from dedalus.tools.parallel import Sync

    args = docopt(__doc__)

    output_path = Path(args['--output']).absolute()
    max_vort = np.float64(args['--max_vort'])
    plot_freq = int(args['--plot_freq'])

    # Create output directory if needed
    with Sync() as sync:
        if sync.comm.rank == 0:
            if not output_path.exists():
                output_path.mkdir()
    post.visit_writes(args['<files>'], main, output=output_path, max_vort=max_vort, plot_freq=plot_freq)
