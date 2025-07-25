'''
Script to perform preliminary analysis tasks on data obtained from running a convection
simulation. Writes total simulation time and Rayleigh and Prandtl numbers to a text file,
and plots the instantaneous Nusselt number and the kinetic energy to help determine when
time averaging should begin. Both a plot over the entire integration time and a zoomed
plot are made.

Usage:
    prelim.py <files>... [--output=<dir>]
    prelim.py <files>...

Options:
    --output=<dir>  Path to analysis output [default: ./prelim]
'''

import h5py # pyright: ignore
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 14})
plt.ioff()
from dedalus.extras import plot_tools

def main(filename, start, count, output):
    with h5py.File(filename, mode='r') as file:
        # get Ra and Pr from analysis file
        Ra = float(file['tasks']['Ra'][-1])
        Pr = float(file['tasks']['Pr'][-1])

        # get time from analysis
        time = np.array(file['scales/sim_time'])
        end_time = time[-1]
        start_time = time[0]

        # get <wT> and KE data
        dset_wT = np.array(file['tasks']['avg_wT'][:])[:,0,0]
        dset_K = np.array(file['tasks']['avg_K'][:])[:,0,0]

        # get instantaneous data
        inst_dset_wT = dset_wT[1:]-dset_wT[:-1]
        inst_dset_K = dset_K[1:]-dset_K[:-1]
        dt_array = time[1:]-time[:-1]
        inst_dset_Nu = 1+ (Ra*Pr)**(1/2)*inst_dset_wT/dt_array
        inst_ave_K = inst_dset_K/dt_array

################################################################################
# generate output of preliminary script
################################################################################

        # write basics to text file
        savename = output.joinpath('prelim.txt')
        info = open(str(savename), 'w')
        info.write('Ra= {:.4e}, Pr= {:.4e}\n'.format(Ra, Pr))
        info.write('End sim time: {:.3f}\n'.format(end_time))
        info.write('Start sim time: {:.3f}\n'.format(start_time))
        info.close()

################################################################################
# plot full plots of instantaneous time averages

        fig=plt.figure(figsize=(50,20))

        ax = fig.add_subplot(2,1,1)
        ax.plot(time[1:], inst_dset_Nu, linewidth=1)
        plt.xlim([time[0],time[-1]])
        plt.title('Instantaneous Nu(t)')
        plt.xlabel('t')
        plt.ylabel('Nu')

        ax = fig.add_subplot(2,1,2)
        ax.plot(time[1:], inst_ave_K, linewidth=1)
        plt.xlim([time[0],time[-1]])
        plt.title('Instantaneous time average of KE')
        plt.xlabel('t')
        plt.ylabel(r'$\frac{1}{100dt}\frac{1}{\Gamma}\int_t^{t+100dt}\int_\Omega u^2+v^2 dxdzd\hat{t}$')

        fig.tight_layout(pad = 1.0)
        savename = output.joinpath('prelim_time_averages.png')
        fig.savefig(str(savename), dpi=400)
        plt.close()

################################################################################
# plot zoomed snapshots of first 1000 time units for Nu and KE
        fig=plt.figure(figsize=(50,20))

        # plot Nusselt number
        ax = fig.add_subplot(2,1,1)
        ax.plot(time[1:], inst_dset_Nu, linewidth=1)
        plt.xlim([0,1000]+start_time)
        plt.title('Instantaneous Nu(t)')
        plt.xlabel('t')
        plt.ylabel('Nu')

        # plot kinetic energy
        ax = fig.add_subplot(2,1,2)
        ax.plot(time[1:], inst_ave_K, linewidth=1)
        plt.xlim([0,1000]+start_time)
        plt.title('Instantaneous time average of KE')
        plt.xlabel('t')
        plt.ylabel(r'$\frac{1}{100dt}\frac{1}{\Gamma}\int_t^{t+100dt}\int_\Omega u^2+v^2 dxdzd\hat{t}$')

        fig.tight_layout(pad = 1.0)
        savename = output.joinpath('prelim_time_averages_zoomed.png')
        fig.savefig(str(savename), dpi=400)
        plt.close()


if __name__ == "__main__":

    import pathlib
    from docopt import docopt
    from dedalus.tools import logging
    from dedalus.tools import post
    from dedalus.tools.parallel import Sync

    args = docopt(__doc__)

    output_path = pathlib.Path(args['--output']).absolute()
    # Create output directory if needed
    with Sync() as sync:
        if sync.comm.rank == 0:
            if not output_path.exists():
                output_path.mkdir()
    post.visit_writes(args['<files>'], main, output=output_path)
