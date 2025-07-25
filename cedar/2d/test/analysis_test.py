"""
Script to perform analysis tasks on data obtained from running a 2D
Rayleigh-Benard convection simulation. Accepts an analysis file as a command
line argument and takes time averages from the specified time. All output is
written to the (optionally) specified file. Outputs plots of the Nusselt number,
various profiles and information about the simulation to a text document.

Usage:
    analysis.py <files>... [--time=<time>] [--output=<dir>]
    analysis.py <files>...

Options:
    --time=<time>   Time at which to begin time averaging [default: 200]
    --output=<dir>  Path to analysis output [default: ./analysis]
"""

import h5py  # pyright: ignore
import numpy as np
from scipy.integrate import cumulative_trapezoid, simpson
from pathlib import Path
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt

plt.rcParams.update({"font.size": 14})
plt.ioff()

from dedalus.extras import flow_tools


def max_difference(arr: np.ndarray) -> float:
    """Computes the percentage difference between the maximum and minimum
    elements in a numpy ndarray.

    Args:
        arr (np.ndarray): Array of real numbers.

    Returns:
        float: The maximum percentage difference.
    """
    ma, mi = np.max(arr), np.min(arr)

    return (ma - mi) / (0.5 * (ma + mi)) * 100


def main(filename, start, count, output, start_ave):
    """Computes helpful quantities and creates plots for a (hopefully)
    converged flow.

    Args:
        filename (str | PathLike): The file path to the analysis file.
        start (Any): Unused, but required for dedalus' post processing function.
        count (Any): Unused, but required for dedalus' post processing function.
        output (str | PathLike): The output directory.
        start_ave (float): The time at which to begin the analysis of the flow.
    """

    with h5py.File(filename, mode="r") as file:
        # Get time from file, then take data only from those times beyond where
        # averaging begins
        full_time = np.array(file["scales/sim_time"])
        start_ind = np.searchsorted(full_time, start_ave)
        time = full_time[start_ind:]
        t_0 = time[0]
        t_f = time[-1]
        total_time = t_f - t_0

        ########################################################################
        # Get data for profile plots
        ########################################################################

        # Get z-axis points
        z = np.array(file["scales/z/1.0"][:])

        # Horizontal temperature profiles, from start_time to end
        dset_temp = np.array(file["tasks"]["avg_T"][:])[start_ind:, 0, :]
        # Cumulative time average of horizontal profile; add background temp
        temp_prof = simpson(dset_temp, time, axis=0) / total_time
        temp_prof = temp_prof + 1 / 2 - z

        # Horizontal kinetic energy and velocity profiles, from start time to
        # end
        dset_u_sq = np.array(file["tasks"]["avg_u_sq"][:])[start_ind:, 0, :]
        dset_w_sq = np.array(file["tasks"]["avg_w_sq"][:])[start_ind:, 0, :]
        dset_K = dset_u_sq + dset_w_sq
        # Cumulative time average for profile
        K_prof = simpson(dset_K, time, axis=0) / total_time
        u_prof = simpson(dset_u_sq, time, axis=0) / total_time
        w_prof = simpson(dset_w_sq, time, axis=0) / total_time

        # Plot the calculated profiles
        # TODO: for loop this
        fig = plt.figure(figsize=(25, 25))

        # Time-averaged horizontal temperature profile
        ax = fig.add_subplot(2, 2, 1)
        ax.plot(temp_prof, z)
        plt.title("Horizontally Averaged Temperature Profile")
        plt.ylim([-0.5, 0.5])
        plt.xlim([0, 1])
        plt.xlabel(r"$\overline{T}$")
        plt.ylabel(r"$z$")

        # Time-averaged horizontal kinetic energy profile
        ax = fig.add_subplot(2, 2, 2)
        ax.plot(K_prof, z)
        plt.title("Horizontally Averaged Kinetic Energy Profile")
        plt.ylim([-0.5, 0.5])
        plt.xlabel(r"$\overline{u^2+v^2}$")
        plt.ylabel(r"$z$")

        ax = fig.add_subplot(2, 2, 3)
        ax.plot(u_prof, z)
        plt.title("Horizontally averaged RMS Horizontal Velocity")
        plt.xlabel(r"$\overline{u^2}$")
        plt.ylabel("z")
        plt.ylim([-0.5, 0.5])

        ax = fig.add_subplot(2, 2, 4)
        ax.plot(np.sqrt(w_prof), z)
        plt.title("Horizontally averaged RMS Vertical Velocity")
        plt.xlabel(r"$\sqrt{\overline{w^2}}$")
        plt.ylabel("z")
        plt.ylim([-0.5, 0.5])

        fig.tight_layout(pad=1.0)
        savename = output.joinpath("profiles.png")
        fig.savefig(str(savename), dpi=400)
        plt.close()


################################################################################

if __name__ == "__main__":

    from docopt import docopt
    from dedalus.tools import logging
    from dedalus.tools import post
    from dedalus.tools.parallel import Sync

    args = docopt(__doc__)

    output_path = Path(args["--output"]).absolute()
    # Create output directory if needed
    with Sync() as sync:
        if sync.comm.rank == 0:
            output_path.mkdir(parents=True, exist_ok=True)
    post.visit_writes(
        args["<files>"],
        main,
        output=output_path,
        start_ave=float(args["--time"]),
    )
