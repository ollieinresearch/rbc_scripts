"""
Script to perform analysis tasks on data obtained from running a 3D (ONLY 3D)
Rayleigh-Benard convection simulation. Accepts an analysis file as a command
line argument and takes time averages from the specified time. All output is
written to the folders 'outputs' and 'preliminary_outputs'. Outputs plots of 
the Nusselt number, various profiles and information about the simulation to a 
text document.

Usage:
    analysis.py <basepath>... [--time=<time>] [--end_time=<end_time>]
    analysis.py <basepath>...

Options:
    --time=<time>   Time at which to begin time averaging [default: 0]
    --end_time=<end_time>   Time at which to end time averaging [default: 0]
"""

import h5py as h5  # pyright: ignore
import numpy as np
from scipy.integrate import cumulative_trapezoid, simpson
from pathlib import Path
import matplotlib
import matplotlib.pyplot as plt
import os
import re
matplotlib.use("Agg")
plt.rcParams.update({"font.size": 14})
plt.ioff()



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



def main(basepath: Path, start_ave: np.float64, end_ave: np.float64):
    """Computes helpful quantities and creates plots for a (hopefully)
    converged flow.

    Args:
        basepath (str | PathLike): The basepath directory.
        start_ave (float): The time at which to begin the analysis of the flow.
    """

    analysis = basepath / "horizontal_analysis"
    output = basepath / 'outputs'
    output.mkdir(exist_ok='true')


    n_secs = 3
    # Get list of all file names in proper numbering
    fs = sorted(
        analysis.glob("horizontal_analysis_s*.h5"),
        key=lambda p: int(re.search(r"horizontal_analysis_s(\d+)\.h5", p.name).group(1))
    )
    fp = fs.pop(0)

    with h5.File(fp, 'r') as f:
        Ra = float(f['tasks']['Ra'][-1])
        Pr = float(f['tasks']['Pr'][-1])

        full_time = f['scales']['sim_time'][:]
        avg_T = f['tasks']['avg_T'][:]
        avg_u_sq = f['tasks']['avg_u_sq'][:]
        avg_w_sq = f['tasks']['avg_w_sq'][:]

        # Currently only 3d, so the try is always accepted.
        try:
            avg_v_sq = f['tasks']['avg_v_sq'][:]
            z_an = f['tasks']['z_an'][0,0,0,:]
            dim = 3
        except:
            z_an = f['tasks']['z_an'][0,0,:]
            dim = 2

        try:
            havg_wT = f['tasks']['havg_wT'][:]
            havg = True
        except:
            havg=False

    
    for fi in fs:
        with h5.File(fi, 'r') as f:
            start_time = full_time[-1]
            if f['scales']['sim_time'][-1] > start_time:
                start_ind = np.searchsorted(f['scales']['sim_time'][:], start_time)

                full_time = np.append(full_time, f['scales']['sim_time'][start_ind:], axis=0)
                avg_T = np.append(avg_T, f['tasks']['avg_T'][start_ind:], axis=0)
                avg_u_sq = np.append(avg_u_sq, f['tasks']['avg_u_sq'][start_ind:], axis=0)
                avg_w_sq = np.append(avg_w_sq, f['tasks']['avg_w_sq'][start_ind:], axis=0)

                # Currently only 3d, so no need to check if this will be in the file or not
                avg_v_sq = np.append(avg_v_sq, f['tasks']['avg_v_sq'][start_ind:], axis=0)
                
                if havg:
                    havg_wT = np.append(havg_wT, f['tasks']['havg_wT'][start_ind:], axis=0)



    increased = full_time[1:] >= full_time[:-1]
    if not np.all(increased):
        print("not all increased")
        indxs = np.argsort(full_time)

        full_time = full_time[indxs] 
        avg_T = avg_T[indxs] 
        avg_u_sq = avg_u_sq[indxs] 
        avg_w_sq = avg_w_sq[indxs] 
        avg_v_sq = avg_v_sq[indxs] 
        
        if havg:
            havg_wT = havg_wT[indxs]
    
    
     
    start_ind = np.searchsorted(full_time, start_ave)
    time = full_time[start_ind:]
    t_0 = time[0]
    t_f = time[-1]
    total_time = t_f - t_0



    ########################################################################
    # Profile plots
    ########################################################################
    z = z_an - 1/2
    # Horizontal temperature profiles, from start_time to end
    if dim == 3:
        avgs = [
            avg_T,
            avg_w_sq,
            avg_u_sq,
            avg_v_sq,
        ]
        dsets = [avg[start_ind:, 0, 0, :] for avg in avgs]

        # Kinetic energy is just the sum of the squared avg velocities
        dsets.append(dsets[1]+dsets[2]+dsets[3])

        horz_tex = r"$\sqrt{\overline{u^2+v^2}}$"
        kin_tex = r"$\overline{u^2+v^2+w^2}$"

    else:
        avgs = [
            avg_T,
            avg_w_sq,
            avg_u_sq,
        ]
        dsets = [avg[start_ind:, 0, :] for avg in avgs]

        # Kinetic energy is just the sum of the squared avg velocities
        dsets.append(dsets[1]+dsets[2])


        horz_tex = r"$\sqrt{\overline{u^2}}$"
        kin_tex = r"$\overline{u^2+w^2}$"


    profs = [simpson(dset, time, axis=0) / total_time for dset in dsets]

    plot_ops = [
        (
            r"T",
            "Horizontally Averaged Temperature Profile"
        ),
        (
            r"$\sqrt{\overline{w^2}}$",
            "Horizontally Averaged RMS Vertical Velocity"
        ),
        (
            horz_tex,
            "Horizontally averaged RMS Horizontal Velocity"
        ),
        (
            kin_tex,
            "Horizontally Averaged Kinetic Energy Profile"
        ),
    ]

    # Plot the calculated profiles
    fig = plt.figure(figsize=(25, 25))
    for ind, (xlabel, title) in enumerate(plot_ops):
        ax = fig.add_subplot(2, 2, ind+1)
        ax.plot(profs[ind], z)
        if ind==0:
            plt.xlim([0,1])
        plt.ylim([-0.5, 0.5])
        plt.xlabel(xlabel)
        plt.ylabel(r"$z$")
        plt.title(title)        

    fig.tight_layout(pad=1.0)
    savename = output.joinpath("profiles.pdf")
    fig.savefig(str(savename), dpi=400)
    plt.close()


    if havg:
        dset = havg_wT[start_ind:, 0,0, :]
        prof = simpson(dset, time, axis=0) / total_time

        # Plot the calculated profiles
        fig = plt.figure(figsize=(12, 12))
        ax = fig.add_subplot()
        ax.plot(prof, z)
        plt.ylim([-0.5, 0.5])
        plt.xlabel(r'$\overline{ \langle wT \rangle}$')
        plt.ylabel(r"$z$")
        plt.title('Horizontally averaged vertical heat transport')        

        fig.tight_layout(pad=1.0)
        savename = output.joinpath("havg_wT.pdf")
        fig.savefig(str(savename), dpi=400)
        plt.close()

################################################################################

if __name__ == "__main__":
    from docopt import docopt

    # Collect arguments
    args = docopt(__doc__)

    basepath = Path(args["<basepath>"][0])
    time = float(args["--time"])
    end_time = float(args["--end_time"])

    main(basepath, time, end_time)
