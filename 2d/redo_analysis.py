"""
Script to perform analysis tasks on data obtained from running a 2D
Rayleigh-Benard convection simulation. Accepts an analysis file as a command
line argument and takes time averages from the specified time. All output is
written to the (optionally)specified file. Outputs plots of the Nusselt number,
various profiles and info about the simulation to a text document.
"""

#TODO: update this with the new analysis.py

import h5py  # pyright: ignore
import numpy as np
import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 14})
plt.ioff()

from dedalus.extras import flow_tools



def normalize_exponent(s):
    return re.sub(r'e([-+]?)(0*)(\d*)', lambda m: f"e{m.group(1)}{m.group(3) if m.group(3).lstrip('0') else '0'}", s)



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


def main(filename, start, count, start_ave,):
    """Computes helpful quantities and creates plots for a (hopefully)
    converged flow.

    Args:
        filename (str | PathLike): The file path to the analysis file.
        start (Any): Unused, but required for dedalus' post processing function.
        count (Any): Unused, but required for dedalus' post processing function.
        output (str | PathLike): The output directory.
        start_ave (float): The time at which to begin the analysis of the flow.
    """

    # From the filename, extract Ra, Pr, and resolution.
    filename = Path(filename)
    fp = filename.parent.stem

    match = re.match(r'.*x(.*)_.*_ra(.*)_pr(.*)_.*', fp)
    
    res = match.group(1)
    ra = match.group(2)
    ra = normalize_exponent(ra)
    pr = match.group(3)
    pr = normalize_exponent(pr)

    output = Path("/project/def-goluskin/ollie/ollie_rb_data/redoing_test")
    output.mkdir(parents=True, exist_ok=True)
    output = output / f"ra{ra}/pr{pr}/{res}_Gam2"
    output.mkdir(parents=True, exist_ok=True)

    # Rayleigh and Prandtl number from the file
    Ra = float(ra)
    Pr = float(pr)
    
    with h5py.File(filename, mode="r") as file:        
        # Get time from file, then take data only from those times beyond where
        # averaging begins
        time = np.array(file["scales/sim_time"])
        start_average = start_ave[str(filename)]
        start_ind = np.searchsorted(time, start_average)
        time = time[start_ind:]
        dt_array = np.diff(time)
        t_0 = time[0]
        t_f = time[-1]
        total_time = t_f - t_0

        # Write the basics to a text file
        savename = output / "info.txt"
        info = open(savename, "w")
        info.write(f"Ra = {Ra:.4e}, Pr = {Pr:.4e}\n")
        info.write(f"Simulation end time: {t_f:.4f}\n")
        info.write(f"Time averaging begins at: {t_0:.4f}\n")
        info.write(
            "Below are calculations of Nu, using various quantities of "
            "interest (QoI).\n"
        )
        info.write("-" * 80)
        info.write("\n")

        # Functions to calculate Nu with different quantities
        def fnu_wT(dset: np.ndarray, time_array) -> np.ndarray:
            """Calculate Nu using the <w*T> formula.

            Args:
                dset (np.ndarray): The w*T quantity from analysis file.
                time_array (np.ndarray): The time array for calculating time
                averages.

            Returns:
                np.ndarray: The Nusselt number throughout time.
            """

            return 1 + ((Ra * Pr) ** (1 / 2) * dset) / time_array

        def fnu_wtheta(dset: np.ndarray, time_array) -> np.ndarray:
            """Calculate Nu using the <w*θ> formula."""

            return 1 + ((Ra * Pr) ** (1 / 2) * dset) / time_array

        def fnu_grad_u_sq(dset: np.ndarray, time_array) -> np.ndarray:
            """Calculate Nu using the <|∇u|^2> formula."""

            return 1 + Pr * dset / time_array

        def fnu_oy_sq(dset: np.ndarray, time_array) -> np.ndarray:
            """Calculate Nu using the <oy^2> formula."""

            return 1 + Pr * dset / time_array

        def fnu_grad_T_sq(dset: np.ndarray, time_array) -> np.ndarray:
            """Calculate Nu using the <|∇T|^2> formula."""

            return dset / time_array

        def fnu_grad_theta_sq(dset: np.ndarray, time_array) -> np.ndarray:
            """Calculate Nu using the <|θ|^2> formula."""

            return 1 + dset / time_array


        def fre_wt(dset: np.ndarray, time_array) -> np.ndarray:
            """Calculate Re using the <w*T> formula.

            Args:
                dset (np.ndarray): The w*T quantity from analysis file.
                time_array (np.ndarray): The time array for calculating time
                averages.

            Returns:
                np.ndarray: The Reynolds number throughout time.
            """

            return np.sqrt( (Ra*dset) / (time_array * Pr))


        # Makes the for loop easier. First entry is the quantity, second is the
        # function to calculate Nu with, third is the LaTeX command for
        # displaying it in the plot.
        qoi_func = [
            (
                "avg_wT",
                fnu_wT,
                r"$\langle u_3 \cdot T \rangle$"),
                fre_wt
        ]


        n = len(qoi_func)

        """ Other unused functions
        ('avg_wtheta', fnu_wtheta, r'$\langle u_3\theta \rangle$'),
        ('avg_grad_theta_sq', fnu_grad_theta_sq, r'$\langle \left| \nabla \theta \right|^2 \rangle$')
        """

        # Plotting Nu - all quantities on 1 figure.
        # 2 columns, the first for instantaneous and the second for cumulative.
        # n rows, where n is the number of quantities we are calculating.
        if n > 1:
            fig, axes = plt.subplots(
                nrows=n,
                ncols=2,
                figsize=(20 * 2, 5 * n),
                layout="constrained",
            )
        else:
            fig, axes = plt.subplots(ncols=2, figsize=(20 * 2, 5), layout="constrained")

        # Indices for splitting into sections for convergence check. Splits the
        # time array into n_secs portions, with each portion covering the same
        # amount of time. These sections don't need to be the same size, if
        # timesteps were variable.
        n_secs = 3
        inds = np.searchsorted(
            time, [t_0 + i * total_time / n_secs for i in range(n_secs)]
        )
        inds = np.append(inds, [-1])

        # Array to hold the different Nu values (varies over section and qoi)
        nus = np.ones((n, n_secs+1))
        res = np.ones((n, n_secs+1))

        # Get time and volume integrals of quantities of interest (QoI) for the
        # correct time range
        for ind, (qoi, nu_func, lab, re_func) in enumerate(qoi_func):
            # The raw analysis file values to obtain <wT>, <|∇u|^2>, etc
            dset_nu = np.array(file["tasks"][qoi][:])[start_ind:, 0, 0]
            dset_re = np.array(file["tasks"]['avg_K'][:])[start_ind:, 0, 0]

            # Instantaneous Nusselt number at each time
            inst_dset_nu = np.diff(dset_nu)
            inst_dset_re = np.diff(dset_re)
            inst_nu = nu_func(inst_dset_nu, dt_array)
            inst_re = re_func(inst_dset_re, dt_array)

            # Cumulative Nusselt number at each point in time and for each
            # section of the sim, taking only the times past start_ave.
            # Currently dividing into thirds (n_secs=3), but could do quarters
            # if you're looking to be even more certain that the convergence
            # isn't a fluke.
            cumu_nus = np.append(1, nu_func(dset_nu[1:] - dset_nu[0], time[1:] - t_0))
            cumu_res = re_func(dset_re[1:] - dset_re[0], time[1:], t_0)

            # Final calculation of Nu
            nus[ind, -1] = cumu_nus[-1]
            res[ind, -1] = cumu_res[-1]

            # Write the data to the info file
            info.write(f"The following data are for the QoI {qoi}:\n")
            z = zip(inds[:-1], inds[1:])

            for i, (ind_1, ind_2) in enumerate(z):

                nus[ind, i] = nu_func(dset[ind_2]-dset[ind_1], time[ind_2]-time[ind_1])
                info.write(
                    f"Nu calculated using data from section {i+1} of {n_secs}:"
                    f" {nus[ind, i]:.4f}\n"
                )

            info.write(
                f"Maximum percent difference in Nusselt over {n_secs} sections:"
                f" {max_difference(nus[ind, :-1]):.4f}%\n"
            )
            info.write(
                "Nu as calculated by cumulative average: " "{:.6f}\n".format(nus[ind, -1])
            )
            info.write(f"Re as calculated by cumulative average: {res[ind, -1]}")
            info.write("-" * 80)
            info.write("\n")

            # Plot the instantaneous Nusselt number
            if n > 1:
                axes_ind = axes[ind]
            else:
                axes_ind = axes

            if len(inst_nu) > 40000:
                axes_ind[0].plot(time[1::10], inst_nu[::10])
            else:
                axes_ind[0].plot(time[1:], inst_nu, linewidth=1)
            axes_ind[0].plot([t_0, t_f], inst_ave_Nu * np.ones(2))
            axes_ind[0].set_xlim([t_0, t_f])
            axes_ind[0].set_title(r"Instantaneous Nu$(t)$ calculated via " + lab)
            axes_ind[0].set_xlabel("t")
            axes_ind[0].set_ylabel("Nu")

            # Plot the cumulative Nusselt number
            if len(cumu_nus) > 40000:
                axes_ind[1].plot(time[::10], cumu_nus[::10], linewidth=1)
            else:
                axes_ind[1].plot(time, cumu_nus, linewidth=1)
            axes_ind[1].plot([t_0, t_f], nus[ind, -1] * np.ones(2))
            axes_ind[1].set_xlim([t_0, t_f])
            axes_ind[1].set_title(r"Cumulative Nu$(t)$ calculated via " + lab)
            axes_ind[1].set_xlabel(r"$t$")
            axes_ind[1].set_ylabel("Nu")

        info.close()

        savename = output.joinpath("Nu_zoomed.png")
        fig.savefig(str(savename), dpi=400)
        plt.close()

        ########################################################################
        # Get data for profile plots
        ########################################################################

        # Get z-axis points
        z = np.array(file["scales/z/1.0"][:])

        # Horizontal temperature profiles, from start_time to end
        dset_temp = np.array(file["tasks"]["avg_T"][:])[start_ind:, 0, :]
        # Cumulative time average of horizontal profile; add background temp
        temp_prof = (dset_temp[-1, :] - dset_temp[0, :]) / (t_f - t_0)
        temp_prof = temp_prof + 1 / 2 - z

        # Horizontal kinetic energy and velocity profiles, from start time to
        # end
        dset_K = np.array(file["tasks"]["h_avg_K"][:])[start_ind:, 0, :]
        dset_u_sq = np.array(file["tasks"]["avg_u_sq"][:])[start_ind:, 0, :]
        dset_w_sq = np.array(file["tasks"]["avg_w_sq"][:])[start_ind:, 0, :]
        # Cumulative time average for profile
        K_prof = (dset_K[-1, :] - dset_K[0, :]) / (t_f - t_0)
        u_prof = (dset_u_sq[-1, :] - dset_u_sq[0, :]) / (t_f - t_0)
        w_prof = (dset_w_sq[-1, :] - dset_w_sq[0, :]) / (t_f - t_0)

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

    from pathlib import Path
    from dedalus.tools import logging
    from dedalus.tools import post
    from dedalus.tools.parallel import Sync
    import re

    base_fp = Path("/project/def-goluskin/RB_data/2D/")


    fp = list(base_fp.rglob("*/*_pr[13]*/*_analysis/analysis.h5"))
    start_times = {}

    for p in fp:
        info_fp = p.parents[1].rglob(f"*/info.txt").__next__()
        with open(info_fp, 'r') as file:
            for line in file:
                match1 = re.search(r'Time averaging.*?begins at:\s*([+-]?\d+(\.\d+)?)', line)
                if match1:
                    start_times[str(p)] = np.float64(match1.group(1))
                    break



    post.visit_writes(
        fp, main, start_ave=start_times,
    )
