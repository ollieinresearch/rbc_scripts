"""
Script to perform analysis tasks on data obtained from running a 3D
Rayleigh-Benard convection simulation. Accepts an analysis file as a command
line argument and takes time averages from the specified time. All output is
written to the folders 'outputs' and 'preliminary_outputs'. Outputs plots of 
the Nusselt number, various profiles and information about the simulation to a 
text document.

Usage:
    analysis.py <files>... [--time=<time>] [--basepath=<dir>]
    analysis.py <files>...

Options:
    --time=<time>   Time at which to begin time averaging [default: 200]
    --basepath=<dir>  Path to parent folder for output [default: ./analysis]
"""

import h5py  # pyright: ignore
import numpy as np
from scipy.integrate import cumulative_trapezoid, simpson
from pathlib import Path
import matplotlib
import matplotlib.pyplot as plt

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



def main(file: Path, basepath: Path, start_ave: np.float64):
    """Computes helpful quantities and creates plots for a (hopefully)
    converged flow.

    Args:
        filename (str | PathLike): The file path to the analysis file.
        basepath (str | PathLike): The basepath directory.
        start_ave (float): The time at which to begin the analysis of the flow.
    """

    output = basepath / 'outputs'
    output.mkdir(exist_ok='true')
    pre_output = basepath / "preliminary_outputs"
    pre_output.mkdir(exist_ok='true')

    n_secs = 3

    with h5py.File(file, mode="r") as f:
        # Rayleigh and Prandtl number from the file
        Ra = float(f["tasks"]["Ra"][-1])
        Pr = float(f["tasks"]["Pr"][-1])

        # Get time from file, then take data only from those times beyond where
        # averaging begins
        full_time = np.array(f["scales/sim_time"])
        start_ind = np.searchsorted(full_time, start_ave)
        time = full_time[start_ind:]
        t_0 = time[0]
        t_f = time[-1]
        total_time = t_f - t_0

        # Preliminary output
        # write basics to text file
        savename = pre_output / 'prelim.txt'
        prelim = open(str(savename), 'w')
        prelim.write(f'Ra= {Ra:.4e}, Pr= {Pr:.4e}\n')
        prelim.write(f'End sim time: {t_f:.3f}\n')
        prelim.write(f'Start sim time: {full_time[0]:.3f}')
        prelim.close()

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

        # Plot full instantaneous KE
        avg_K = np.array(f["tasks"]["avg_K"])
        dim = len(avg_K.shape) - 1


        ########################################################################
        # Calculating Nu with different quantities        
        ########################################################################
        
        # Different functions for different datasets
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


        def fnu_grad_T_sq(dset: np.ndarray, time_array) -> np.ndarray:
            """Calculate Nu using the <|∇T|^2> formula."""

            return dset / time_array

        
        def fnu_oy_sq(dset: np.ndarray, time_array) -> np.ndarray:
            """Calculate Nu using the <oy^2> formula."""

            return 1 + Pr * dset / time_array



        def fnu_grad_u_sq(dset: np.ndarray, time_array) -> np.ndarray:
            """Calculate Nu using the <|∇u|^2> formula."""

            return 1 + Pr * dset / time_array



        # Makes the for loop easier. First entry is the quantity, second is the
        # function to calculate Nu with, third is the LaTeX command for
        # displaying it in the plot.
        qoi_func = [
            (
                "avg_wT",
                fnu_wT,
                r"$\langle u_3 \cdot T \rangle$",
            ),
            (
                "avg_grad_T_sq",
                fnu_grad_T_sq,
                r"$\langle \left| \nabla T \right|^2 \rangle$",
            ),
        ]

        if dim == 3:
            qoi_func.append((
                "avg_grad_u_sq",
                fnu_grad_u_sq,
                r"$\langle \left| \nabla u \right|^2 \rangle$",
            ))
        else:
            qoi_func.append((
                "avg_oy_sq",
                fnu_oy_sq,
                r"$\langle \left| \omega \right|^2 \rangle$",
            ))


        n = len(qoi_func)

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
            pre_fig, pre_axes = plt.subplots(
                nrows=n+2,
                figsize=(20, 5 * (n + 1)),
                layout="constrained",
            )

        else:
            fig, axes = plt.subplots(ncols=2, figsize=(20 * 2, 5), layout="constrained")

        
        inst_ke_ylab = r'$\frac{1}{\Gamma}\int_\Omega u^2+w^2 dxdz$'
        time_avg_ylab = r'$\frac{1}{100dt}\frac{1}{\Gamma}\int_t^{t+100dt}\int_\Omega u^2+w^2 dxdzd\hat{t}$'

        if dim == 3:
            inst_ke_ylab = r'$\frac{1}{\Gamma}\int_\Omega u^2+v^2+w^2 dxdydz$'
            time_avg_ylab = r'$\frac{1}{100dt}\frac{1}{\Gamma}\int_t^{t+100dt}\int_\Omega u^2+v^2+w^2 dxdydzd\hat{t}$'

        inst_K = np.ravel(avg_K)
        pre_axes[0].plot(full_time, inst_K, linewidth=1)
        pre_axes[0].hlines(np.mean(inst_K), full_time[0], t_f, color='orange')   
        pre_axes[0].set_xlim([full_time[0], t_f])
        pre_axes[0].set_title('Instantaneous KE')
        pre_axes[0].set_xlabel(r'$t$')
        pre_axes[0].set_ylabel(inst_ke_ylab)

        inst_avg_K = cumulative_trapezoid(inst_K, full_time) / (full_time[1:] - full_time[0])
        pre_axes[1].plot(full_time[1:], inst_avg_K, linewidth=1)
        pre_axes[1].hlines(inst_avg_K[-1], full_time[0], t_f, color='orange')   
        pre_axes[1].set_xlim([full_time[0], t_f])
        pre_axes[1].set_title('Instantaneous time average of KE')
        pre_axes[1].set_xlabel(r'$t$')
        pre_axes[1].set_ylabel(time_avg_ylab)


        # Indices for splitting into sections for convergence check. Splits the
        # time array into n_secs portions, with each portion covering the same
        # amount of time. These sections don't need to be the same size, if
        # timesteps were variable.
        
        inds = np.searchsorted(
            time, [t_0 + i * total_time / n_secs for i in range(n_secs)]
        )
        # We also want the endpoint as an index
        inds = np.append(inds, [-1])

        # Array to hold the different Nu values (varies by section and quantity)
        nus = np.ones((n, n_secs+1))

        # Get time and volume integrals of quantities of interest (QoI) for the
        # correct time range
        for ind, (qoi, nu_func, lab) in enumerate(qoi_func):
            # The raw analysis file values to obtain <wT>, <|∇u|^2>, etc
            dset_full = np.ravel(np.array(f["tasks"][qoi]))
            dset = dset_full[start_ind:]

            # Instantaneous Nusselt number at each time
            inst_nu_full = nu_func(dset_full, 1.0)
            inst_ave_Nu_full = np.mean(inst_nu_full)
            inst_nu = inst_nu_full[start_ind:]
            inst_ave_Nu = np.mean(inst_nu)

            ####################################################################
            # Plot full plots of instantaneous Nu
            pre_axes[ind+2].plot(full_time, inst_nu_full, linewidth=1)
            pre_axes[ind+2].hlines(inst_ave_Nu_full, full_time[0], t_f, color='orange')
            pre_axes[ind+2].set_xlim([full_time[0], t_f])
            #pre_axes[ind+2].set_ylim([np.floor(np.min(inst_nu_full)), np.ceil(np.max(inst_nu_full))])
            pre_axes[ind+2].set_title(r'Instantaneous Nu($t$) via ' + lab)
            pre_axes[ind+2].set_xlabel(r'$t$')
            pre_axes[ind+2].set_ylabel('Nu')

            # Cumulative Nusselt number at each point in time and for each
            # section of the sim, taking only the times past start_ave.
            # Currently dividing into thirds (n_secs=3), but could do quarters
            # if you're looking to be even more certain that the convergence
            # isn't a fluke.

            # Integrate over time axis and save the integral value at each time
            cumu_dset = cumulative_trapezoid(dset, time)
            cumu_nus = nu_func(cumu_dset, time[1:] - t_0)

            # Final cumulative Nu is the longest time average (so hopefully 
            # best converged)
            nus[ind, -1] = nu_func(simpson(dset, time, axis=0), total_time)

            # Write the Nu data to the info file
            info.write(f"The following data are for the QoI {qoi}:\n")
            
            # Need adjacent indices to figure out where the sections are
            z = zip(inds[:-1], inds[1:])
            for i, (ind_1, ind_2) in enumerate(z):
                # Save the Nu for the specific section
                nus[ind, i] = nu_func(cumu_dset[ind_2]-cumu_dset[ind_1], time[ind_2]-time[ind_1])
                info.write(
                    f"Nu calculated using data from section {i+1} of {n_secs}:"
                    f" {nus[ind, i]:.4f}\n"
                )

            info.write(
                f"Maximum percent difference in Nusselt over {n_secs} sections:"
                f" {max_difference(nus[ind, :-1] - 1):.4f}%\n"
            )
            info.write(
                f"Nu as calculated by average of instantaneous Nu:"
                f" {inst_ave_Nu:.6f}\n"
            )
            info.write(
                f"Nu as calculated by cumulative average:"
                f" {nus[ind, -1]:.6f}\n"
            )
            info.write("-" * 80)
            info.write("\n")

            # Check which index should be used (one QoI vs many)
            if n > 1:
                axes_ind = axes[ind]
            else:
                axes_ind = axes

            # Plot the instantaneous Nusselt number
            max_points = 40000
            skip = 10
            if len(inst_nu) > max_points:
                axes_ind[0].plot(time[::skip], inst_nu[::skip])
            else:
                axes_ind[0].plot(time, inst_nu, linewidth=1)
            # Line to show average instantaneous Nu
            axes_ind[0].hlines(inst_ave_Nu, t_0, t_f, color='orange')
            axes_ind[0].set_xlim([t_0, t_f])
            axes_ind[0].set_title(r"Instantaneous Nu$(t)$ calculated via " + lab)
            axes_ind[0].set_xlabel(r"$t$")
            axes_ind[0].set_ylabel("Nu")

            # Plot the cumulative Nusselt number
            if len(cumu_nus) > max_points:
                axes_ind[1].plot(time[1::skip], cumu_nus[::skip], linewidth=1)
            else:
                axes_ind[1].plot(time[1:], cumu_nus, linewidth=1)
            # Line to show the final average
            axes_ind[1].hlines(nus[ind, -1], t_0, t_f, color='orange')
            axes_ind[1].set_xlim([t_0, t_f])
            axes_ind[1].set_title(r"Cumulative Nu$(t)$ calculated via " + lab)
            axes_ind[1].set_xlabel(r"$t$")
            axes_ind[1].set_ylabel("Nu")

        # Save the cumu/inst Nu plot, and full time plot
        fig.savefig(output / "Nu.png", dpi=400)
        pre_fig.savefig(pre_output / 'prelim_time_averages.png', dpi=400)

        # Save the info file
        info.write(
            f"Maximum percent difference in Nusselt over all QoI:"
            f" {max_difference(nus[:, -1]-1):.4f}%\n"
        )
        info.close()


        ########################################################################
        # Profile plots
        ########################################################################

        # Get z-axis points
        z = np.array(f["scales/z/1.0"][:])

        # Horizontal temperature profiles, from start_time to end
        if dim == 3:
            avgs = [
                "avg_T",
                "avg_w_sq",
                "avg_u_sq",
                "avg_v_sq",
            ]
            dsets = np.array([np.array(f["tasks"][avg][start_ind:, 0, 0, :]) for avg in avgs])
            print(dsets.shape)
            # Kinetic energy is just the sum of the squared avg velocities
            dsets = np.append([dsets, dsets[1]+dsets[2]+dsets[3]], axis=0)

            horz_tex = r"$\sqrt{\overline{u^2+v^2}}$"
            kin_tex = r"$\overline{u^2+v^2+w^2}$"

        else:
            avgs = [
                "avg_T",
                "avg_w_sq",
                "avg_u_sq",
            ]
            dsets = np.array([np.array(f["tasks"][avg][start_ind:, 0, :]) for avg in avgs])
            print(dsets.shape)

            # Kinetic energy is just the sum of the squared avg velocities
            dsets = np.append([dsets, dsets[1]+dsets[2]], axis=0)


            horz_tex = r"$\sqrt{\overline{u^2}}$"
            kin_tex = r"$\overline{u^2+w^2}$"

        print(f"dsets shape = {dsets.shape}")
        profs = [simpson(dset, time, axis=0) / total_time for dset in dsets]
        print(f"cumu_shape = {cumu_dsets.shape}")

        plot_ops = [
            (
                r"$\langle u_3 \cdot T \rangle$",
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
            if title == "Horizontally Averaged Temperature Profile":
                plt.xlim([0,1])
            plt.ylim([-0.5, 0.5])
            plt.xlabel(xlabel)
            plt.ylabel(r"$z$")
            plt.title(title)        

        fig.tight_layout(pad=1.0)
        savename = output.joinpath("profiles.png")
        fig.savefig(str(savename), dpi=400)
        plt.close()

################################################################################

if __name__ == "__main__":
    from docopt import docopt

    # Collect arguments
    args = docopt(__doc__)

    file = Path(args["<files>"][0])
    basepath = Path(args["--basepath"])
    time = float(args["--time"])

    main(file, basepath, time)
