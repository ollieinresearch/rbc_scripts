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

    analysis = basepath / "analysis"
    output = basepath / 'outputs'
    output.mkdir(exist_ok='true')


    n_secs = 3
    fs = analysis.glob('*.h5')
    fp = next(fs)

    with h5.File(fp, 'r') as f:
        Ra = float(f['tasks']['Ra'][-1])
        Pr = float(f['tasks']['Pr'][-1])

        full_time = f['scales']['sim_time'][:]
        avg_K = f['tasks']['avg_K'][:]
        avg_wT = f['tasks']['avg_wT'][:]
        avg_vorticity_sq = f['tasks']['avg_vorticity_sq'][:]
        avg_grad_T_sq = f['tasks']['avg_grad_T_sq'][:]
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

    for fi in fs:
        with h5.File(fi, 'r') as f:
            start_time = full_time[-1]
            if f['scales']['sim_time'][-1] > start_time:
                start_ind = np.searchsorted(f['scales']['sim_time'][:], start_time)
                full_time = np.append(full_time, f['scales']['sim_time'][start_ind:], axis=0)
                avg_K = np.append(avg_K, f['tasks']['avg_K'][start_ind:], axis=0)
                avg_wT = np.append(avg_wT, f['tasks']['avg_wT'][start_ind:], axis=0)
                avg_vorticity_sq = np.append(avg_vorticity_sq, f['tasks']['avg_vorticity_sq'][start_ind:], axis=0)
                avg_grad_T_sq = np.append(avg_grad_T_sq, f['tasks']['avg_grad_T_sq'][start_ind:], axis=0)
                avg_T = np.append(avg_T, f['tasks']['avg_T'][start_ind:], axis=0)
                avg_u_sq = np.append(avg_u_sq, f['tasks']['avg_u_sq'][start_ind:], axis=0)
                avg_w_sq = np.append(avg_w_sq, f['tasks']['avg_w_sq'][start_ind:], axis=0)

                # Currently only 3d, so no need to check if this will bei n the file or not
                avg_v_sq = np.append(avg_v_sq, f['tasks']['avg_v_sq'][start_ind:], axis=0)

    

    avg_K = np.ravel(avg_K)
    avg_wT = np.ravel(avg_wT)
    avg_vorticity_sq = np.ravel(avg_vorticity_sq)
    avg_grad_T_sq = np.ravel(avg_grad_T_sq)


    increased = full_time[1:] >= full_time[:-1]
    if not np.all(increased):
        print("not all increased")
        indxs = np.argsort(full_time)

        full_time = full_time[indxs] 
        avg_K = avg_K[indxs] 
        avg_wT = avg_wT[indxs] 
        avg_vorticity_sq = avg_vorticity_sq[indxs] 
        avg_grad_T_sq = avg_grad_T_sq[indxs] 
        avg_T = avg_T[indxs] 
        avg_u_sq = avg_u_sq[indxs] 
        avg_w_sq = avg_w_sq[indxs] 
        avg_v_sq = avg_v_sq[indxs] 
    
    
     
    start_ind = np.searchsorted(full_time, start_ave)
    time = full_time[start_ind:]
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

    
    def fnu_vorticity_sq(dset: np.ndarray, time_array) -> np.ndarray:
        """Calculate Nu using the <ω^2> formula."""

        return 1 + Pr * dset / time_array



    def fnu_grad_u_sq(dset: np.ndarray, time_array) -> np.ndarray:
        """Calculate Nu using the <|∇u|^2> formula."""

        return 1 + Pr * dset / time_array



    def fre_K(dset: np.ndarray, time_array) -> np.ndarray:
        """Calculate Re using the kinetic energy formula.

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
            avg_wT,
            "avg_wT",
            fnu_wT,
            r"$\langle u_3 \cdot T \rangle$",
        ),
        (
            avg_grad_T_sq,
            "avg_grad_T_sq",
            fnu_grad_T_sq,
            r"$\langle \left| \nabla T \right|^2 \rangle$",
        ),
        (
            avg_vorticity_sq,
            "avg_vorticity_sq",
            fnu_vorticity_sq,
            r"$\langle \left| \omega \right|^2 \rangle$",
        ),
    ]

    n = len(qoi_func)

    # Plotting Nu - all quantities on 1 figure.
    # 2 columns, the first for instantaneous and the second for cumulative.
    # n+1 rows, where n is the number of quantities we're calculating.
    # n rows are for Nu, 1 is for KE/Re
    if n > 1:
        fig, axes = plt.subplots(
            nrows=n+1,
            ncols=2,
            figsize=(12 * 2, 3 * (n+1)),
            layout="constrained",
        )

    else:
        fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(20 * 2, 5 * 2), layout="constrained")
    
    ke_ax = axes[0]
    
    inst_ke_ylab = r'$\frac{1}{\Gamma}\int_\Omega u^2+w^2 dxdz$'
    re_ylab = r'$\sqrt{\frac{Ra}{Pr}}\frac{1}{\Gamma}{\int_t^{t+100dt}\int_\Omega u^2+w^2 dxdzd\hat{t}}^{1/2}$'

    if dim == 3:
        inst_ke_ylab = r'$\frac{1}{\Gamma}\int_\Omega u^2+v^2+w^2 dxdydz$'
        re_ylab = r'$\sqrt{\frac{Ra}{Pr}}\frac{1}{\Gamma}{\int_t^{t+100dt}\int_\Omega u^2+v^2+w^2 dxdydzd\hat{t}}^{1/2}$'
    
    dset_re = avg_K[start_ind:]

    cumu_dset_re = cumulative_trapezoid(dset_re, time)
    cumu_re = fre_K(cumu_dset_re, time[1:] - t_0)
    final_re = fre_K(simpson(dset_re, time), total_time)

    # instantaneous KE
    ke_ax[0].plot(full_time, avg_K, linewidth=1)
    ke_ax[0].hlines(np.mean(avg_K), full_time[0], t_f, color='orange')   
    ke_ax[0].set_xlim([full_time[0], t_f])
    ke_ax[0].set_title('Instantaneous KE')
    ke_ax[0].set_xlabel(r'$t$')
    ke_ax[0].set_ylabel(inst_ke_ylab)

    # Cumulative Reynolds
    ke_ax[1].plot(time[1:], cumu_re, linewidth=1)
    ke_ax[1].hlines(final_re, full_time[0], t_f, color='orange')   
    ke_ax[1].set_xlim([t_0, t_f])
    ke_ax[1].set_title('Re')
    ke_ax[1].set_xlabel(r'$t$')
    ke_ax[1].set_ylabel(re_ylab)


    # Indices for splitting into sections for convergence check. Splits the
    # time array into n_secs portions, with each portion covering the same
    # amount of time. These sections don't need to be the same size, if
    # timesteps were variable.
    sec_times = [t_0 + i * total_time / n_secs for i in range(n_secs)]
    inds = np.searchsorted(
        time, sec_times
    )
    # We also want the endpoint as an index
    sec_times = np.append(sec_times, time[-1])
    inds = np.append(inds, [-1])

    # Array to hold the different Nu values (varies by section and quantity)
    nus = np.ones((n, n_secs+1))

    # Get time and volume integrals of quantities of interest (QoI) for the
    # correct time range
    for ind, (qoi, qoi_name, nu_func, lab) in enumerate(qoi_func):
        # The raw analysis file values to obtain <wT>, <|∇u|^2>, etc
        dset = qoi[start_ind:]

        # Instantaneous Nusselt number at each time
        inst_nu = nu_func(dset, 1.0)

        # Cumulative Nusselt number at each point in time and for each
        # section of the sim, taking only the times past start_ave.
        # Currently dividing into thirds (n_secs=3), but could do quarters
        # if you're looking to be even more certain that the convergence
        # isn't a fluke.

        """
        # Integrate over time axis and save the integral value at each time
        mask = ~np.isnan(dset)
        dset = dset[mask]
        time = time[mask]
        """

        cumu_dset = cumulative_trapezoid(dset, time, axis=0)
        cumu_nus = nu_func(cumu_dset, time[1:] - t_0)

        # Final cumulative Nu is the longest time average (so hopefully 
        # best converged)
        nus[ind, -1] = nu_func(simpson(dset, time, axis=0), total_time)

        # Write the Nu data to the info file
        info.write(f"The following data are for the QoI {qoi_name}:\n")
        
        # Need adjacent indices to figure out where the sections are
        z = zip(inds[:-1], inds[1:])
        for i, (ind_1, ind_2) in enumerate(z):
            # Save the Nu for the specific section
            nus[ind, i] = nu_func(simpson(dset[ind_1:ind_2], time[ind_1:ind_2], axis=0), time[ind_2]-time[ind_1])
            info.write(
                f"Nu calculated using data from section {i+1} ({time[ind_1]:.1f} to {time[ind_2]:.1f}):"
                f" {nus[ind, i]:.4f}\n"
            )

        info.write(
            f"Maximum percent difference in Nusselt-1 over {n_secs} sections:"
            f" {max_difference(nus[ind, :-1] - 1):.4f}%\n"
        )
        info.write(
            f"Nu as calculated by cumulative average:"
            f" {nus[ind, -1]:.6f}\n"
        )
        info.write("-" * 80)
        info.write("\n")

        axes_ind = axes[ind+1]

        # Plot the instantaneous Nusselt number
        max_points = 40000
        skip = 10
        if len(inst_nu) > max_points:
            axes_ind[0].plot(time[::skip], inst_nu[::skip])
        else:
            axes_ind[0].plot(time, inst_nu, linewidth=1)
        # Line to show final average
        for i, (ind_1, ind_2) in enumerate(zip(sec_times[:-1], sec_times[1:])):
            # Display the Nu for the specific section
            axes_ind[0].hlines(nus[ind, i], ind_1, ind_2)
        axes_ind[0].hlines(nus[ind,-1], t_0, t_f, color='orange')
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
    fig.savefig(output / "info.pdf", dpi=400)

    # Save the info file
    info.write(
        f"Maximum percent difference in Nusselt-1 over all QoI:"
        f" {max_difference(nus[:, -1]-1):.4f}%\n"
    )
    info.write(
        f"Final Nusselt number:"
        f" {nus[0, -1]:.6f}\n"
    )
    info.write(
        f"Final Reynolds number:"
        f" {final_re:.6f}\n"
    )
    info.close()


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
        shps = [x.shape for x in dsets]

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

    shps = [x.shape for x in profs]
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

################################################################################

if __name__ == "__main__":
    from docopt import docopt

    # Collect arguments
    args = docopt(__doc__)

    basepath = Path(args["<basepath>"][0])
    time = float(args["--time"])
    end_time = float(args["--end_time"])

    main(basepath, time, end_time)
