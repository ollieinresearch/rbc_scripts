import h5py  # pyright: ignore
import numpy as np
from scipy.integrate import cumulative_trapezoid, simpson
from pathlib import Path
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use("Agg")
plt.rcParams.update({"font.size": 14})
plt.ioff()

def main(file: Path, basepath: Path, start_ave: np.float64, end_ave: np.float64):
    """Computes helpful quantities and creates plots for a (hopefully)
    converged flow.

    Args:
        filename (str | PathLike): The file path to the analysis file.
        basepath (str | PathLike): The basepath directory.
        start_ave (float): The time at which to begin the analysis of the flow.
        end_ave (float): Time at which to end analysis; 0 implies the end.
    """

    output = basepath / 'outputs'
    output.mkdir(exist_ok='true')


    n_secs = 3

    with h5py.File(file, mode="r") as f:
        # Rayleigh and Prandtl number from the file
        Ra = float(f["tasks"]["Ra"][-1])
        Pr = float(f["tasks"]["Pr"][-1])

        # Get time from file, then take data only from those times beyond where
        # averaging begins
        full_time = np.array(f["scales"]["sim_time"])
        start_ind = np.searchsorted(full_time, start_ave)
        time = full_time[start_ind:]
        t_0 = time[0]
        t_f = time[-1]
        total_time = t_f - t_0

        avg_K = np.array(f["tasks"]["avg_K"])
        dim = len(avg_K.shape) - 1

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
            dsets = [np.array(f["tasks"][avg][start_ind:, 0, 0, :]) for avg in avgs]
            dsets[0] += 1/2-z
            # Kinetic energy is just the sum of the squared avg velocities
            dsets.append(dsets[1]+dsets[2]+dsets[3])

            horz_tex = r"$\sqrt{\overline{u^2+v^2}}$"
            kin_tex = r"$\overline{u^2+v^2+w^2}$"

        else:
            avgs = [
                "avg_T",
                "avg_w_sq",
                "avg_u_sq",
            ]
            dsets = [np.array(f["tasks"][avg][start_ind:, 0, :]) for avg in avgs]
            dsets[0] += 1/2-z
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

        grad_ke = np.gradient(profs[-1])
        quad_interp = np.poly1d(np.polynomial.polyfit(z, grad_ke, 2))

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

        ax.plot(quad_interp(z), z, 'r')        

        fig.tight_layout(pad=1.0)
        savename = output.joinpath("profiles.png")
        fig.savefig(str(savename), dpi=400)
        plt.close()



if __name__ == "__main__":

    fp = Path()
    basepath = fp.parent
    main(fp, basepath, start_time, 0)