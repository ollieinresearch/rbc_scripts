import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


def main(basepath: Path):

    output_dir = basepath.parents[1] / "rbc_simulations_output"
    output_dir.mkdir(exist_ok=True)
    
    
    for dim in range(2,4):
        fp = basepath / f"{dim}d/rapr.csv"

        data = pd.read_csv(fp)
        nrows, ncols = data.shape

        ra0 = 1e5
        pr0 = 1e0

        ra_values = ra0 * 10**np.arange(ncols)
        pr_values = pr0 * 10**(-0.5 * np.arange(nrows))

    
        def get_diagonal_slice(pr_0, alpha):
            nu_vals = []
            ra_vals = []
            pr_vals = []

            i = 0  # column index (r)
            j = pr_0

            while i < ncols and j < nrows:
                nu_vals.append(data[j, i])
                ra_vals.append(ra_values[i])
                pr_vals.append(pr_values[j])

                i += 1
                # Move vertically: Δj = log10(p_step) / log10(p_ratio)
                j += alpha

            return np.array(ra_vals), np.array(nu_vals)

        # === Plotting ===
        fig, ax = plt.subplots()

        alphas = [1,2,3,4]  # You can adjust this list

        for a in alphas:
            slices = []
            for pr in range(nrows-3):
                if pr + 2*a <= nrows:
                    slices.append(get_diagonal_slice(pr, a))
            for r, n in slices:
                plt.loglog(r, n, label=f"a = {a}")

        ax.set_xscale('log')
        ax.set_xlabel("r")
        ax.set_ylabel("Data value")
        ax.set_title("Diagonal slices through CSV")
        ax.legend()
        plt.show()



if __name__ == "__main__":

    main(Path(__file__))