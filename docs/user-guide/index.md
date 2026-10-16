# User guide

How the pieces of the workflow fit together, and how to use each one. The life of a run:

1. **Set up** a run directory with `params.conf` and a job script
   ({doc}`directories`, {doc}`parameters`).
2. **Submit** from a login node with `submit.sh`, which queues the simulation and a
   dependent post-processing job ({doc}`running`).
3. **Run.** The simulation starts fresh or from a checkpoint ({doc}`restarts`), writes
   output, and stops itself before the time limit with a final checkpoint
   ({ref}`graceful-stop`).
4. **Post-process.** One-node job: Nu, Re and convergence, profiles, spectra, optional
   3D renders ({doc}`post-processing`).
5. **Back up and prune.** Outputs and the newest checkpoints are copied to home and project
   space; old checkpoints are deleted from scratch ({doc}`storage`).
6. **Collect.** Every run becomes a row of `results.csv` ({doc}`cases-and-results`).
7. **Repeat** 2–6 until the run is converged and resolved
   ({doc}`../theory/convergence`, {doc}`../theory/resolution`).
8. **Figures** for the paper ({doc}`figures`) and visualizations ({doc}`visualization`).

Fir and Trillium differ in a few places, all collected in {doc}`clusters`.

```{toctree}
:maxdepth: 1

workflow
directories
parameters
running
restarts
storage
post-processing
cases-and-results
figures
visualization
clusters
```
