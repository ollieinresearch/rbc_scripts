# Rayleigh–Bénard Convection Simulation Scripts

Direct numerical simulations of three-dimensional Rayleigh–Bénard convection at small
Prandtl number with [Dedalus](https://dedalus-project.readthedocs.io), run on the clusters
of the Digital Research Alliance of Canada.

This repository holds the simulation, the job scripts for Fir and Trillium, the `rbc`
package (restarts, checkpoint management, analysis, convergence and resolution checks,
results collection) and every figure of the paper.

**Documentation: https://ollieinresearch.github.io/rbc_scripts/** (built from `docs/`; to build it
locally, `pip install -e ".[docs]" && sphinx-build -b html docs docs/_build/html`).

## Install on a cluster

```bash
cd $SCRATCH && git clone https://github.com/ollieinresearch/rbc_scripts.git
cd rbc_scripts && bash shell/install_envs.sh     # creates $SCRATCH/ded3 and $SCRATCH/ded3_vis
```

## Run

```bash
mkdir -p $SCRATCH/ollie_rb_data/3d/ra1e5/pr1e0/48_Gam2_2 && cd $_
cp $SCRATCH/rbc_scripts/shell/{params.conf,run_sim.sh} .      # run_sim_trillium.sh on Trillium
# edit params.conf and the #SBATCH lines, then, from a login node:
bash $SCRATCH/rbc_scripts/shell/submit.sh
```

The documentation's *Getting started* section walks through this first run in full.

## Layout

```text
sim/rayleigh_benard.py   the Dedalus simulation
shell/                   job scripts and workflow functions
rbc/                     analysis, restarts, collection, figures, visualization
data/cases.csv           registry of parameter points (results.csv is kept out of git)
docs/                    the documentation
```

## License

GPL-3.0-or-later (see `LICENSE`), as for Dedalus.
