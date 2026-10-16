---
html_theme.sidebar_secondary.remove: true
---

# Rayleigh–Bénard Convection Simulation Scripts

Direct numerical simulations of three-dimensional {term}`Rayleigh–Bénard convection`
at small {term}`Prandtl numbers <Prandtl number>`, run with [Dedalus](https://dedalus-project.readthedocs.io)
on the clusters of the [Digital Research Alliance of Canada](https://alliancecan.ca), and
everything around them: job scripts, restarts, checkpoint management, analysis,
resolution and convergence checks, and the figures of the paper.

The project measures how the heat transport, the {term}`Nusselt number`, depends on the
{term}`Rayleigh number` and the Prandtl number, in particular along {term}`slices <slice>`
$Pr \propto Ra^{-1/2}$ of parameter space, where the theory of Shishkina and Lohse
makes a definite prediction ([Shishkina & Lohse 2024](https://arxiv.org/abs/2407.16573); see {doc}`about/references`). {doc}`about/history` explains how the project started and
why the main choices were made.

::::{grid} 1 2 2 2
:gutter: 3

:::{grid-item-card} {fas}`rocket` Getting started
:link: getting-started/index
:link-type: doc

New to the project? Accounts, installing the environments on a cluster, and a complete
first run from submission to figures.
:::

:::{grid-item-card} {fas}`book-open` User guide
:link: user-guide/index
:link-type: doc

How the workflow fits together: parameters, running and restarting, storage and
backups, post-processing, the results table, figures and visualizations.
:::

:::{grid-item-card} {fas}`square-root-variable` Theory and methods
:link: theory/index
:link-type: doc

The equations and units, the Nusselt-number identities, boundary layers, statistical
convergence, spectral resolution, scaling exponents, and the parallel-I/O details.
:::

:::{grid-item-card} {fas}`code` Reference
:link: reference/index
:link-type: doc

Every command and its options, every Python function, the output files and the columns
of `results.csv`, and the {doc}`glossary`.
:::

:::{grid-item-card} {fas}`images` Gallery
:link: gallery/index
:link-type: doc

Every paper figure, generated from synthetic data with the command that makes it, and
pictures and movies from the simulations.
:::

:::{grid-item-card} {fas}`life-ring` Troubleshooting
:link: troubleshooting
:link-type: doc

Symptoms, causes and fixes for the problems already met: hangs, restarts, quota,
timesteps, MPI binding and more.
:::
::::

## Five commands to know

```bash
bash $SCRIPTS_PATH/shell/submit.sh            # queue a run and its post-processing (login node)
python -m rbc.collect                         # gather every run into data/results.csv
python -m rbc.plots.all --png                 # make every paper figure
python -m rbc.analysis.nusselt RUN_DIR        # Nu, Re and convergence for one run
scancel --signal=USR1 <jobid>                 # stop a run cleanly, writing a checkpoint
```

Hover over underlined terms and links for a preview; click a {fas}`chevron-down`
**Derivation** box to expand the mathematics behind a statement.

```{toctree}
:hidden:
:maxdepth: 2

getting-started/index
user-guide/index
theory/index
reference/index
gallery/index
troubleshooting
development
about/index
glossary
```
