# Getting started

These pages take you from no account to a finished first simulation with its figures.
Work through them in order:

1. {doc}`accounts`: an Alliance account, multifactor authentication, SSH keys and access
   to the group's allocations.
2. {doc}`install`: the repository and the two Python environments on a cluster.
3. {doc}`first-run`: a small simulation ($Ra = 10^5$, $Pr = 1$, one node, one hour),
   followed from submission to Nusselt number, restart and figures.

After that, the {doc}`../user-guide/index` explains each part of the workflow in depth,
and {doc}`../theory/index` explains the physics and numerics behind it.

## What you are assumed to know

Python, and calculus up to partial differential equations. Nothing else: the pages explain
or link to what you need about the Linux shell, git, the Slurm scheduler, spectral methods
and convection. If something is unfamiliar, these are good starting points:

| Topic | Where to learn it |
|---|---|
| The Linux shell | [Alliance: Linux introduction](https://docs.alliancecan.ca/wiki/Linux_introduction) |
| Running jobs with Slurm | [Alliance: Running jobs](https://docs.alliancecan.ca/wiki/Running_jobs) |
| Python on the clusters | [Alliance: Python](https://docs.alliancecan.ca/wiki/Python) |
| Storage on the clusters | [Alliance: Storage and file management](https://docs.alliancecan.ca/wiki/Storage_and_file_management) |
| git | [Pro Git, chapters 1–3](https://git-scm.com/book/en/v2) |
| Dedalus | [Dedalus tutorials](https://dedalus-project.readthedocs.io/en/latest/pages/tutorials.html) |
| Rayleigh–Bénard convection | {doc}`../theory/equations`, then the references in {doc}`../about/references` |

```{tip}
Keep the {doc}`../glossary` open in another tab. Every term underlined in these pages
links to it, and hovering shows the definition.
```

```{toctree}
:hidden:

accounts
install
first-run
```
