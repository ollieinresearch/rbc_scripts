# Directories and naming

## The repository

```text
rbc_scripts/
├── sim/rayleigh_benard.py     the Dedalus simulation
├── shell/                     job scripts and workflow functions (see the workflow page)
├── rbc/                       the Python package
│   ├── cases.py               naming conventions, slices, loading cases.csv and results.csv
│   ├── io.py                  reading file-handler output across restarts
│   ├── collect.py             all runs -> results.csv
│   ├── resources.py           job history and core-hours of a run
│   ├── restart.py             choosing, keeping and pruning checkpoints
│   ├── regrid.py              changing the resolution of a checkpoint
│   ├── shutdown.py            collective graceful stop
│   ├── vast.py                file-handler fixes for Trillium, fast restart loading
│   ├── style/                 Matplotlib style sheets, colours and markers
│   ├── analysis/              nusselt, profiles, exponents
│   ├── spectra/               core, cumulative, snapshots, compare
│   ├── plots/                 one module per paper figure, and all
│   └── vis/                   fields, slices, jointpdf, render3d
├── data/
│   ├── cases.csv              the registry of parameter points (in git)
│   └── results.csv            written by rbc.collect (not in git until publication)
├── docs/                      these pages
├── requirements.txt           simulation environment
└── requirements_vis.txt       visualization environment
```

The old scripts and where their code went are listed in {ref}`old-to-new`.

## The run tree

Every simulation is one {term}`run directory` on scratch:

```text
$SCRATCH/ollie_rb_data/3d/ra1e8/pr1e0/108_Gam2_2
                      │   │     │     └── label: <N_z>_Gam<Γ>_2
                      │   │     └──────── Pr = 1   (pr<mantissa>e<exponent>)
                      │   └────────────── Ra = 10^8
                      └────────────────── dimension
```

Ra and Pr are written as `1e<k>` on whole decades and `3e<k>` on half-decades, so
$Pr = 10^{-2.5} \approx 3.2\times10^{-3}$ is `pr3e-3`. {func}`rbc.cases.parse_run_path`
reads them back as exact half-decade exponents ($\log_{10}$ values), which is how every
script identifies a run; {func}`rbc.cases.folder_value` goes the other way.

The label starts with the number of Chebyshev modes $N_z$. `Gam2_2` is the
{term}`aspect ratio` Γ = (2,2) (so $N_x = 2N_z, N_y = 2N_z$).


The analysis takes the true grid sizes from the files, not from the label, so a label
that does not follow the pattern only affects sorting.

## Inside a run directory

```text
108_Gam2_2/
├── params.conf               settings of this run
├── run_sim.sh                job script (run_sim_trillium.sh on Trillium)
├── post.sh, submit.sh        post-processing job and submission (copied by sim_ready.sh)
├── jobs/                     job_sim_<jobid>.out and job_post_<jobid>.out, one log per job
├── stop_reason.txt           why the last job stopped
├── state/                    checkpoints: state_s<N>.h5 (+ state_s<N>/ process files)
├── restart/restart.h5        link to (or regridded copy of) the checkpoint to restart from
├── analysis/                 volume averages, every A_FREQ time units
├── horizontal_analysis/      horizontally averaged profiles, every A_FREQ
├── spectra/                  spectra computed in the run, every SPEC_FREQ
├── snapshots/                full 3D fields (only with SNAPSHOTS=1)
├── outputs/                  everything the analysis writes (below)
├── visualization/            3D renders and movies
└── watchdog/                 stack dumps, only if a job hung
```

`outputs/` holds `nusselt.json`, `info.txt`, `info.pdf`, `profiles.json`,
`profiles.npz`, `profiles.pdf`, `spectra.json`, `spec_range.txt`, the spectra figures and
`jobs.csv`; each is described in {doc}`../reference/files`.

Run directories made with the earlier scripts also contain `plotting.sh` and
`post_process.sh`, and their job logs at the top level; those scripts are no longer used.
`sim_ready.sh` is current again, as the helper that copies the scripts into a run
directory ({doc}`../getting-started/first-run`).

## Where copies live

The same tree is mirrored, path for path below `$SCRATCH`, into home and project space by
{doc}`storage`:

| Copy | Path of the example run | Contents |
|---|---|---|
| working | `/scratch/ollie/ollie_rb_data/3d/ra1e8/pr1e0/108_Gam2_2` | everything |
| home | `$HOME/ollie_rb_data/3d/ra1e8/pr1e0/108_Gam2_2` | outputs, analysis, spectra, logs, 2 newest checkpoints |
| project | `/project/def-goluskin/ollie/ollie_rb_data/3d/ra1e8/pr1e0/108_Gam2_2` | the same |
