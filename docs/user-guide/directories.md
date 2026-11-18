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
                      │   │     │     └── label: <N_z>[_<N_x>]_Gam2_2[_<stepper>][_safety<k>]
                      │   │     └──────── Pr = 1   (pr<mantissa>e<exponent>)
                      │   └────────────── Ra = 10^8
                      └────────────────── dimension
```

Ra and Pr are written as `1e<k>` on whole decades and `3e<k>` on half-decades, so
$Pr = 10^{-2.5} \approx 3.2\times10^{-3}$ is `pr3e-3`. {func}`rbc.cases.parse_run_path`
reads them back as exact half-decade exponents ($\log_{10}$ values), which is how every
script identifies a run; {func}`rbc.cases.folder_value` goes the other way.

The label starts with the number of Chebyshev modes $N_z$. `Gam2_2` is the
{term}`aspect ratio` Γ = (2,2). The baseline run has $N_x = N_y = 2N_z$, RK222 and the
default CFL safety factor; anything else is written into the label
({func}`rbc.cases.parse_label`):

| Label | Run |
|---|---|
| `192_Gam2_2` | $N_z = 192$, $N_x = N_y = 384$, RK222, default safety (the baseline) |
| `192_384_Gam2_2` | the same, with $N_x$ written out; any other $N_x$ is written this way, e.g. `192_512_Gam2_2` |
| `192_Gam2_2_rk443` | timestepper RK443 (any Dedalus scheme name, lower case: `rk111`, `sbdf2`, `cnab2`, ...) |
| `192_Gam2_2_safety2` | CFL safety factor divided by 2 (`safety4`: by 4; `safety1p5`: by 1.5) |
| `192_Gam2_2_rk443_safety2` | both |

$N_x$ must be even (RealFourier) and should factor into small primes for the FFTs.

The label is a name; `params.conf` is what the run used (`RES_Z`, `RES_HOR`, `STEPPER`,
`CFL_SAFETY`), so it wins where the two disagree, and the disagreement is reported
(`config_mismatch` in `results.csv`, a warning in `variants.txt`). The grid sizes used by
the analysis come from the output files themselves.

Runs in one parameter folder that differ in resolution, timestepper or safety factor are
compared automatically ({mod}`rbc.analysis.variants`, {ref}`resolution-error-budget`). Only
baseline runs (no stepper or safety suffix) are candidates for the {term}`production run`.

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
