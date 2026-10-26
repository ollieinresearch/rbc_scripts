# Post-processing a run

`post.sh` does this automatically after every simulation job; each step can also be run
by hand (with the simulation environment active).

## Nusselt and Reynolds numbers

```bash
python -m rbc.analysis.nusselt RUN --t0 auto
```

reads `analysis/` and writes `outputs/nusselt.json`, `outputs/info.txt` and
`outputs/info.pdf`:

- Nu from the three exact identities ({doc}`../theory/nusselt`) and from the plate
  fluxes;
- Re;
- for each, the statistics behind {doc}`../theory/convergence`: section means, the
  batch-means standard error (and the autocorrelation one for comparison), autocorrelation
  time, effective samples, and the section and trend tests;
- the plate balance, reported but not part of the convergence rule;
- `converged` (the original section rule) and `converged_stat` (the statistical rule).

`--t0` is the start of the averaging window: a time, or `auto` for the MSER-5 estimate of
the end of the initial transient. `info.pdf` shows the scan of the mean against the
window start, which is the check that the choice does not matter, and the running
autocorrelation time with the block means behind the standard error.

## Profiles and boundary layers

```bash
python -m rbc.analysis.profiles RUN --t0 T0
```

reads `horizontal_analysis/` and writes the time-averaged profiles (`profiles.npz`), the
boundary-layer thicknesses and node counts (`profiles.json`) and a figure
(`profiles.pdf`). The total heat flux $Nu(z)$ is plotted as a check: it is independent of
$z$ in a converged run, so its spread measures how far from stationary the window is
({doc}`../theory/boundary-layers`).

## Spectra

```bash
python -m rbc.spectra.cumulative RUN --time T0
python -m rbc.spectra.compare RUN/..
```

The first writes the time-averaged spectra of one run, the decades table
`spec_range.txt` and `spectra.json`; the second overlays every resolution of the same
(Ra, Pr) in `outputs/spectra_compare.pdf` of the parameter folder
({doc}`../theory/resolution`).

## Choosing the averaging start

`AVG_TIME` in `params.conf` sets it for the automatic post-processing and for
{mod}`rbc.collect`. `auto` is a sensible default; a number is reproducible and can be
chosen from `info.pdf`. Profiles and spectra always use the start that the Nusselt
analysis settled on, so the three agree.

(post-processing-cost)=
## Cost and history of a run

```bash
python -m rbc.resources RUN --jobs
```

prints the number of jobs, simulated time, iterations, mean timestep, wall-clock hours and
core-hours. Jobs since `outputs/jobs.csv` was introduced are exact. Older ones are
reconstructed from the output files, as explained in {mod}`rbc.resources`: the wall time
is a slight underestimate (start-up is missing) and the number of ranks is taken from the
current `params.conf`.

Slurm's accounting database may still know older jobs:

```bash
sacct -u $USER -S 2020-01-01 -X --format=JobID,JobName%20,Elapsed,NCPUS,CPUTime,State,WorkDir%90
```

How far back it goes depends on each cluster's retention settings, so this may return only
recent jobs.
