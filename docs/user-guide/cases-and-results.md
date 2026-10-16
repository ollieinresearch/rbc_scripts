# `cases.csv` and `results.csv`

Two tables drive every figure. They join on `(dim, log_ra, log_pr)`.

## `cases.csv`: what should exist

`data/cases.csv` is edited by hand. One row per parameter point, whatever resolutions
were run there:

```text
dim, log_ra, log_pr, status,   res_check, nz, notes
3,   7,      -1.5,   finished, 1,         ,
3,   8,      -2.5,   running,  0,         , restarted from 6e-3 run
```

| Column | Meaning |
|---|---|
| `dim` | 3 (2 for the earlier two-dimensional runs) |
| `log_ra`, `log_pr` | $\log_{10} Ra$ and $\log_{10} Pr$, multiples of 0.5 |
| `status` | `finished`, `running` or `planned` |
| `res_check` | 1 if a {term}`resolution study` was done at this point |
| `nz` | optional: pins the {term}`production run` to this $N_z$ |
| `notes` | free text |

Lines starting with `#` are ignored, so a point can be switched off without deleting it.
{func}`rbc.cases.load_cases` checks for unknown statuses and duplicated points.

Nothing else needs editing when points are added: slices are found automatically as
every line through at least three points ({func}`rbc.cases.slices`), and the parameter
map ({mod}`rbc.plots.param_grid`) draws filled circles for finished points, open circles
for running ones and stars for resolution studies.

## `results.csv`: what was computed

Written by {mod}`rbc.collect`, one row per run directory:

```bash
python -m rbc.collect                     # all runs under $DATA_ROOT
python -m rbc.collect --recalculate       # redo every run's analysis (e.g. after changing AVG_TIME)
python -m rbc.collect /project/def-goluskin/ollie/ollie_rb_data/3d   # from the project backup
```

For each run it reuses `outputs/*.json` if they are newer than the data, and otherwise
runs the analysis first. Runs are processed in parallel (`--jobs`). Paths matching
`--exclude` (default `_old` or `test`) are skipped.

The table is also written as `results.parquet` when `pyarrow` is importable: the same data
with column types preserved, convenient for pandas (`pd.read_parquet`).

### Production runs

Where a point has several resolutions, one row is flagged `production`: the $N_z$ pinned in
`cases.csv`, else the finest run marked `converged`, else the finest run. Figures use
production rows unless they need all of them (the resolution study).

### Combining clusters

Each row records the `cluster` it was collected on, and collecting keeps rows from other
clusters. To combine Fir and Trillium runs, copy the other cluster's table and merge it:

```bash
scp fir:/path/to/rbc_scripts/data/results.csv results_fir.csv     # from Trillium, say
python -m rbc.collect --merge results_fir.csv
```

(results-columns)=
### Columns

| Column | Meaning |
|---|---|
| `cluster`, `path` | where and which run directory |
| `dim`, `log_ra`, `log_pr`, `label` | from the path |
| `nx`, `ny`, `nz` | grid (modes), from the files |
| `t0`, `t1`, `T_avg`, `n_samples` | averaging window and its length |
| `Nu`, `Nu_se` | Nu from $\langle wT\rangle$ and its batch-means standard error |
| `Nu_gradT`, `Nu_vort` | Nu from $\langle\lvert\nabla T\rvert^2\rangle$ and from $\langle\lvert\omega\rvert^2\rangle$ |
| `Nu_plates`, `plate_imbalance_z` | mean of the two plate fluxes; their difference in standard errors |
| `section_spread_pct`, `qoi_spread_pct` | largest spread of $Nu-1$ over thirds of the window, and over the three identities (%) |
| `converged` | sections within 1% for all three identities (the original rule) |
| `Nu_se_tau`, `tau_int`, `n_eff` | standard error from the autocorrelation time, the time itself, effective samples |
| `p_sections`, `drift_z` | section test and trend over the window (worst of the three identities) |
| `t0_mser` | MSER-5 estimate of the end of the transient |
| `converged_stat` | the statistical rule ({doc}`../theory/convergence`) |
| `Re`, `Re_se`, `Re_rms` | Reynolds number from the volume average and from the profiles |
| `delta_T`, `delta_T_rms`, `delta_u`, `delta_u_slope` | boundary-layer thicknesses ({doc}`../theory/boundary-layers`) |
| `n_delta_*` | Chebyshev nodes inside each layer |
| `Nu_flux_spread_pct` | spread of the total heat flux over $z$ (%) |
| `D_<field>_<panel>` | spectral decades of the time-mean spectrum; field `vel` or `temp`, panel `nz`, `h`, `perp`, `3d`, `meanz` ({doc}`../theory/resolution`) |
| `Dmin_<field>_<panel>` | smallest decades of any single write |
| `n_jobs`, `t_run`, `iterations`, `dt_mean` | job history: jobs, time simulated, steps, mean timestep |
| `wall_hours`, `ranks`, `core_hours`, `core_hours_source` | cost; the source says whether it is logged or estimated |
| `production` | the row used for this point's results |
