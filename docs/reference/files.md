# Output files

Everything a run writes, in its {term}`run directory`. HDF5 files from Dedalus share one
layout: `scales/` holds per-write metadata (`sim_time`, `iteration`, `timestep`,
`wall_time`, `write_number`) and grids, and `tasks/` holds one dataset per output, with
time as the first axis. Read them with {func}`rbc.io.load_series`, which joins sets across
restarts correctly.

## Written by the simulation

`analysis/analysis_s<N>.h5`, every `A_FREQ` (volume averages, scalars):

| Task | Quantity |
|---|---|
| `avg_wT` | $\langle wT\rangle$ |
| `avg_grad_T_sq` | $\langle\lvert\nabla T\rvert^2\rangle$ |
| `avg_vorticity_sq` | $\langle\lvert\boldsymbol\omega\rvert^2\rangle$ |
| `avg_K` | $\langle\lvert\boldsymbol u\rvert^2\rangle$ (twice the kinetic energy) |
| `Ra`, `Pr` | the parameters, for files without a parseable path |

`horizontal_analysis/horizontal_analysis_s<N>.h5`, every `A_FREQ` (horizontal averages,
profiles in $z$ on the dealiased Chebyshev grid, attached as a dimension scale):

| Task | Quantity |
|---|---|
| `avg_T` | $\langle T\rangle_h$ |
| `avg_u_sq`, `avg_v_sq`, `avg_w_sq` | $\langle u^2\rangle_h$, $\langle v^2\rangle_h$, $\langle w^2\rangle_h$ |
| `avg_T_sq` | $\langle T^2\rangle_h$ |
| `havg_wT` | $\langle wT\rangle_h$ |

`state/state_s<N>.h5` (+ folder of process files in `virtual` mode): the full solution,
in coefficient space by default; one write per set.

`snapshots/snapshots_s<N>.h5`, every `SNAPS_FREQ` with `SNAPSHOTS=1`, 30 writes per set:
`temperature`, `u` (3 components), `vorticity` ($\lvert\boldsymbol\omega\rvert^2$), on the
dealiased grid.

`spectra/spectra_s<N>.h5`, every `SPEC_FREQ`: groups `vel/` and `temp/` with datasets
`E_kx`, `E_ky`, `E_nz`, `E_perp`, `E_3d`, `mean_E_nz` (time × mode); `modes/` with the
mode-number axes and shell counts; attributes with the grid and box. See
{mod}`rbc.spectra.core`.

`outputs/jobs.csv`: one row per job:

| Column | Meaning |
|---|---|
| `finished` | date and time the job ended |
| `job_id`, `cluster`, `nodes` | from Slurm |
| `ranks`, `cpus_per_rank`, `mesh` | parallel layout |
| `nx`, `ny`, `nz` | grid (modes) |
| `t_start`, `t_end`, `iter_start`, `iter_end` | stretch of simulation covered |
| `wall_s`, `setup_s`, `loop_s` | wall-clock seconds: total, before the loop, in the loop |
| `core_hours` | `wall_s × ranks × cpus_per_rank / 3600` |
| `stop_reason` | `finished`, `deadline`, `SIGUSR1`, ... |

`stop_reason.txt`: the last job's stop reason. `watchdog/`: stack dumps, kept only if a job
hung.

## Written by the analysis

`outputs/nusselt.json` (and the readable `info.txt`), from {func}`rbc.analysis.nusselt.compute`:
`Ra`, `Pr`, the window (`t_start`, `t_end`, `t0`, `t1`, `T_avg`, `n_samples`), `t0_mser`,
`criteria`, `convergence_qoi`; per identity (`avg_wT`, `avg_grad_T_sq`, and
`avg_vorticity_sq` if it was written) and for `Re`:
the value, `sections`, `spread_pct`, `se`, `se_tau`, `tau_int`, `n_eff`, `rel_se`,
`p_sections`, `drift_pct`, `drift_z`; `plates` (`Nu_bot`, `Nu_top`, `imbalance_pct`,
`imbalance_z`); `Nu`, `Nu_se`, `qoi_spread_pct`, `vort_diff_pct`, `converged`,
`converged_stat`.

`outputs/profiles.npz`, from {func}`rbc.analysis.profiles.compute`: arrays `z` (centred,
$-1/2 \dots 1/2$), `T`, `w_rms`, `uh_rms`, `u2`, `theta_rms`, `nu_z`.

`outputs/profiles.json`: `Ra`, `Pr`, `nz`, the window, `Nu_wall_bot`, `Nu_wall_top`,
`Re_rms`, the boundary-layer thicknesses and node counts, `Nu_flux_mean`,
`Nu_flux_spread_pct`.

`outputs/spectra.json`, from {func}`rbc.spectra.cumulative.write_json`: the window, the
grid, and `D_`, `Dmin_`, `pileup_` for each field and panel, and `D_<field>_meanz`.

`outputs/spec_range.txt`: the same decades as a table, with the time statistics of $D(t)$.

Figures: `info.pdf`, `profiles.pdf`, `cumulative_spectra_velocity.pdf`,
`cumulative_spectra_temperature.pdf`, `decades_timeseries.pdf`; in the parameter folder,
`outputs/spectra_compare.pdf`; renders and movies in `visualization/`.

## Tables in the repository

`data/cases.csv` and `data/results.csv`: {doc}`../user-guide/cases-and-results`.
