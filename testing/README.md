# Resolution testing

Scripts for testing how well Nu is resolved, kept out of the main `rbc` package until
they are validated. Run them from the repository root (`python -m testing.<module>`);
they only read `outputs/nusselt.json`, `outputs/resolution.json` (optional) and
`params.conf`, so they are light enough for a login node.

## The idea

At a parameter point with several resolutions, the finest runs settle on one value of
Nu (the **plateau**). How far each run is from that value, compared with its statistical
noise, says whether it is resolved. Over many points, the calibration asks which
resolution measure (k_max η, k_max η_T, decades, decades of n²E, the vorticity identity)
separates resolved runs from the rest with one threshold that holds at every (Ra, Pr).
That threshold can then be applied to runs without a resolution study.

## Steps

1. **Every run has `outputs/nusselt.json`** (and, for the decades measures,
   `outputs/resolution.json`). Both are written by `post.sh` / `python -m rbc.collect`;
   for one run: `python -m rbc.analysis.nusselt RUN --t0 auto`, then
   `python -m rbc.analysis.resolution RUN`.
2. **One point** (login node):

   ```bash
   python -m testing.convergence $SCRATCH/ollie_rb_data/3d/ra1e6/pr1e-2
   ```

   writes `ra1e6/pr1e-2/outputs/testing/convergence.{txt,pdf,json}`: the plateau, the
   converged Nu and each run's difference from it, with its class (good, likely good,
   likely bad, bad).
3. **Calibration** over every point with two or more runs (login node):

   ```bash
   python -m testing.calibrate $SCRATCH/ollie_rb_data/3d --out $SCRATCH/resolution_testing
   ```

   writes `calibration.{txt,csv,pdf,json}` there (and step 2's files for every point).
   Read `calibration.txt`: measures are listed best first, with their thresholds overall
   and per point. A good measure rejects few runs that were actually resolved, and has
   similar per-point thresholds.
4. **Check other runs** against the thresholds (login node):

   ```bash
   python -m testing.check $SCRATCH/ollie_rb_data/3d/ra1e7 --calibration $SCRATCH/resolution_testing/calibration.json --out $SCRATCH/resolution_testing/check
   ```

   writes `check.{txt,csv,pdf}`. A margin above 1 passes. Runs outside the calibrated
   (Ra, Pr) range are marked `*`: their verdict is an extrapolation.

Options: `--tol` (percent of Nu−1, default 1), `--alpha` (plateau test level, default
0.05), `--include-variants` (keep `_rk443` / `_safety2` runs; off by default),
`--no-plot`; `check` also takes `--threshold strict|safe|bad`.

## Definitions

- **Plateau**: sort the runs by grid size; the plateau starts at the coarsest run from
  which every set "this run and all finer ones" passes a χ² test of one common Nu at
  level `alpha`. Its inverse-variance weighted mean is the converged Nu.
- **Difference**: δ = Nu − Nu_ref and σ = √(SE² + SE_ref²), in % of Nu−1. A plateau
  run is compared with the mean of the *other* plateau runs. Bound = |δ| + 2σ.
- **Classes**: good if bound ≤ tol; bad if |δ| − 2σ > tol; otherwise likely good or
  likely bad according to |δ| ≤ tol.
- **Thresholds** (calibrate): the best value of a measure among likely-bad and bad runs
  (`t_safe`, default), among bad runs only (`t_bad`), or among every run not shown good
  (`t_strict`).
- **k_max η** = πη/Δ, with η = (Pr²/(Ra(Nu−1)))^(1/4) from the exact dissipation and Δ
  the largest spacing (L/N horizontally, the Chebyshev mid-plane gap vertically);
  k_max η_T uses η Pr^(−3/4).
