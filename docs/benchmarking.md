# Choosing resolution, process mesh and node count for Dedalus RBC runs on Trillium

Oliver Brown (UVic, supervisor David Goluskin) · account rrg-goluskin · October 2026

This document explains how we choose the configuration of our 3D Rayleigh–Bénard
convection runs: the vertical resolution `RES_Z`, the MPI process mesh `MESHX × MESHY`,
and the node count. We make these choices with short benchmark jobs rather than by
assumption. Every rule is derived below, the measurements so far are included, and the
code that implements each step is quoted in §8.

## 1. The simulation and what limits it

- **Solver.** Dedalus v3, pseudo-spectral: Fourier in x and y, Chebyshev in z, 3/2
  dealiasing, RK222 IMEX time stepping, adaptive dt from a CFL condition. The aspect ratio
  is 2, so `RES_HOR = 2 RES_Z`.
- **Grid sizes.** At `RES_Z = 256` the coefficients are 512 × 512 × 256 and the dealiased
  grid is 768 × 768 × 384.
- **Parallel layout.** Dedalus distributes the first two axes over a 2D process mesh
  ("pencils"). Each step does forward and backward transforms, global transposes between
  layouts, and banded LU solves, one per horizontal mode.
- **Launch.** One MPI rank per core, packed 8 per L3 cache with `mpirun --map-by
  ppr:8:l3cache:pe=1 --bind-to core`, on whole 192-core nodes.

## 2. The benchmark itself

Each configuration runs as its own job (`shell/bench_one.sh`) on exactly the nodes its mesh
fills, so no allocated core sits idle. The job makes two launches of
`sim/rayleigh_benard.py`, with the production equations, process mesh and rank placement.

**Launch 1: rehearsal.** The script builds everything and takes 10 steps. Nothing is
recorded. Its purpose is to warm the nodes (libraries and file-system caches), so the
measured launch doesn't pay first-use costs. §3 shows those costs biased our first
measurements.

**Launch 2: a short production job.** Starting from a small random temperature
perturbation, with no restart, the script:

1. builds the **production outputs** read from `params.conf` (checkpoint, snapshots,
   analysis, horizontal averages, in-run spectra, flow statistics) before the first step,
   as a production run does, so step times include the output writes and memory includes
   the output buffers;
2. steps at a **fixed** dt for 15 minutes (`BENCH_SECONDS`), in blocks of about 20 s, each
   between two barriers so its time is the slowest rank's. Rank 0 decides after each
   block whether to continue and broadcasts it, so every rank stops on the same step (the
   same mechanism as Dedalus' own wall-time stop; the script needs the block boundaries
   as well). There is no `mpirun` timeout;
3. leaves the first 5 minutes (`BENCH_DISCARD`) out of the timing: they hold the first LU
   factorisations, FFTW planning and the first write of every handler. The time per step
   is the total time of the remaining blocks over their steps (500 to 1500 steps at the
   resolutions in use). The noise estimate is the spread of the means of the four equal
   parts of that 10 minute window, (max − min) / mean;
4. after every block, reads `MemTotal − MemAvailable` on every node and logs the fullest
   (`bench_mem.csv`). If a node passes 95% of its memory, stepping stops and the result is
   recorded as too small (`mem_stop`), instead of the out-of-memory killer losing the job;
5. forces every handler to write once more (as a job's final checkpoint), then records the
   memory two ways: each process's maximum resident set size summed per node (comparable
   to `sacct` MaxRSS), and `MemTotal − MemAvailable` on each node. The larger of the two,
   on the fullest node, is used;
6. appends one row to a CSV file. The job then deletes the output files it wrote.

The job requests `BENCH_SECONDS` + 6 minutes (21 minutes) by default, for the rehearsal, two
solver builds and the final writes. If that turns out too short, stepping stops early
enough before the allocation ends to save the result (`truncated` in the CSV): the job's
end time comes from Slurm. Billing is for the time used either way, but shorter requests
are easier for the scheduler to backfill.

**Why 15 minutes.** Production memory logs of 24 h jobs show that after 15 minutes a job
has reached 80 ± 2.5% of its peak memory, so the peak can be projected from a short run
(§5.2). The earlier 4 minute benchmark built the outputs only at the end and under-read
production memory about 4× (§6.1). Ten minutes of timed steps also average over many
output writes.

**Why a fixed dt.** Every change of dt makes Dedalus refactorise every LU matrix, so CFL
adaptation would add noise that depends on the flow rather than on the configuration.

**Why a random start, and why Ra and Pr don't matter.** A step costs the same whatever the
field values: the transforms and transposes have fixed sizes, and the solves use
already-factorised matrices. There is no iterative solver whose iteration count depends on
the state. The same holds for Ra and Pr, so a benchmark at a given resolution applies to
every run at that resolution.

**What it does not measure:**
- **Memory growth over hours** directly: it is projected from 15 minutes (§5.2).
- **The timestep.** dt depends on the flow, so dt-dependent comparisons use dt ∝ 1/RES_Z
  (§4.2). Wall-time projections use the mean timestep of an earlier run at the same Ra and
  Pr, scaled the same way.

## 3. Load balance of the process mesh

Dedalus splits a distributed axis of `G` groups over `P` ranks into blocks of `⌈G/P⌉`; the
last ranks get a short block or nothing. A step waits for the fullest rank, so the
load-balance efficiency of one mesh dimension is

    e(P) = min over G of  G / (P ⌈G/P⌉),

over the group counts that dimension may have to split:

- `RES_Z` Fourier cos/sin pairs;
- `3 RES_Z` Fourier grid points;
- `RES_Z` Chebyshev coefficients;
- `1.5 RES_Z` Chebyshev grid points.

The mesh's efficiency is `e(MESHX) e(MESHY)`. A dimension divides all four counts exactly
when it divides `RES_Z/2`: with `RES_Z = 2m`, it must divide both `2m` and `3m`, and
`gcd(2m, 3m) = m`.

**A consequence specific to Trillium.** With one rank per core, a job has
`N = nodes × 192 = nodes × 2⁶ × 3` ranks, so the mesh always contains a factor of 3. A
perfectly balanced mesh is possible only when `RES_Z/2` contains a 3:
- possible: 240, 288, 192, 180;
- impossible: 256 (`256/2 = 2⁷`), 280.

**Evidence that the model holds:**
- **Memory pattern.** A production run at 256 on 5 nodes with mesh 30 × 32 had
  `e(30) = 256/270 = 0.948`: 28 rows of 9 Fourier pairs, one row of 4, one empty row. The
  last node therefore holds `(4·9 + 4)/(6·9) = 0.74` of a full node's matrices. It used
  530 GiB against about 692 GiB on the other nodes (0.77). The small excess over 0.74 is
  fixed per-rank overhead.
- **Step cost.** 30 × 32 on 5 nodes took 4.12 node-seconds per step, against 3.57 for
  128 × 9 on 6 nodes (e = 0.981).

**Mesh orientation: open question.** Our first benchmark ran all configurations in one
allocation, and in each orientation pair the slower run was the one that first used a
fresh node:

| Pair, first run listed first | Slower / faster |
|---|---|
| 9×128 then 128×9 (256; 9×128 was the first job on the 6th node) | 1.17 |
| 24×40 then 40×24 (240; first job of the allocation) | 1.04 |
| 36×48 then 48×36 (288; 36×48 was the first job on the 9th node) | 1.004 |

Run order and orientation are confounded, so we don't yet know whether orientation
matters. Every job now has a rehearsal launch, and stage 1 benchmarks both orientations.
Stage 2 uses whichever orientation won.

**cpus-per-task.** On the same 8 nodes at 256, `C = 3` (16 × 24 mesh, 2 ranks per L3
cache) took 1.087 s per step against 0.512 s for `C = 1` (24 × 64), 2.1× slower. Each `C = 3`
rank did 4× the work in 2.1× the time, so having more memory bandwidth per rank helps, but
not enough to break even.

A production comparison at RES_Z = 288 points the same way:

| Configuration | Per step |
|---|---|
| 8 nodes, C = 2, 24 × 32 | 8 × 0.9786 s = 7.83 node-seconds |
| 10 nodes, C = 1, 48 × 40 | 10 × 0.6526 s = 6.53 node-seconds |

C = 1 was 1.20× cheaper per step, despite its worse balance (0.90 against 0.96). All runs
now use one core per rank.

## 4. Choosing the resolution (stage 1)

The physically required resolution comes first, from our spectral-decay criteria. The
benchmark only chooses between nearby resolutions that we would accept.

### 4.1 Transform sizes

FFTW is fastest for sizes `2ᵃ 3ᵇ 5ᶜ 7ᵈ 11ᵉ 13ᶠ` with `e + f ≤ 1`, and especially for powers
of 2; other sizes fall back to a slower general algorithm [1]. The Chebyshev transforms
are DCT-II/III, whose speed depends on the factorisation of the transform length [2].
Dedalus recommends "powers of two, or powers of two multiplied by other small factors" [3].
The transform lengths `3 RES_Z` and `1.5 RES_Z` have the primes of `RES_Z` plus a 3, so
classifying `RES_Z` classifies every transform:

| Tier | RES_Z | Use |
|---|---|---|
| A | 3-smooth (`2ᵃ3ᵇ`) | preferred |
| B | 5-smooth | acceptable |
| C | one factor of 7 | only if a benchmark shows it pays |
| — | 11, 13 or larger, 7², or odd | rejected |

No published source quantifies the per-point penalty for a 5 or a 7 at these sizes. That
is why the decision rests on measurement rather than on a penalty table.

### 4.2 Cost of a resolution

The work per step is proportional to the number of grid points times a logarithm,
`RES_Z³ log RES_Z`. The advective CFL limit makes `dt ∝ 1/RES_Z`, so the number of steps
to reach a fixed simulation time is `∝ RES_Z`. Diffusion is implicit, so it does not limit
dt. The cost to reach a given simulation time is therefore

    C(N) ∝ N⁴ log N.

Break-even: going from `N₁` down to `N₂` saves `1 − C(N₂)/C(N₁)`. The smaller size is
worth it unless its factorisation or balance penalty on the step time exceeds
`C(N₁)/C(N₂) − 1`:

| Change | (N₂/N₁)⁴ | log ratio | Cost ratio | Saving | Break-even penalty |
|---|---|---|---|---|---|
| 288 → 280 | 0.8934 | 0.9950 | 0.889 | 11.1% | 12.5% |
| 256 → 240 | 0.7725 | 0.9884 | 0.763 | 23.7% | 31% |
| 192 → 180 | 0.7725 | 0.9877 | 0.763 | 23.7% | 31% |

So 192 → 180 or 256 → 240 almost certainly pays, provided the lower resolution is
physically adequate. 288 → 280 is too close to call: 280 adds a 7 and cannot be balanced
on 192-rank nodes. It is the kind of case worth benchmarking.

### 4.3 Which resolutions to benchmark

The lowest acceptable resolution, `RES_FLOOR`, comes from our physics criterion: a short
run whose time-averaged spectra span at least 4 decades. Below the floor nothing is
considered. Above it, a higher resolution can only win if its better factorisation or mesh
balance outweighs its extra work.

Requiring candidates to be "at least X apart" or "at most X apart" doesn't answer that
question. We use a **predicted-cost band** instead: a candidate is worth measuring only if
its predicted cost is within the plausible penalty of the cheapest, so that measurement
could change the answer. The planner:

1. takes every even `RES_Z` from the floor `F` up to `1.125 F`, and discards rejected tiers;
2. gives each one the most balanced node count between its memory minimum (§5.2, from a
   prior per-cell memory before anything is measured) and 1.35× that;
3. predicts its relative cost as `(N/F)⁴ log N / log F × prior / e`. The priors are 1.00,
   1.03 and 1.08 for tiers A, B and C. They are assumptions, used only for this
   preselection;
4. keeps the floor, the best-predicted resolution of each tier, and nothing predicted
   more than 15% above the cheapest; at most 4 resolutions.

For a floor of 256, every higher candidate is predicted at least 35% dearer (e.g. 270, 280,
288), so only 256 is run, in both orientations. For a floor of 192, 216 is predicted 78%
dearer. Going up in resolution almost never pays, because the cost grows as `N⁴`.

**Horizontal and vertical resolution separately.** Sometimes the vertical boundary layers
limit the resolution, not the bulk. A resolution is therefore a pair: `RES_HOR` modes in x
and in y, and `RES_Z` Chebyshev modes. By default `RES_HOR` follows `RES_Z` at the ratio in
`params.conf`. Setting `RES_HOR_FLOOR` in `bench.conf` lets the two vary independently;
`RES_HOR_ABOVE=0` keeps the horizontal fixed and refines only in z. With independent axes:
- the work per step is `RES_HOR² RES_Z log(RES_HOR² RES_Z)`;
- steps per unit time scale as `max(RES_HOR/Lx, RES_Z)`, the finer of the horizontal
  spacing and the mean vertical spacing;
- a pair is classified by the worse tier of its two axes.

With `RES_HOR = Lx × RES_Z`, all of these reduce to the formulas above.

### 4.4 The decision

The cost per unit of simulated time is `node-seconds per step / dt`, and `dt ∝ 1/RES_Z`,
so it is proportional to `nodes × s_per_step × RES_Z`. The lowest-cost resolution at or
above the floor is chosen. If the runner-up is within the measured noise (at least 3%), the
**higher** resolution is taken: equal cost, better resolved. The result is written to
`decision_resolution.txt` with a figure.

## 5. Choosing the node count (stage 2)

### 5.1 Objective

More nodes shorten the wall time but usually raise node-hours. We minimise

    J = (node-hours) × (wall time)^k = N T^(1+k),     k = 1 by default,

which weights a fractional saving in wall time equally with the same fractional saving in
node-hours. It is the analogue of the energy-delay product.

Write `r = N/N₀`, `S = T₀/T` and `E = S/r`, the efficiency relative to the smallest count
tested. Then

    J/J₀ = r · (T/T₀)^(1+k) = r (r E)^-(1+k) = r^-k E^-(1+k),

so a larger count beats the baseline exactly when

    E > r^(-k/(1+k))       (k = 1:  E > r^-1/2).

| Node ratio r | Efficiency needed (k = 1) |
|---|---|
| 1.25 | 0.89 |
| 1.5 | 0.82 |
| 2 | 0.71 |

These are close to common practice. Princeton's scaling guide [6], for example, settles on
two nodes at 90% efficiency and rejects a larger count at 67%.

Example: going from 8 to 12 nodes at E = 0.9 gives `S = 1.5 × 0.9 = 1.35`. That is 26% less
wall time for 11% more node-hours, and `J/J₀ = 0.82`, so it is accepted.

Counts whose J is within the measured noise of the minimum count as ties, and the fewest
nodes among them is chosen. The decision file also reports the marginal efficiency
`(T_prev/T − 1)/(N/N_prev − 1)`, the share of ideal speed-up delivered by the last nodes
added. Minimising J with k = 1 keeps it above about 0.45.

**Faster-than-ideal jumps.** If a step up in node count is faster than ideal by more than
the noise, the decision file warns. The code cannot scale better than ideal (cache effects
give a few percent), so the smaller counts were slowed by something that went away: memory
pressure (a nearly full node without swap evicts program code and reads it again from the
file system) or the network (MPI falling back to TCP). The warning lists the measured memory
on both sides of the jump. Times at the slow counts should not be trusted until it is
explained.

### 5.2 Range swept

**Memory model.** Memory per node is a fixed part (libraries, MPI, Python; about 90 GiB at
192 ranks per node, measured on a 96³-sized problem where the data is negligible) plus a
part that scales with resolution and is spread over the nodes:

    peak per node ≈ M_fixed + V / nodes

`V` comes from stage 1's measurement: with `m` the memory per node at the end of the 15
minute benchmark on `n` nodes and `f = 0.80` the share of a 24 h job's peak reached after 15
minutes (`BENCH_MEM_FRACTION`),

    V = n × (m / f − M_fixed)

so on the benchmarked count the projected peak is `m / f`. The ±2.5% spread of `f` is ±3% on
the projection, inside the 15% headroom. The minimum node count is the smallest that keeps
the projected peak below `747 GiB × 0.85`. 747 GiB is the `sacct` MaxRSS at which a run on
Trillium was killed for running out of memory (another cluster needs its own value), and
0.85 is headroom. Results of the earlier short benchmark are still read, with `MEM_FACTOR`
and the growth allowance `MEM_GROWTH` instead of `f`; the decision files warn when the two
kinds are mixed, because their step times are not comparable.

**History.** Every 15 minute result is added to a history file
(`~/rbc_bench_history.csv`) with its cluster, cpus-per-task and snapshot setting; the
stage 2 choice is marked. A later benchmark on the same cluster uses it twice:
- **stage 1** starts each candidate resolution at the node count its memory needs, with
  `V` scaled from the nearest resolution measured before by the ratio of
  `(RES_HOR/2)² RES_Z` (field and transform buffers scale with the grid; the fixed part
  per node is not scaled). For 288 → 384 this is (384/288)³ = 2.37;
- **stage 2** also sweeps the count chosen at the nearest resolution, scaled by the same
  ratio (work per step scales with the grid size, up to a logarithm), extending the range
  if needed.

Only results with the same cpus-per-task are used: C = 2 used about twice the memory of
C = 1.

**Range.** From that minimum up to 2× it, keeping node counts whose best mesh has balance
of at least `min(0.93, best in range − 0.03)`, at most 6 of them, always including both
ends. A count whose projected end-of-job memory exceeds the budget is never chosen.

**Expected scaling regime.** Burns et al. [4] report Dedalus' parallel efficiency staying
above 50% down to about 8 pencils per core. At 256 on 6 nodes there are 512² / 1152 ≈ 228
pencils per core, so this range is far from that limit.

The result is written to `decision_nodes.txt` with a figure (`decision_nodes.pdf`, vector,
and a 200 dpi `.png`) of four panels, each comparing against the choice:
- (a) step time against nodes, log-log, with perfect scaling and the equal-J curve through
  the choice (a point below that curve would have lower J);
- (b) projected memory per node at the end of a production job: fixed part, fields at the
  start and growth, with the measured benchmark value, the budget and the node's memory;
- (c) J relative to its minimum, with the noise band of a tie, each count's change in
  node-hours and wall time against the choice, and the marginal efficiency of each step;
- (d) node-hours against wall time, each relative to the best feasible count: the
  time-to-solution/cost Pareto plot common in HPC scaling studies, with curves of equal J
  (with k = 1, J is cost × delay, the same form as the energy-delay product used in
  computer architecture).

## 6. Limitations and open questions

### 6.1 Memory

Three observations led to the 15 minute benchmark with outputs built from the start.

- **The compute-only benchmark under-read production by about 3×.** At 256 on 5 nodes it
  peaked at about 1,200 GiB in total. Production at the same configuration already used
  about 3,300 GiB 10 minutes into the run.
- **Two production runs at the same resolution and node count (288, 8 nodes) peaked at 635
  and about 300 GB per node.** Output settings are the first suspect.

- **Growth.** That production run grew from about 3,280 GiB at 10 minutes to about
  3,540 GiB at 3.6 hours, roughly 2.3% per hour; a run at 288 on 8 nodes grew much less
  over 12 hours. A fixed growth rate (`MEM_GROWTH`) could not fit both. The per-node
  memory logs of production runs (`rbc/memguard.py`) then showed that 15 minutes in, a
  24 h job is at 80 ± 2.5% of its peak, which replaces the growth rate.

The 0.80 comes from production jobs that start from a checkpoint and adapt dt; the
benchmark starts from rest at a fixed dt. Its 15 minutes also start after the solver build,
a little later in the job than a production run's, which errs towards a higher reading.
Comparing a production run's first 15 minutes against its benchmark (`bench_mem.csv`)
would check the transfer.

**Measurement.** Summing `ru_maxrss` over a node's ranks counts shared pages once per rank:
for one run it gave 240 GiB against 207 GiB from `sacct`. `MemTotal − MemAvailable`
includes the operating system. The tool takes the larger of the two.

**Guard in production runs.** Nodes without swap give no gradual warning as memory fills:
they run at full speed until close to the limit, then slow sharply (the kernel evicts
program code and libraries and reads them back) and the out-of-memory killer ends the job.
So production runs now check the
memory directly (`rbc/memguard.py`). About every 5 minutes of wall time (at iterations all ranks agree on), each node reads
`MemTotal − MemAvailable`, and all ranks take the maximum over nodes; the fullest node
decides.
- **Above 95% on any node:** the run stops cleanly, writes its checkpoint and records
  `stop_reason` "memory". The next job resumes from that checkpoint in fresh processes, so
  the memory that built up is released.
- **Above 85% in the first 20 minutes:** a warning is logged that the configuration starts
  too close to the limit.

### 6.2 Other gaps

- **Cost model.** It assumes `dt ∝ 1/RES_Z`. If near-wall Chebyshev clustering sets dt
  instead, lower resolutions save even more.
- **Noise.** Each configuration is timed in many 25-step blocks; the spread of the block
  means is the noise estimate used for ties. A configuration can also be repeated.

### 6.3 Questions for SciNet

1. What is the node-failure rate for 48-node, 12-hour jobs? It sets the checkpoint interval
   (currently moving from 60 to 15 minutes of wall time).
2. Is `--requeue` recommended for automatic restarts on Trillium?
3. Is there a recommended rank placement across the CCDs beyond 8 ranks per L3 cache?

## 7. Measurements so far

First benchmark, job 2697353, Ra = 10⁶, Pr = 10⁻⁴, 100 steps at dt = 5 × 10⁻⁴:

| RES_Z | Nodes | Mesh | C | Balance e | s/step | Node-s/step | Peak GiB/node (max/min) |
|---|---|---|---|---|---|---|---|
| 240 | 5 | 24×40 | 1 | 1.000 | 0.634 | 3.17 | 201/198 |
| 240 | 5 | 40×24 | 1 | 1.000 | 0.607 | 3.04 | 196/192 |
| 256 | 5 | 30×32 | 1 | 0.948 | 0.823 | 4.12 | 240/196 |
| 256 | 6 | 9×128 | 1 | 0.981 | 0.695 | 4.17 | 239/216 |
| 256 | 6 | 128×9 | 1 | 0.981 | 0.595 | 3.57 | 190/185 |
| 256 | 8 | 24×64 | 1 | 0.970 | 0.512 | 4.10 | 177/151 |
| 256 | 8 | 16×24 | 3 | 0.970 | 1.087 | 8.69 | 117/114 |
| 288 | 9 | 36×48 | 1 | 1.000 | 0.665 | 5.98 | 200/196 |
| 288 | 9 | 48×36 | 1 | 1.000 | 0.662 | 5.96 | 193/189 |

These runs predate the rehearsal launch (§3), so the first run on each fresh node may be biased slow; the best configuration of each resolution below was in each case a second run.

The stage 1 rule (§4.4) applied to these rows, with every resolution eligible, relative to 256:

```
CHOICE: RES_Z = 240 on 5 nodes, mesh 40x24 (C = 1)
Reason: lowest node-hours per unit of simulated time.

 RES_Z factors      nodes     mesh    eff   s/step  node-s  cost/t  wall/t   work
   240 2^4 3 5          5    40x24  1.000   0.6073    3.04   0.797   0.957  1.044  <-
   256 2^8              6    128x9  0.981   0.5949    3.57   1.000   1.000  1.000
   288 2^5 3^2          9    48x36  1.000   0.6622    5.96   1.878   1.252  1.148
```

Here `cost/t` is node-hours per unit of simulated time and `wall/t` wall time per unit of
simulated time, both relative to 256. `work` is node-seconds per `RES_Z³ ln RES_Z`,
relative to 256.

The cost model of §4.2, with the balance factors of §3, gets the ordering right and most
of the magnitude:

| RES_Z | Predicted cost/t | Measured |
|---|---|---|
| 240 | 0.7635 × 0.981 = 0.749 | 0.797 |
| 288 | 1.636 × 0.981 = 1.605 | 1.878 |

The `work` column is the residual: per unit of work, 240 ran 4% and 288 ran 15% slower than
256. 288 is 3-smooth and perfectly balanced, so its residual is more likely communication
or cache effects at the larger size and node count than factorisation. Residuals of this
size are why the decision uses measurements and treats the model only as a preselection.

With a floor of 256, the choice is 256 on 6 nodes with mesh 128 × 9, which is 13% fewer
node-seconds per step (and 28% less time per step) than the 30 × 32 mesh we had been
using. The stage 2 node sweep is the next job; its results will follow.

## 8. The code in action

**Commands.** Everything runs from a benchmark folder holding `params.conf` (physics,
output settings and paths) and `bench.conf` (what to benchmark: the floor resolution,
timing and memory settings), on a login node:

```bash
python -m rbc.bench          # stage 1: shows the plan, asks, submits one job per configuration
python -m rbc.bench          # later: reports progress; once done, decides the resolution,
                             #        plans stage 2 and asks before submitting it
python -m rbc.bench          # after stage 2: decides the node count, prints the production lines
python -m rbc.bench status   # progress only
```

**One configuration per job, on exactly its nodes** (`shell/bench_one.sh`, submitted
with `--nodes`, `--ntasks-per-node` and `--time` set by `rbc.bench`): a rehearsal launch,
then the measured one.

```bash
mpirun -n $(( mx * my )) --bind-to core --map-by "ppr:$(( 8 / c )):l3cache:pe=$c" \
    -x HDF5_USE_FILE_LOCKING -x OMP_NUM_THREADS -x RBC_WATCHDOG \
    python3 "$SCRIPTS_PATH/sim/rayleigh_benard.py" \
    --Ra "$RA" --Pr "$PR" --Lx "$LX" --Ly "$LY" --nz "$res" --meshx "$mx" --meshy "$my" \
    --stepper "$STEPPER" --dt "$DT" --sim_time 1e9 --basepath "$run" \
    --bench "${BENCH_ITERS:-100}" --bench_csv "$run/bench.csv" < /dev/null
exit $?
```

**The timed loop** (`run_bench` in `sim/rayleigh_benard.py`):

```python
for i in range(BENCH_WARMUP):                 # 20 steps: first LU factorisation, FFTW plans
    solver.step(dt)
sizes = [iters // BENCH_BLOCKS + (b < iters % BENCH_BLOCKS) for b in range(BENCH_BLOCKS)]
block_means = []
for n in (k for k in sizes if k):
    comm.Barrier()                            # every block timed between barriers:
    t0 = time.time()                          # the slowest rank's time
    for i in range(n):
        solver.step(dt)                       # fixed dt: no refactorisation
    comm.Barrier()
    block_means.append((time.time() - t0) / n)
loop = sum(m * n for m, n in zip(block_means, (k for k in sizes if k)))
spread = (max(block_means) - min(block_means)) / (loop / iters)
rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20      # GiB, per rank
node_rss = node.allreduce(rss, op=MPI.SUM)                            # summed per node
```

**Load balance** (`rbc/plan.py`):

```python
def groups(res_z, gamma=2.0):
    nh = round(gamma * res_z)
    return nh // 2, 3 * nh // 2, res_z, 3 * res_z // 2

def efficiency(p, sizes):
    return min(g / (p * math.ceil(g / p)) for g in sizes)
```

**The node decision** (`decide_nodes` in `rbc/bench.py`):

```python
for p in pts:
    p["r"] = p["nodes"] / base["nodes"]
    p["S"] = base["s"] / p["s"]
    p["E"] = p["S"] / p["r"]
    p["thr"] = p["r"] ** (-k / (1 + k))          # break-even efficiency against the baseline
    p["J"] = p["nodes"] * p["s"] ** (1 + k) / (base["nodes"] * base["s"] ** (1 + k))
jmin = min(p["J"] for p in pts)
top = min(pts, key=lambda p: p["J"])
ties = [p for p in pts if p["J"] <= jmin * (1 + (1 + k) * noise_of(p, top))]
best = min(ties, key=lambda p: p["nodes"])      # ties within noise go to fewer nodes
```

## Sources

1. FFTW 3 manual, Real-data DFTs: https://www.fftw.org/fftw3_doc/Real_002ddata-DFTs.html
2. FFTW 3 manual, Real-to-Real Transforms: https://www.fftw.org/fftw3_doc/Real_002dto_002dReal-Transforms.html
3. Dedalus documentation, Performance tips: https://dedalus-project.readthedocs.io/en/latest/pages/performance_tips.html
4. K. J. Burns, G. M. Vasil, J. S. Oishi, D. Lecoanet, B. P. Brown, "Dedalus: A flexible framework for numerical simulations with spectral methods", Phys. Rev. Research 2, 023068 (2020): https://link.aps.org/doi/10.1103/PhysRevResearch.2.023068
5. SciNet, Trillium Quickstart: https://docs.scinet.utoronto.ca/index.php/Trillium_Quickstart
6. Princeton Research Computing, Scaling analysis: https://researchcomputing.princeton.edu/support/knowledge-base/scaling-analysis
