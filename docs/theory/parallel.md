# Parallel runs: stopping, I/O and hangs

Production runs use thousands of MPI ranks. Most problems at that scale come from one rank
doing something the others do not, which leaves them all waiting forever, or from the
filesystem being asked the same question thousands of times at once. This page explains
the mechanisms the code relies on and the problems they solve; the symptoms are in
{doc}`../troubleshooting`.

## Collective operations

Dedalus distributes each field over a 2D {term}`process mesh` and moves data between
layouts (grid space, coefficient space, and the transposes between them) with collective
MPI calls: every rank of a communicator must make the same call in the same order. A rank
that skips one, or makes a different one, leaves the others blocked inside MPI, consuming
their allocation and printing nothing. Each fix below is a way of keeping every rank on the
same path.

(collective-stop)=
## Graceful stop

**Problem.** A signal (or the job's deadline) reaches the ranks at slightly different
times. A rank that leaves the time loop on its own blocks the others in the next transpose,
and the final checkpoint, itself a collective write, never happens.

**Mechanism** ({class}`rbc.shutdown.GracefulStop`). Signal handlers only set a local flag
(calling MPI or HDF5 inside a handler is unsafe). Every 10 iterations all ranks combine
their flags with one `Allreduce(MAX)`. If any rank has a flag, every rank sees it on the
same iteration, leaves the loop together, and writes the checkpoint collectively. The
deadline check is part of the same flag, using a time computed once from Slurm's end time
(`set_deadline`). SIGTERM is treated differently: it means Slurm will kill the job in about
30 seconds, too soon for a large checkpoint, and a checkpoint cut off midway would become
the newest, corrupt, restart file; so none is written.

The cost is one tiny collective every 10 steps, negligible against the transposes of a
step.

(watchdog-mechanism)=
## Watchdog

A hang shows nothing in the log. The {ref}`watchdog <watchdog>` uses Python's
`faulthandler`, whose timer runs in a C thread and can therefore dump the stack of a thread
blocked inside MPI or HDF5, which a Python-level timer cannot. It is re-armed every 10
iterations; the stage being armed (restart load, first step, final checkpoint) is written
to the same file, so a dump also says how far the run got. One rank per node dumps (to keep
the number of files small), and those ranks then exit, which ends the job.

(handler-modes)=
## File-handler modes

Dedalus writes output in one of three ways (`PARA`, `PAR_SMALL`):

| Mode | How | Good for |
|---|---|---|
| `virtual` | every rank writes its own process file in a set folder; a small joint file maps them into one dataset (HDF5 virtual datasets) | large output (checkpoints, snapshots): fast, no coordination |
| `mpio` | all ranks write one shared file with collective MPI-IO | when creating thousands of files is the problem (Lustre metadata storms) |
| `gather` | data is gathered to rank 0, which writes one file | small output (`analysis`, `horizontal_analysis`) |

(vast-hang)=
## The VAST hang at file creation

**Problem.** On Trillium's {term}`VAST` filesystem (NFS), a node can cache the answer "this
file does not exist" for a short time after another node has created it. Dedalus asks, on
every rank and at every write, whether the current set file exists, and creates a new set
inside MPI barriers. A rank whose node holds a stale "does not exist" takes a different
branch from the others, and the job deadlocks. A related race: each rank creates its
process file inside a folder that rank 0 has just made, which other nodes may still see as
missing.

**Earlier workaround.** `fixed_get_file.py`, a wrapper by Mike, made rank 0 decide whether
the file exists and broadcast the answer.

**Fix** ({mod}`rbc.vast`, applied by default). No rank asks the filesystem at all: each
handler remembers which set number it has created, which is the same on every rank by
construction. That removes $N$ `stat` calls per write as well as the race. The set folder
is then created once per node (rank 0 of each node calls `mkdir` with `exist_ok`, which
refreshes that node's cache), followed by a barrier within the node, before any rank
creates its process file. It is harmless on Lustre, so it is on everywhere; `VAST_PATCH=0`
turns it off.

(n2-opens)=
## Slow restarts: N² file opens

**Problem.** A `virtual` checkpoint's time axis is unlimited, so to read any part of it
HDF5 must open every process file to learn its size. When every rank loads its own slab
through the joint file, $N$ ranks open $N$ files each: for 5184 ranks, about
$2.7\times10^7$ opens, which kept large runs in their restart load for a long time.

**Fix** ({func}`rbc.vast.load_state`). Rank 0 reads the joint file's virtual mapping,
which records which block of the global array each process file holds, without opening any
process file. Each rank then opens only the process files that overlap its own slab and
assembles the slab from them: one file when the process mesh is unchanged, a few when it
changed. Dedalus' own loader is used, on every rank together, only for a plain file (a
regridded restart), a change of resolution, a mapping that leaves part of some rank's slab
uncovered, or an error while reading; the log names the reason. `--verify_load` compares
the two loaders on a small case.

(fftw-deadlock)=
### The restart deadlock it once caused

The first version of the fast loader skipped fields for which a rank holds no data (tau
fields stored in grid space are empty on most ranks). Those ranks left the fields in their
default layout while the others changed it, so on the first timestep the ranks disagreed
about which fields to transpose together. The transposes build their FFTW plans
collectively, and the job hung inside the plan creation; the watchdog's stacks all ended in
Dedalus' `distributor.py`, in `_plan`. Every rank now sets every field's layout and scales,
with or without data.

(export-env)=
## Environment on other nodes (Trillium)

Trillium's `sbatch` implies `--export=NONE`, which Slurm records as `SLURM_EXPORT_ENV=NONE`
in the job. Every `srun` inside the job inherits that, including the one `mpirun` uses to
start its daemons on the other nodes, so ranks there would start without the modules,
virtual environment and variables loaded by the job script. `workflow.sh` sets
`SLURM_EXPORT_ENV=ALL`
([SchedMD bug 6772](https://support.schedmd.com/show_bug.cgi?id=6772)), and
`trill_run_sim` forwards the variables the simulation reads with `-x`.

## Rank placement

On Trillium, ranks are placed with `--map-by ppr:8:l3cache:pe=1`: eight ranks per
{term}`L3 cache`, each bound to its own core. Spectral codes spend much of their time in
FFTs and transposes, which are limited by memory bandwidth, and spreading ranks evenly
over the caches uses all of it. With Open MPI, check the placement with
`--report-bindings`; a run whose ranks are unbound or share cores is noticeably slower.
