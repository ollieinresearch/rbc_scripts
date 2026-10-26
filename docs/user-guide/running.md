# Running and stopping

## Submitting

Always submit from a {term}`login node`, inside the run directory:

```bash
bash $SCRIPTS_PATH/shell/submit.sh                    # job script chosen from $CC_CLUSTER
bash $SCRIPTS_PATH/shell/submit.sh run_sim.sh --time=03:00:00   # explicit script, extra sbatch options
```

`submit.sh` makes the `jobs/` folder for the logs, submits the job script, then `post.sh`
with `--dependency=afterany:<job>`, so
post-processing starts when the simulation ends whether it finished, stopped early or
failed. Plain `sbatch run_sim.sh` also works; the simulation job then tries to queue the
post-processing itself, which Trillium does not allow ({doc}`clusters`).

(sbatch-lines)=
## The `#SBATCH` lines

| Line | Rule |
|---|---|
| `--nodes`, `--ntasks-per-node` | nodes × tasks per node = `N` = `MESHX × MESHY`. Both clusters have 192 cores per node. |
| `--cpus-per-task` (Trillium) | = `C`. |
| `--time` | At most the cluster's limit. Longer jobs wait longer in the queue; every job ends with a checkpoint, so splitting a run into several jobs costs only the start-up time. |
| `--job-name` | Something like `ra8pr-2.5`, to recognize it in `squeue`. |
| `--account` | The allocation to charge. |
| `--output` | `jobs/job_sim_%j.out`: one log per job in `jobs/`, `%j` being the job number. Slurm does not create the folder, so submit with `submit.sh` (or `mkdir jobs` first). |

## Watching

```bash
squeue -u $USER                       # state of your jobs (PD pending, R running)
squeue -u $USER --start               # estimated start times
tail -f jobs/job_sim_<jobid>.out      # the log
sacct -j <jobid> --format=JobID,Elapsed,State,ExitCode   # after it ends
```

The log reports the iteration, time, timestep and Re every 500 iterations, and how often
the timestep changed. Signs of trouble: Re that is NaN (the run stops itself without a
checkpoint), a timestep shrinking steadily, or a log that stops advancing (see
{doc}`../troubleshooting`).

(graceful-stop)=
## Stopping cleanly

A run should always end with a checkpoint, so the next job continues exactly where it
stopped. Three things end a run cleanly:

1. reaching its end time (`SIM_TIME` after the start);
2. the deadline, `STOP_MARGIN` seconds before Slurm's time limit (`set_deadline` reads the
   limit with `squeue` when the job starts);
3. the signal `SIGUSR1`, sent by hand:

   ```bash
   scancel --signal=USR1 <jobid>
   ```

In the last two cases the simulation finishes its current step, all ranks agree to stop,
and it writes the checkpoint. `stop_reason.txt` records which case it was. Plain
`scancel <jobid>` sends SIGTERM, after which Slurm kills the job about 30 seconds later:
too soon for a large checkpoint, so none is attempted, and the next job restarts from the
previous one.

Why "all ranks agree" matters, and how it is done, is in {ref}`collective-stop`. To test
without any of this, set `GRACEFUL_STOP=0`.

(watchdog)=
## The watchdog

If a job hangs (all ranks waiting on each other forever), it would otherwise burn its
whole allocation. With `WATCHDOG=1200`, if no rank makes progress for 20 minutes, the
first rank on each node writes its Python stack to
`watchdog/rank<rank>_<host>.txt` and the job ends. The stacks show where the run was
stuck, which is usually enough to find the cause ({ref}`fftw-deadlock` is an example).
Set `WATCHDOG_ALL=1` to have every rank write its stack. The files are deleted when a
job ends normally.

## Interactive debugging

For quick tests, an interactive job avoids the queue wait for the batch system, e.g. on
Trillium `debugjob` (see each cluster's page) or `salloc` elsewhere. Then, inside the run
directory:

```bash
source params.conf && source $SETUP_PATH && source $WORKFLOW_PATH
load_modules && determine_ic && srun_sim       # or trill_run_sim
```

Every workflow function can be called the same way.
