"""How many frames to render at once.

The visualizations render one frame per worker process. The right number of workers is
set by three limits, and the smallest wins:

* the frames themselves (no point in more workers than frames);
* the cores this job may use (``os.sched_getaffinity``, which follows Slurm's allocation;
  rendering a frame is single-threaded Python and Matplotlib);
* memory: the memory the job may still use (its cgroup limit, or what the node has free),
  times a safety fraction, divided by an estimate of one frame's peak memory.

``--jobs 0`` (the default of every visualization command) chooses this way and prints
the reasoning; a positive ``--jobs`` is used as given, but never more than the number of
frames.
"""
from __future__ import annotations

import os
from pathlib import Path

GB = 1024 ** 3
#: Fraction of the available memory the workers may use together.
MEMORY_FRACTION = 0.7


def available_cpus() -> int:
    """Cores this process may run on (the Slurm allocation inside a job).

    Returns
    -------
    int
    """
    try:
        return len(os.sched_getaffinity(0))
    except AttributeError:                       # not Linux
        return os.cpu_count() or 1


def _cgroup_limit() -> float | None:
    """Memory limit of this process's cgroup in bytes, or None if there is none."""
    candidates = []
    try:
        for line in Path("/proc/self/cgroup").read_text().splitlines():
            rel = line.split(":", 2)[2]
            candidates += [Path("/sys/fs/cgroup") / rel.lstrip("/") / "memory.max",             # v2
                           Path("/sys/fs/cgroup/memory") / rel.lstrip("/") / "memory.limit_in_bytes"]  # v1
    except (OSError, IndexError):
        pass
    candidates += [Path("/sys/fs/cgroup/memory.max"), Path("/sys/fs/cgroup/memory/memory.limit_in_bytes")]
    for c in candidates:
        try:
            v = c.read_text().strip()
        except OSError:
            continue
        if v and v != "max" and int(v) < 2 ** 60:
            return float(v)
    return None


def available_memory() -> float:
    """Memory this job can still use, in bytes.

    The smaller of the node's available memory (``MemAvailable``) and the cgroup limit
    Slurm sets for the job (``--mem``; ``--mem=0`` means the whole node).

    Returns
    -------
    float
    """
    free = None
    try:
        for line in Path("/proc/meminfo").read_text().splitlines():
            if line.startswith("MemAvailable:"):
                free = float(line.split()[1]) * 1024
    except OSError:
        pass
    if free is None:
        free = float(os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_AVPHYS_PAGES"))
    lim = _cgroup_limit()
    return min(free, lim) if lim else free


def auto_jobs(n_tasks: int, mem_per_task: float, requested: int = 0, what: str = "frames") -> int:
    """Number of worker processes for ``n_tasks`` independent tasks.

    Parameters
    ----------
    n_tasks : int
        Number of tasks (frames, snapshots).
    mem_per_task : float
        Estimated peak memory of one task, in bytes.
    requested : int, optional
        A positive value is used as given (capped at ``n_tasks``); 0 chooses automatically.
    what : str, optional
        Name of the tasks, for the printed explanation.

    Returns
    -------
    int
        At least 1.
    """
    n_tasks = max(1, int(n_tasks))
    if requested > 0:
        return min(requested, n_tasks)
    cpus = available_cpus()
    mem = available_memory()
    by_mem = max(1, int(MEMORY_FRACTION * mem // max(mem_per_task, 1)))
    jobs = max(1, min(n_tasks, cpus, by_mem))
    limit = {n_tasks: what, cpus: "cores", by_mem: "memory"}[jobs] if jobs in (n_tasks, cpus, by_mem) else "?"
    print(f"{jobs} workers for {n_tasks} {what} (limited by {limit}: {cpus} cores, "
          f"{mem / GB:.0f} GB available, about {mem_per_task / GB:.2f} GB per task)", flush=True)
    return jobs


def worker_threads(jobs: int) -> int:
    """Cores left for each worker's own threads (software OpenGL renders with several).

    Parameters
    ----------
    jobs : int
        Number of workers.

    Returns
    -------
    int
    """
    return max(1, available_cpus() // max(jobs, 1))


def render_all(fn, items, jobs, retries=2, spawn=False, tasks_per_child=None, before_retry=None) -> list:
    """Run ``fn(args)`` for every item whose output file is missing, in parallel.

    Each item is ``(path, args)``; ``fn`` must write ``path`` (best under a temporary name,
    renamed when complete, so a killed worker leaves nothing that looks finished). Items
    whose file exists are skipped, so an interrupted movie resumes where it stopped. If a
    worker dies (killed for memory, or a lost X connection), the pool is rebuilt and the
    missing items are tried again, up to ``retries`` more rounds.

    Parameters
    ----------
    fn : callable
        Module-level function taking ``args`` (it is pickled to the workers).
    items : list of tuple
        ``(path, args)`` pairs.
    jobs : int
        Worker processes.
    retries : int, optional
        Extra rounds for items that failed.
    spawn : bool, optional
        Start workers as fresh interpreters (needed with OpenGL state); default fork.
    tasks_per_child : int, optional
        Replace each worker after this many items (Python 3.11+; implies spawn).
    before_retry : callable, optional
        Called before each extra round (to restart a virtual display, say).

    Returns
    -------
    list of pathlib.Path
        Outputs still missing after the last round (empty on success).
    """
    import multiprocessing as mp
    import sys
    from concurrent.futures import ProcessPoolExecutor, as_completed
    from concurrent.futures.process import BrokenProcessPool

    kw = {}
    if spawn or tasks_per_child:
        kw["mp_context"] = mp.get_context("spawn")
    if tasks_per_child and sys.version_info >= (3, 11):
        kw["max_tasks_per_child"] = tasks_per_child

    def missing():
        return [(p, a) for p, a in items if not Path(p).exists()]

    for attempt in range(retries + 1):
        todo = missing()
        if not todo:
            if not attempt:
                print(f"all {len(items)} already rendered", flush=True)
            break
        if attempt:
            print(f"round {attempt + 1}: {len(todo)} still missing", flush=True)
            if before_retry:
                before_retry()
        else:
            print(f"rendering {len(todo)} of {len(items)} ({len(items) - len(todo)} already done)", flush=True)
        try:
            with ProcessPoolExecutor(max(1, min(jobs, len(todo))), **kw) as ex:
                futs = {ex.submit(fn, a): p for p, a in todo}
                for f in as_completed(futs):
                    try:
                        f.result()
                    except BrokenProcessPool:
                        raise
                    except Exception as e:
                        print(f"{Path(futs[f]).name}: {e!r}", flush=True)
        except BrokenProcessPool as e:
            print(f"a worker died ({e}); rebuilding the pool", flush=True)
    return [Path(p) for p, _ in missing()]


def make_movie(fdir, pattern, out, fps) -> Path | None:
    """Assemble numbered PNG frames into an MP4 with ffmpeg (if installed).

    Parameters
    ----------
    fdir : pathlib.Path
        Frame folder.
    pattern : str
        printf pattern of the frame names, e.g. ``frame_%05d.png``.
    out : pathlib.Path
        MP4 to write.
    fps : int
        Frame rate.

    Returns
    -------
    pathlib.Path or None
        ``out``, or None without ffmpeg.
    """
    import shutil
    import subprocess
    if not shutil.which("ffmpeg"):
        print("ffmpeg not found: frames written, no movie")
        return None
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-r", str(fps), "-i", str(Path(fdir) / pattern),
                    "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2", "-pix_fmt", "yuv420p", str(out)], check=True)
    return Path(out)


def prepare_frames(fdir, settings: dict, n_frames: int, redo: bool = False) -> None:
    """Make a frame folder ready for :func:`render_all`.

    Removes half-written frames (``*.part.png``) and frames beyond ``n_frames``, and all
    frames if ``redo`` or if they were made with other settings: the settings are kept in
    ``settings.json`` in the folder, so changing the runs or options of a movie never
    mixes old frames into it.

    Parameters
    ----------
    fdir : pathlib.Path
        Frame folder (created if missing).
    settings : dict
        Everything that changes a frame's content (JSON-serializable after ``str``).
    n_frames : int
        Number of frames of the movie.
    redo : bool, optional
        Remove every existing frame.
    """
    import json
    fdir = Path(fdir)
    fdir.mkdir(parents=True, exist_ok=True)
    for stale in fdir.glob("*.part.png"):
        stale.unlink()
    manifest = fdir / "settings.json"
    new = json.dumps(settings, sort_keys=True, default=str, indent=1)
    old = manifest.read_text() if manifest.exists() else None
    if old is not None and old != new and not redo:
        print(f"{fdir.name}: settings changed since the existing frames were made; rendering all again")
    clear = redo or (old is not None and old != new)
    for f in fdir.glob("frame_*.png"):
        if clear or int(f.stem.split("_")[1]) >= n_frames:
            f.unlink()
    manifest.write_text(new)
