"""Choose the restart checkpoint, and keep or prune old checkpoints.

Prints ``<index> <sim_time> <path>`` for ``determine_ic`` in ``shell/workflow.sh`` and
points ``restart/restart.h5`` at the chosen set. State sets that cannot be read (for
example from a job killed mid-write) are skipped. Usage::

    python -m rbc.restart RUN_DIR [--time T]     # T = 0 -> latest write
    python -m rbc.restart RUN_DIR --from OTHER_RUN --shape NX NY NZ
        seed a new run from another run's checkpoint, regridded into
        restart/restart.h5 if its (coefficient) shape differs. Once RUN_DIR has
        checkpoints of its own, --from is ignored, so RESTART_FROM can stay set.
        A change of resolution is refused within a run: make a new run folder.
    python -m rbc.restart RUN_DIR --keep 2 [--prune]
        print the 2 newest readable state sets (paths of set folder and .h5); with
        --prune, delete every older set. Nothing is deleted unless 'keep' readable
        sets exist, and the set restart/restart.h5 points to is never deleted.
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import h5py
import numpy as np

from rbc.io import set_files


def readable(path: Path) -> np.ndarray | None:
    """Cheap integrity check of a state set.

    For virtual sets every process file must exist (a missing one reads back as
    zeros, not as an error); for other sets one element of each task is read.

    Parameters
    ----------
    path : pathlib.Path
        ``state_s<N>.h5``.

    Returns
    -------
    numpy.ndarray or None
        The set's write times, or None if it is empty or unreadable.
    """
    try:
        with h5py.File(path, "r") as f:
            t = f["scales/sim_time"][:]
            if t.size == 0:
                return None
            have = None
            for name in f["tasks"]:
                d = f["tasks"][name]
                if d.is_virtual:
                    need = {Path(s.file_name).name for s in d.virtual_sources()}
                    if have is None:
                        have = {p.name for p in (path.parent / path.stem).iterdir()}
                    if need - have:
                        return None
                else:
                    d[(-1,) + (0,) * (d.ndim - 1)]
            return t
    except Exception:
        return None


def choose(run: Path, time: float = 0.0):
    """The checkpoint to restart from.

    Parameters
    ----------
    run : pathlib.Path
        Run directory with a ``state/`` folder.
    time : float, optional
        Simulation time to restart from (0: the newest readable write).

    Returns
    -------
    index : int
        Write index within the set.
    sim_time : float
        Its simulation time.
    path : pathlib.Path
        The set file.

    Raises
    ------
    FileNotFoundError
        If no set is readable.
    """
    for p in reversed(set_files(run / "state", "state")):
        t = readable(p)
        if t is None:
            print(f"restart: skipping unreadable {p.name}", file=sys.stderr)
            continue
        if time == 0:
            return len(t) - 1, float(t[-1]), p
        if t[0] <= time + 1e-12 or p == set_files(run / "state", "state")[0]:
            i = int(np.argmin(np.abs(t - time)))
            return i, float(t[i]), p
    raise FileNotFoundError(f"no readable state files in {run / 'state'}")


def has_checkpoint(run: Path) -> bool:
    """Whether a run has at least one readable state set of its own.

    Parameters
    ----------
    run : pathlib.Path
        Run directory.

    Returns
    -------
    bool
    """
    state = run / "state"
    return state.is_dir() and any(readable(p) is not None for p in set_files(state, "state"))


def keep_and_prune(run: Path, keep: int, prune: bool) -> list[Path]:
    """The newest readable state sets, optionally deleting all older ones.

    Nothing is deleted unless ``keep`` readable sets exist, and the set that
    ``restart/restart.h5`` points to is never deleted.

    Parameters
    ----------
    run : pathlib.Path
        Run directory.
    keep : int
        Number of sets to keep.
    prune : bool
        Delete sets older than the oldest kept one (folder and ``.h5``).

    Returns
    -------
    list of pathlib.Path
        Kept set files, newest first.
    """
    sets = set_files(run / "state", "state")
    good = [p for p in reversed(sets) if readable(p) is not None][:keep]
    if len(good) < keep:
        print(f"only {len(good)} readable state sets in {run}; nothing pruned", file=sys.stderr)
        return good
    link = run / "restart" / "restart.h5"
    protect = {link.resolve()} if link.exists() else set()
    if prune:
        oldest_kept = min(int(p.stem.rsplit("_s", 1)[1]) for p in good)
        for p in sets:
            if int(p.stem.rsplit("_s", 1)[1]) < oldest_kept and p.resolve() not in protect:
                shutil.rmtree(p.with_suffix(""), ignore_errors=True)
                p.unlink(missing_ok=True)
                print(f"pruned {p.name}", file=sys.stderr)
    return good


def get_parser() -> argparse.ArgumentParser:
    """Command-line options (see the module docstring).

    Returns
    -------
    argparse.ArgumentParser
    """
    a = argparse.ArgumentParser(prog="python -m rbc.restart",
                                description="Choose the restart checkpoint; keep or prune old ones.")
    a.add_argument("run", type=Path, help="run directory")
    a.add_argument("--time", type=float, default=0.0, help="restart time (0 = newest checkpoint)")
    a.add_argument("--keep", type=int, help="print the N newest readable sets (folder, then .h5)")
    a.add_argument("--prune", action="store_true", help="with --keep: delete older sets")
    a.add_argument("--from", dest="src", type=Path, help="take the checkpoint from this run")
    a.add_argument("--shape", type=int, nargs=3, metavar=("NX", "NY", "NZ"),
                   help="target resolution; the checkpoint is regridded if it differs")
    return a


def main(argv=None) -> None:
    """Choose the checkpoint and link it (or regrid it) to ``restart/restart.h5``.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments (default: ``sys.argv[1:]``).
    """
    a = get_parser().parse_args(argv)
    run = a.run.resolve()
    if a.keep:
        for p in keep_and_prune(run, a.keep, a.prune):
            print(p.with_suffix(""))
            print(p)
        return
    src = run
    if a.src and a.src.resolve() != run:
        if has_checkpoint(run):
            print(f"restart: {run.name} has its own checkpoints; ignoring --from {a.src}", file=sys.stderr)
        else:
            src = a.src.resolve()
    i, t, p = choose(src, a.time)
    old = None
    if a.shape:
        from rbc import regrid
        old = regrid.coeff_shape(p)
    change = old is not None and old != tuple(a.shape)
    if change and src == run:
        sys.exit(f"restart: checkpoint {p.name} has shape {old} but RES_HOR/RES_Z give "
                 f"{tuple(a.shape)}. Change resolution in a new run folder with "
                 f"RESTART_FROM={run}")
    link = run / "restart" / "restart.h5"
    link.parent.mkdir(exist_ok=True)
    if link.is_symlink() or link.exists():
        link.unlink()
    if change:
        print(f"restart: regridding {old} -> {tuple(a.shape)} from {p}", file=sys.stderr)
        regrid.regrid(p, link, tuple(a.shape), i)
        print(0, repr(t), link)
        return
    link.symlink_to(p.resolve())
    print(i, repr(t), p)

if __name__ == "__main__":
    main()