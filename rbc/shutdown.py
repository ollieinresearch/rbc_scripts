"""Collective, deadlock-free early stop for Dedalus time loops.

Signals reach MPI ranks at slightly different times; a rank that leaves the loop
alone hangs the others in the next transpose. Here a signal (or passing the
deadline) only sets a local flag. Every `cadence` iterations all ranks combine
flags with one Allreduce(MAX), so they all leave on the same iteration and can
write a final checkpoint collectively.

In ``sim/rayleigh_benard.py``::

    from rbc.shutdown import GracefulStop
    stop = GracefulStop(MPI.COMM_WORLD, logger=logger)   # deadline from $RBC_DEADLINE
    while solver.proceed:
        dt = CFL.compute_timestep()
        solver.step(dt)
        ...
        if stop.check(solver.iteration):
            break
    if stop.safe_to_write:                               # normal end, deadline or SIGUSR1
        solver.evaluate_handlers([checkpoints], dt=dt)   # resume exactly here next job
    stop.record(basepath)                                # writes stop_reason.txt (rank 0)

Slurm: `#SBATCH --signal=USR1@<s>` (no "B:") signals every task directly;
`scancel --signal=USR1 <jobid>` stops a run cleanly by hand.
RBC_GRACEFUL=0 disables all of this (no handlers, no deadline, no collective check).
"""
from __future__ import annotations

import os
import signal
import time
from pathlib import Path

_REASONS = {0: None, 1: "deadline", int(signal.SIGUSR1): "SIGUSR1",
            int(signal.SIGTERM): "SIGTERM", int(signal.SIGINT): "SIGINT"}


class GracefulStop:
    """Collective early stop for a Dedalus time loop.

    Parameters
    ----------
    comm : mpi4py.MPI.Comm
        Communicator of the simulation.
    deadline : float, optional
        Unix time at which to stop; default ``$RBC_DEADLINE`` (set by
        ``set_deadline`` in ``shell/workflow.sh``), or none.
    cadence : int, optional
        Iterations between the collective checks.
    signals : tuple of signal.Signals, optional
        Signals that request a stop.
    logger : logging.Logger, optional
        Where to report arming and stopping (rank 0 only).

    Attributes
    ----------
    enabled : bool
        False when ``RBC_GRACEFUL=0``: no handlers, no deadline, no collective check.
    reason : str or None
        Why the loop stopped (``"deadline"``, ``"SIGUSR1"``, ...), None while running.
    """

    def __init__(self, comm, deadline: float | None = None, cadence: int = 10,
                 signals=(signal.SIGUSR1, signal.SIGTERM, signal.SIGINT), logger=None):
        from mpi4py import MPI

        self.comm, self.cadence, self.logger, self._max = comm, max(1, cadence), logger, MPI.MAX
        self.enabled = os.environ.get("RBC_GRACEFUL", "1") != "0"   # RBC_GRACEFUL=0: plain behaviour
        self._local = 0
        self.reason = None
        if not self.enabled:
            self.deadline = None
            if self.logger and comm.rank == 0:
                self.logger.info("GracefulStop disabled (RBC_GRACEFUL=0)")
            return
        env = os.environ.get("RBC_DEADLINE")
        self.deadline = deadline if deadline is not None else (float(env) if env else None)
        for s in signals:
            signal.signal(s, self._on_signal)
        if self.logger and comm.rank == 0 and self.deadline:
            self.logger.info(f"GracefulStop armed: if still running in {self.deadline - time.time():.0f} s, "
                             "stop and write a final checkpoint")

    def _on_signal(self, signum, frame):
        """Signal handler: only set a flag (MPI and HDF5 calls are not safe here)."""
        self._local = max(self._local, int(signum))

    def check(self, iteration: int) -> bool:
        """Whether every rank should leave the loop now.

        Collective: every rank must call this every iteration (it only communicates
        every ``cadence`` iterations, with one Allreduce).

        Parameters
        ----------
        iteration : int
            Current solver iteration.

        Returns
        -------
        bool
            True on the same iteration on every rank once any rank has seen a signal
            or passed the deadline.
        """
        if not self.enabled or iteration % self.cadence:
            return False
        if self.deadline is not None and time.time() >= self.deadline:
            self._local = max(self._local, 1)
        flag = self.comm.allreduce(self._local, op=self._max)
        if flag:
            self.reason = _REASONS.get(flag, f"signal {flag}")
            if self.logger and self.comm.rank == 0:
                self.logger.warning(f"Stopping at iteration {iteration}: {self.reason}")
        return bool(flag)

    @property
    def safe_to_write(self) -> bool:
        """Whether a final checkpoint can be written.

        False after SIGTERM or SIGINT: SIGTERM means Slurm's hard kill is about 30 s
        away (KillWait), and a large checkpoint cut off midway would leave a corrupt
        newest state file.
        """
        return self.reason not in ("SIGTERM", "SIGINT")

    def record(self, basepath) -> None:
        """Write the stop reason to ``stop_reason.txt`` (rank 0).

        Parameters
        ----------
        basepath : str or pathlib.Path
            Run directory.
        """
        if self.comm.rank == 0:
            Path(basepath, "stop_reason.txt").write_text(f"{self.reason or 'finished'}\n")
