"""Stop a run cleanly before any node runs out of memory.

Linux nodes without swap do not slow down as memory fills; the job runs at full speed
until the out-of-memory killer ends it, losing everything since the last checkpoint.
So the guard watches the memory directly. About every ``interval`` seconds of wall time
the first rank of each node reads ``MemTotal - MemAvailable`` (what the killer sees)
and all ranks take the maximum over nodes with one Allreduce: the fullest node decides,
since one node running out ends the whole job.

The checks happen at iterations every rank agrees on, without signals: the first is
100 iterations after the start, and at each check rank 0 turns the measured time per
iteration into the number of iterations until the next one, which is broadcast with
the result.

* Above ``stop_frac`` of a node's memory, :meth:`MemoryGuard.check` returns True; the
  caller leaves the loop and writes its final checkpoint (``stop_reason`` "memory").
  The next job starts from that checkpoint in fresh processes, with the memory that
  built up during the run released.
* Above ``warn_frac`` during the first ``early_s`` seconds, a warning is logged once:
  the configuration starts close to the limit, so its jobs will be cut short. More
  nodes would help.

Environment overrides: ``RBC_MEM_STOP`` (default 0.95; 0 disables the guard),
``RBC_MEM_WARN`` (0.85), ``RBC_MEM_EARLY`` (1200 s), ``RBC_MEM_INTERVAL`` (300 s).
"""
from __future__ import annotations

import os
import socket
import time
from pathlib import Path


def used_fraction() -> tuple[float, float]:
    """Memory in use on this node.

    Returns
    -------
    fraction : float
        ``(MemTotal - MemAvailable) / MemTotal`` (0 if ``/proc/meminfo`` is unreadable).
    used : float
        GiB in use.
    """
    try:
        info = dict(line.split(":", 1) for line in Path("/proc/meminfo").read_text().splitlines())
        total = float(info["MemTotal"].split()[0])
        avail = float(info["MemAvailable"].split()[0])
        return (total - avail) / total, (total - avail) / 2**20
    except (OSError, KeyError, ValueError):
        return 0.0, 0.0


class MemoryGuard:
    """Collective memory check; every rank must call :meth:`check` every iteration.

    Parameters
    ----------
    comm : mpi4py.MPI.Comm
        All ranks.
    first : int, optional
        Iterations from the start to the first check.
    logger : logging.Logger, optional
        For warnings and :meth:`status`.
    """

    def __init__(self, comm, first: int = 100, logger=None):
        from mpi4py import MPI

        self.comm, self.first, self.logger, self._max = comm, max(1, first), logger, MPI.MAX
        self.interval = float(os.environ.get("RBC_MEM_INTERVAL", "300") or 300)
        self.next = None
        self.mark = None
        self.stop_frac = float(os.environ.get("RBC_MEM_STOP", "0.95") or 0)
        self.warn_frac = float(os.environ.get("RBC_MEM_WARN", "0.85") or 0)
        self.early_s = float(os.environ.get("RBC_MEM_EARLY", "1200") or 0)
        node = comm.Split_type(MPI.COMM_TYPE_SHARED)
        self.leader = node.rank == 0
        node.Free()
        self.t0 = time.time()
        self.warned = False
        self.last = (0.0, 0.0, "")
        self.peak = 0.0

    def _measure(self):
        frac, gib = used_fraction() if self.leader else (0.0, 0.0)
        frac, gib, host = self.comm.allreduce((frac, gib, socket.gethostname() if self.leader else ""),
                                              op=self._max)
        return frac, gib, host

    def check(self, iteration: int) -> bool:
        """Whether every rank should stop now because a node is nearly full.

        Parameters
        ----------
        iteration : int
            Solver iteration (identical on every rank).

        Returns
        -------
        bool
            The same on every rank.
        """
        if not self.stop_frac:
            return False
        if self.next is None:
            self.next, self.mark = iteration + self.first, (iteration, time.time())
        if iteration < self.next:
            return False
        frac, gib, host = self.last = self._measure()
        now = time.time()
        per_iter = (now - self.mark[1]) / max(iteration - self.mark[0], 1)
        step = self.comm.bcast(int(min(max(self.interval / max(per_iter, 1e-6), 10), 10**6)), root=0)
        self.next, self.mark = iteration + step, (iteration, now)
        self.peak = max(self.peak, frac)
        if (not self.warned and self.warn_frac and frac > self.warn_frac
                and time.time() - self.t0 < self.early_s):
            self.warned = True
            if self.logger and self.comm.rank == 0:
                self.logger.warning(f"memory: {host} at {frac:.0%} ({gib:.0f} GiB) within "
                                    f"{self.early_s / 60:.0f} min of the start; jobs will be cut short "
                                    f"by the memory guard. More nodes would help.")
        if frac > self.stop_frac:
            if self.logger and self.comm.rank == 0:
                self.logger.warning(f"memory: {host} at {frac:.0%} ({gib:.0f} GiB) > {self.stop_frac:.0%}; "
                                    f"stopping and writing a checkpoint")
            return True
        return False

    def status(self) -> str:
        """Last measurement, for the periodic log line.

        Returns
        -------
        str
        """
        frac, gib, host = self.last
        return f"mem max {frac:.0%} ({gib:.0f} GiB, {host}), peak {self.peak:.0%}" if host else "mem n/a"
