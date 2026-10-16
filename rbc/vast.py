"""Dedalus v3 file-handler fixes for NFS-like filesystems (VAST on Trillium).
Harmless on Lustre. Replaces the ``fixed_get_file.py`` wrapper. Usage::

    import rbc.vast; rbc.vast.patch()          # before building file handlers
    write, dt = rbc.vast.load_state(solver, restart_path, index)

1. get_file (the hang). Dedalus asks every rank, at every write, whether the
   current set file exists; a new set is created inside MPI barriers. On NFS a
   rank on another node can hold a cached "does not exist" answer, enters the
   create barrier alone, and the job deadlocks. Mike's fix made rank 0 decide.
   Here no rank asks the filesystem: the handler remembers which set number it
   has created, which is identical on every rank by construction. No stat()
   storm (N ranks x every write), no broadcast.

2. Set folders. Each rank creates its process file inside a folder that rank 0
   just made; other nodes may still have the folder cached as missing. Dedalus'
   workaround touches and deletes one temp file per rank (2N metadata ops in one
   directory). The cache is per node, so one mkdir(exist_ok) per node refreshes
   it: ~N/192 operations.

3. Restart load (the slow start). The joint state file is a virtual dataset
   whose time axis is unlimited in every mapping, so HDF5 opens every source
   file to resolve its size, even to read one rank's slab. N ranks x N process
   files = N^2 opens (5184^2 = 2.7e7). Each rank instead reads its own process
   file (1 open), after checking it holds exactly the slab the current
   decomposition needs; any mismatch on any rank falls back to Dedalus' loader.
   Every rank must still set every field's layout and scales, even with no local
   data, or the first transposes deadlock (see :func:`load_state`).
"""
from __future__ import annotations

import os
from pathlib import Path

import h5py
import numpy as np

_NODE = {}


def _node_comm(comm):
    """Communicator of the ranks on this node (cached per parent communicator)."""
    key = id(comm)
    if key not in _NODE:
        from mpi4py import MPI
        _NODE[key] = comm.Split_type(MPI.COMM_TYPE_SHARED)
    return _NODE[key]


def patch() -> None:
    """Replace Dedalus' set-file creation with the node-aware, stat-free version.

    Patches ``H5FileHandlerBase.get_file`` and ``H5VirtualFileHandler.create_current_file``
    in ``dedalus.core.evaluator``. Call before any file handler is created.
    """
    from dedalus.core import evaluator as ev
    from dedalus.tools.parallel import Sync

    def get_file(self, **kw):
        """Open the current set file, creating it once per set without asking the
        filesystem whether it exists."""
        if getattr(self, "_rbc_created_set", None) != self.set_num:
            self.create_current_file()
            self._rbc_created_set = self.set_num
        return self.open_file(**kw)

    def create_virtual(self):
        """Create the set folder once per node, then this rank's process file."""
        ev.H5FileHandlerBase.create_current_file(self)
        with Sync(self.comm):
            if self.comm.rank == 0:
                self.current_path.mkdir(exist_ok=True)
        node = _node_comm(self.comm)
        if node.rank == 0:
            os.makedirs(self.current_path, exist_ok=True)   # refreshes this node's cache
        node.Barrier()
        if not self.empty:
            with h5py.File(str(self.current_process_file), "w-") as f:
                self.setup_process_file(f)

    ev.H5FileHandlerBase.get_file = get_file
    ev.H5VirtualFileHandler.create_current_file = create_virtual


def _expected(field, layout, scales):
    """Global shape, local start and local shape that this rank should hold."""
    dim = field.dist.dim
    els = layout.local_elements(field.domain, scales)
    gshape = layout.global_shape(field.domain, scales)
    start = tuple(int(e[0]) if e.size else int(g) for e, g in zip(els, gshape))
    shape = tuple(int(e.size) for e in els)
    return tuple(int(g) for g in gshape), start[-dim:], shape[-dim:]


def load_state(solver, path, index=-1, logger=None):
    """Load a checkpoint reading only this rank's process file.

    Drop-in replacement for ``solver.load_state(path, index)``. Falls back to Dedalus'
    loader (on every rank) if the file is not a virtual set with matching tasks or if
    any rank's process file does not hold exactly the slab the current
    decomposition needs.

    Parameters
    ----------
    solver : dedalus.core.solvers.InitialValueSolver
        The solver whose state is loaded.
    path : str or pathlib.Path
        ``restart/restart.h5`` (a link to a ``state_s<N>.h5``).
    index : int, optional
        Write within the file.
    logger : logging.Logger, optional
        For progress and fallback messages.

    Returns
    -------
    write : int
        Write number of the loaded state.
    dt : float
        Timestep at that write (the restart's first timestep with CFL).
    """
    comm = solver.dist.comm_cart          # process files are named by this rank
    path = Path(path).resolve()
    names = [f.name for f in solver.state]
    meta = None
    if comm.rank == 0:
        with h5py.File(path, "r") as f:
            sc = f["scales"]
            meta = dict({k: sc[k][index] for k in ("write_number", "timestep", "iteration", "sim_time")},
                        tasks={n: (bool(np.all(f["tasks"][n].attrs["grid_space"])),
                                   tuple(f["tasks"][n].attrs["scales"]))
                               for n in names if n in f["tasks"]},
                        virtual=all(f["tasks"][n].is_virtual for n in names if n in f["tasks"]))
    meta = comm.bcast(meta, root=0)

    plan, ok = [], meta["virtual"] and set(meta["tasks"]) == set(names)
    pfile = path.with_suffix("") / f"{path.stem}_p{comm.rank}.h5"
    if ok:
        try:
            h = h5py.File(pfile, "r") if pfile.exists() else None
            for field in solver.state:
                grid, scales = meta["tasks"][field.name]
                layout = solver.dist.grid_layout if grid else solver.dist.coeff_layout
                scales = solver.dist.remedy_scales(scales if grid else 1)
                gshape, start, shape = _expected(field, layout, scales)
                if not np.prod(shape):
                    plan.append((field, layout, scales, grid, None))   # still set layout/scales below
                    continue
                d = h["tasks"][field.name]
                a = d.attrs
                if (tuple(a["global_shape"][-len(gshape):]) != gshape
                        or tuple(a["local_start"][-len(start):]) != start
                        or tuple(a["local_shape"][-len(shape):]) != shape):
                    ok = False
                    break
                plan.append((field, layout, scales, grid, d[index]))
            if h is not None:
                h.close()
        except Exception as e:
            if logger:
                logger.warning(f"vast.load_state rank {comm.rank}: {e!r}")
            ok = False
    from mpi4py import MPI
    if not comm.allreduce(ok, op=MPI.LAND):
        if logger and comm.rank == 0:
            logger.warning("vast.load_state: decomposition changed; using Dedalus load_state")
        return solver.load_state(path, index=index)

    # Every rank must leave every field in the same layout and scales, including ranks
    # holding none of it (e.g. tau fields in grid space): which fields get transposed
    # together decides the collective FFTW transpose plans, so a mismatch deadlocks.
    for field, layout, scales, grid, data in plan:
        if grid:
            field.preset_scales(scales)
        field.preset_layout(layout)
        if data is not None:
            np.copyto(field.data, data)
        if grid:
            field.change_scales(field.domain.dealias)
    solver.iteration = solver.initial_iteration = meta["iteration"]
    solver.sim_time = solver.initial_sim_time = meta["sim_time"]
    if logger and comm.rank == 0:
        logger.info(f"Loaded state from process files: write {meta['write_number']}, "
                    f"t = {meta['sim_time']}, dt = {meta['timestep']}")
    return meta["write_number"], meta["timestep"]


def verify_load(solver, path, index=-1, logger=None) -> float:
    """Load a checkpoint with both loaders and compare.

    Meant to be run once at small scale (``--verify_load``) after a change to
    :func:`load_state`.

    Parameters
    ----------
    solver : dedalus.core.solvers.InitialValueSolver
        The solver.
    path : str or pathlib.Path
        Checkpoint.
    index : int, optional
        Write within the file.
    logger : logging.Logger, optional
        Reports the result on rank 0.

    Returns
    -------
    float
        Largest absolute difference over all coefficients and ranks (0 when correct).
    """
    load_state(solver, path, index)
    fast = {f.name: np.copy(f["c"]) for f in solver.state}
    solver.load_state(path, index=index)
    diff = max(float(np.max(np.abs(fast[f.name] - f["c"]), initial=0.0)) for f in solver.state)
    from mpi4py import MPI
    diff = solver.dist.comm_cart.allreduce(diff, op=MPI.MAX)
    if logger and solver.dist.comm_cart.rank == 0:
        logger.info(f"verify_load: max |fast - dedalus| over all coefficients = {diff:.3e}")
    return diff
