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
   file to resolve its size, even to read one rank's slab. N ranks x M process
   files = N*M opens (5184^2 = 2.7e7). Instead rank 0 reads the virtual mapping
   (which process file holds which block; no process file is opened) and each
   rank opens only the process files that overlap its own slab: one when the
   process mesh is unchanged, a few when it changed. Dedalus' loader is used
   only for plain files (regridded restarts), a change of resolution, or a gap
   in the mapping. Every rank must still set every field's layout and scales,
   even with no local data, or the first transposes deadlock (see :func:`load_state`).
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


def _regions(dset, base):
    """Process files of a virtual task and the block of the global array each holds.

    Reads only the dataset's virtual mapping, so no process file is opened.

    Parameters
    ----------
    dset : h5py.Dataset
        Virtual task of a joint set file.
    base : pathlib.Path
        Folder of the joint file; relative source names resolve against it.

    Returns
    -------
    gshape : tuple of int
        Global shape without the write axis (tensor components, then space).
    regions : list of tuple
        ``(file, dataset, start, stop)`` for each process file, with ``start`` and
        ``stop`` over the axes of ``gshape``.

    Raises
    ------
    ValueError
        If a mapping is not one contiguous block, or does not span every tensor component.
    """
    gshape, regions = None, []
    for s in dset.virtual_sources():
        gshape = tuple(int(n) for n in s.vspace.shape[1:])
        start, stride, count, block = (v[1:] for v in s.vspace.get_regular_hyperslab())
        if any(c > 1 and st != b for st, c, b in zip(stride, count, block)):
            raise ValueError(f"{dset.name}: mapping of {s.file_name} is not a contiguous block")
        lo = tuple(int(a) for a in start)
        hi = tuple(int(a + (c - 1) * st + b) for a, st, c, b in zip(start, stride, count, block))
        regions.append((str(base / s.file_name), s.dset_name, lo, hi))
    return gshape, regions


def _read_slab(regions, ncomp, start, shape, index, files):
    """Assemble one rank's slab of a task from the process files that overlap it.

    Parameters
    ----------
    regions : list of tuple
        From :func:`_regions`.
    ncomp : int
        Number of leading tensor-component axes (read whole).
    start, shape : tuple of int
        The rank's slab over the spatial axes.
    index : int
        Write within the set (non-negative).
    files : dict
        Open process files by path; files opened here are added.

    Returns
    -------
    numpy.ndarray or None
        The slab (components, then space), or None if the mapped blocks leave part
        of it uncovered.

    Raises
    ------
    ValueError
        If a process file's dataset does not have the shape its mapping says.
    """
    lo = np.array(start)
    hi = lo + np.array(shape)
    buf, covered = None, 0
    for fname, dname, rlo, rhi in regions:
        slo, shi = np.array(rlo[ncomp:]), np.array(rhi[ncomp:])
        a, b = np.maximum(lo, slo), np.minimum(hi, shi)
        if np.any(a >= b):
            continue
        if fname not in files:
            files[fname] = h5py.File(fname, "r")
        d = files[fname][dname]
        if tuple(d.shape[1:]) != tuple(h - l for l, h in zip(rlo, rhi)):
            raise ValueError(f"{fname}: {dname} has shape {d.shape[1:]}, mapping says {rlo}..{rhi}")
        if buf is None:
            buf = np.empty(tuple(d.shape[1:ncomp + 1]) + tuple(shape), d.dtype)
        src = (index,) + (slice(None),) * ncomp + tuple(slice(x - y, z - y) for x, y, z in zip(a, slo, b))
        dst = (Ellipsis,) + tuple(slice(x - y, z - y) for x, y, z in zip(a, lo, b))
        buf[dst] = d[src]
        covered += int(np.prod(b - a))
    return buf if covered == int(np.prod(shape)) else None


_FALLBACK = {1: "not a virtual set with the solver's tasks (e.g. a regridded restart, "
                "which Dedalus reads with one open per rank)",
             2: "resolution differs from the checkpoint",
             3: "the process files do not cover every rank's slab",
             4: "error while reading process files (see the warnings above)"}


def load_state(solver, path, index=-1, logger=None):
    """Load a checkpoint, each rank reading only the process files that overlap its slab.

    Drop-in replacement for ``solver.load_state(path, index)``. Works whether or not
    the process mesh changed since the checkpoint was written. Falls back to Dedalus'
    loader (on every rank) for a plain file, a change of resolution, or a gap in the
    mapping.

    Parameters
    ----------
    solver : dedalus.core.solvers.InitialValueSolver
        The solver whose state is loaded.
    path : str or pathlib.Path
        ``restart/restart.h5`` (a link to a ``state_s<N>.h5``, or a regridded file).
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
    from mpi4py import MPI

    comm = solver.dist.comm_cart
    path = Path(path).resolve()
    names = [f.name for f in solver.state]
    meta = None
    if comm.rank == 0:
        with h5py.File(path, "r") as f:
            sc = f["scales"]
            idx = index % len(sc["sim_time"])
            meta = {k: sc[k][idx] for k in ("write_number", "timestep", "iteration", "sim_time")}
            meta["index"] = idx
            tasks = f["tasks"]
            meta["virtual"] = (set(names) <= set(tasks)
                               and all(tasks[n].is_virtual for n in names))
            if meta["virtual"]:
                meta["tasks"] = {n: (bool(np.all(tasks[n].attrs["grid_space"])),
                                     tuple(tasks[n].attrs["scales"]),
                                     *_regions(tasks[n], path.parent))
                                 for n in names}
    meta = comm.bcast(meta, root=0)

    plan, files, code = [], {}, 0 if meta["virtual"] else 1
    if code == 0:
        try:
            for field in solver.state:
                grid, scales, gshape_file, regions = meta["tasks"][field.name]
                layout = solver.dist.grid_layout if grid else solver.dist.coeff_layout
                scales = solver.dist.remedy_scales(scales if grid else 1)
                gshape, start, shape = _expected(field, layout, scales)
                ncomp = len(gshape_file) - len(gshape)
                if gshape_file[ncomp:] != gshape:
                    code = 2
                    break
                data = None
                if np.prod(shape):
                    data = _read_slab(regions, ncomp, start, shape, meta["index"], files)
                    if data is None:
                        code = 3
                        break
                plan.append((field, layout, scales, grid, data))
        except Exception as e:
            if logger:
                logger.warning(f"vast.load_state rank {comm.rank}: {e!r}")
            code = 4
        finally:
            for h in files.values():
                h.close()
    code = comm.allreduce(code, op=MPI.MAX)
    if code:
        if logger and comm.rank == 0:
            logger.warning(f"vast.load_state: using Dedalus load_state ({_FALLBACK[code]})")
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
    most = comm.allreduce(len(files), op=MPI.MAX)
    if logger and comm.rank == 0:
        logger.info(f"Loaded state from process files (at most {most} per rank): write "
                    f"{meta['write_number']}, t = {meta['sim_time']}, dt = {meta['timestep']}")
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
    load_state(solver, path, index, logger)
    fast = {f.name: np.copy(f["c"]) for f in solver.state}
    solver.load_state(path, index=index)
    diff = max(float(np.max(np.abs(fast[f.name] - f["c"]), initial=0.0)) for f in solver.state)
    from mpi4py import MPI
    diff = solver.dist.comm_cart.allreduce(diff, op=MPI.MAX)
    if logger and solver.dist.comm_cart.rank == 0:
        logger.info(f"verify_load: max |fast - dedalus| over all coefficients = {diff:.3e}")
    return diff