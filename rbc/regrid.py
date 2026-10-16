"""Change the resolution of a coefficient-space checkpoint (exact spectral interpolation).

Dedalus coefficients do not depend on N: RealFourier stores unit-amplitude
(cos k, -sin k) pairs at element 2k, 2k+1 (Nyquist dropped) and ChebyshevT uses
the orthonormal basis, coefficient n for degree n (see rbc.spectra.core). So moving
to a new resolution is zero-padding (refine) or truncation (coarsen) of every axis;
the represented field is unchanged when refining, low-pass filtered when coarsening.
Axes of size 1 (constant directions of tau fields) are left alone.

Usage::

    python -m rbc.regrid SRC.h5 OUT.h5 --nx 256 --nz 128 [--ny 256] [--index -1]

OUT.h5 holds one write as plain (non-virtual) datasets, readable by Dedalus'
load_state on any process mesh. Grid-space checkpoints need no regrid: Dedalus
interpolates those itself when the shape changes.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np

SCALES = ("sim_time", "timestep", "wall_time", "iteration", "write_number")


def is_coeff(path: Path) -> bool:
    """Whether every task of a checkpoint is stored in coefficient space.

    Parameters
    ----------
    path : pathlib.Path
        Checkpoint (``state_s<N>.h5``).

    Returns
    -------
    bool
    """
    with h5py.File(path, "r") as f:
        return not any(np.all(f["tasks"][n].attrs["grid_space"]) for n in f["tasks"])


def coeff_shape(path: Path) -> tuple[int, int, int] | None:
    """Coefficient shape of a checkpoint.

    Parameters
    ----------
    path : pathlib.Path
        Checkpoint file.

    Returns
    -------
    tuple of int or None
        (Nx, Ny, Nz) of its largest task, or None for a grid-space checkpoint.
    """
    with h5py.File(path, "r") as f:
        shapes = [f["tasks"][n].shape[-3:] for n in f["tasks"]]
    return tuple(int(v) for v in max(shapes, key=np.prod)) if is_coeff(path) else None


def regrid(src: Path, out: Path, shape: tuple[int, int, int], index: int = -1) -> None:
    """Write one write of a coefficient checkpoint at a new resolution.

    Every axis is zero-padded (refining) or truncated (coarsening). The file is
    written to ``<out>.part`` and renamed at the end, so a failure leaves no
    half-written restart file.

    Parameters
    ----------
    src : pathlib.Path
        Coefficient-space checkpoint.
    out : pathlib.Path
        New file (plain datasets, readable by Dedalus on any process mesh).
    shape : tuple of int
        Target (Nx, Ny, Nz); Nx and Ny must be even.
    index : int, optional
        Write within ``src``.

    Raises
    ------
    ValueError
        For odd Fourier sizes or a grid-space task.
    """
    nx, ny, nz = shape
    if nx % 2 or ny % 2:
        raise ValueError("Fourier sizes must be even")
    tmp = out.with_suffix(".h5.part")
    with h5py.File(src, "r") as f, h5py.File(tmp, "w") as g:
        g.attrs.update(dict(f.attrs))
        g.attrs["writes"] = 1
        g.attrs["regridded_from"] = str(Path(src).resolve())
        for k in SCALES:
            if k in f["scales"]:
                g.create_dataset(f"scales/{k}", data=f["scales"][k][[index]])
        g.create_group("tasks")
        for name in f["tasks"]:
            d = f["tasks"][name]
            if np.all(d.attrs["grid_space"]):
                raise ValueError(f"{name} is in grid space; Dedalus can load it at a new resolution directly")
            old = d.shape[-3:]
            new = tuple(o if o == 1 else t for o, t in zip(old, shape))
            lead = d.shape[1:-3]                       # tensor components
            o = g.create_dataset(f"tasks/{name}", shape=(1,) + lead + new, dtype=d.dtype, fillvalue=0)
            for a, v in d.attrs.items():
                o.attrs[a] = v
            m = tuple(slice(0, min(a, b)) for a, b in zip(old, new))
            for comp in np.ndindex(*lead) if lead else [()]:
                o[(0,) + comp + m] = d[(index,) + comp + m]   # one component at a time (memory)
    tmp.rename(out)


def get_parser() -> argparse.ArgumentParser:
    """Command-line options.

    Returns
    -------
    argparse.ArgumentParser
    """
    p = argparse.ArgumentParser(prog="python -m rbc.regrid", description=__doc__.split("\n")[0])
    p.add_argument("src", type=Path, help="coefficient-space checkpoint")
    p.add_argument("out", type=Path, help="file to write")
    p.add_argument("--nx", type=int, required=True, help="Fourier modes in x")
    p.add_argument("--ny", type=int, help="Fourier modes in y (default: nx)")
    p.add_argument("--nz", type=int, required=True, help="Chebyshev modes in z")
    p.add_argument("--index", type=int, default=-1, help="write within the source file")
    return p


def main(argv=None) -> None:
    """Regrid one checkpoint from the command line.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments (default: ``sys.argv[1:]``).
    """
    a = get_parser().parse_args(argv)
    regrid(a.src, a.out, (a.nx, a.ny or a.nx, a.nz), a.index)
    print(a.out)


if __name__ == "__main__":
    main()
