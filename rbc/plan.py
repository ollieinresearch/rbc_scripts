"""Choose node counts, ranks per node and process meshes that divide the problem evenly.

Dedalus splits each distributed axis into blocks of ``ceil(G / P)`` groups, so when the
mesh dimension ``P`` does not divide the group count ``G`` the trailing ranks get a short
block or nothing, and every step waits for the fullest rank. Ranks with no Fourier groups
also hold no subproblem matrices, which is why they show up as a node using less memory.

Group counts per mesh dimension (``RES_HOR = gamma * RES_Z``, 3/2 dealiasing):

* ``RES_HOR / 2``  RealFourier cos/sin pairs (coefficient space; sets the solve)
* ``3 RES_HOR / 2`` Fourier grid points
* ``RES_Z``        Chebyshev coefficients
* ``3 RES_Z / 2``  Chebyshev grid points

Both mesh dimensions are checked against all four, so the result does not depend on
which axis each dimension lands on in a given layout. For gamma = 2 a mesh dimension
divides all four exactly when it divides ``RES_Z / 2``.

The efficiency of a mesh is ``e(P0) * e(P1)`` with ``e(P) = min over G of G / (P ceil(G/P))``:
the ideal work per rank over the work of the fullest rank (1 = balanced).

Usage::

    python -m rbc.plan --res 256                      # 1-12 nodes, C = 1, 2, 3
    python -m rbc.plan --res 240 256 288 --nodes 4 10 --all
    python -m rbc.plan --res 288 --mem-ref 256 3300   # scale memory from a measured run
    python -m rbc.plan --res 240 256 288 --nodes 5 9 --cpt 1 --all --top 2 --emit > bench.conf
        configurations for ``trill_bench``: one ``RES_Z NODES C MESHX MESHY`` line each

Ranks per node follow ``trill_run_sim``'s packing: ``8 // C`` ranks on each 8-core L3
cache, so C = 3 gives 2 ranks per cache (48 per node, 2 cores idle per cache), not 64.
``--l3 0`` gives ``cores // C`` instead (srun on the Lustre clusters).
"""
from __future__ import annotations

import argparse
import math


def groups(res_z: int, gamma: float = 2.0) -> tuple[int, ...]:
    """Group counts that a mesh dimension may have to split.

    Parameters
    ----------
    res_z : int
        Chebyshev modes (``RES_Z``).
    gamma : float, optional
        Aspect ratio; ``RES_HOR = gamma * RES_Z``.

    Returns
    -------
    tuple of int
        Fourier pairs, Fourier grid points, Chebyshev modes, Chebyshev grid points.
    """
    nh = round(gamma * res_z)
    return nh // 2, 3 * nh // 2, res_z, 3 * res_z // 2


def efficiency(p: int, sizes: tuple[int, ...]) -> float:
    """Load balance of splitting each size into ``p`` Dedalus blocks.

    Parameters
    ----------
    p : int
        Mesh dimension.
    sizes : tuple of int
        Group counts (from :func:`groups`).

    Returns
    -------
    float
        Worst ``G / (p * ceil(G / p))`` over the sizes (1 when ``p`` divides all).
    """
    return min(g / (p * math.ceil(g / p)) for g in sizes)


def meshes(n: int, sizes: tuple[int, ...]):
    """Every factorisation ``n = p0 * p1`` with its efficiency.

    Parameters
    ----------
    n : int
        Number of ranks.
    sizes : tuple of int
        Group counts.

    Yields
    ------
    tuple
        ``(p0, p1, efficiency)``.
    """
    for p0 in range(1, n + 1):
        if n % p0 == 0:
            p1 = n // p0
            yield p0, p1, efficiency(p0, sizes) * efficiency(p1, sizes)


def ranks_per_node(c: int, cores: int = 192, l3: int = 8) -> int:
    """Ranks per node for ``cpus-per-task = c``.

    Parameters
    ----------
    c : int
        Cores per rank (``C``).
    cores : int, optional
        Cores per node.
    l3 : int, optional
        Cores per L3 cache when ranks are packed per cache (``trill_run_sim``); 0 for
        plain packing (``cores // c``).

    Returns
    -------
    int
    """
    return (cores // l3) * (l3 // c) if l3 else cores // c


def plan(res_list, nodes=(1, 12), cpt=(1, 2, 3), cores=192, gamma=2.0,
         mem_ref=None, show_all=False, l3=8, top=1):
    """Rows of the planning table, best meshes per (resolution, nodes, ranks per node).

    Parameters
    ----------
    res_list : list of int
        ``RES_Z`` values.
    nodes : tuple of int, optional
        Inclusive range of node counts.
    cpt : tuple of int, optional
        ``cpus-per-task`` values (``C``); ranks per node from :func:`ranks_per_node`.
    cores : int, optional
        Cores per node.
    gamma : float, optional
        Aspect ratio.
    l3 : int, optional
        Cores per L3 cache for packed placement (0: plain packing).
    mem_ref : tuple, optional
        ``(RES_Z, total GiB)`` of a measured run; memory is scaled as ``RES_Z**3``.
    show_all : bool, optional
        Include unbalanced configurations (otherwise only efficiency 1).
    top : int, optional
        Meshes kept per configuration, best first (highest efficiency, then most square).

    Returns
    -------
    list of dict
    """
    rows = []
    for res in res_list:
        sizes = groups(res, gamma)
        for nn in range(nodes[0], nodes[1] + 1):
            for c in cpt:
                rpn = ranks_per_node(c, cores, l3)
                if rpn == 0:
                    continue
                n = nn * rpn
                ranked = sorted(meshes(n, sizes), key=lambda m: (-m[2], abs(m[0] - m[1]), m[0]))
                for p0, p1, e in ranked[:top]:
                    if e < 1 and not show_all:
                        continue
                    row = dict(RES_Z=res, nodes=nn, C=c, per_node=rpn, N=n, mesh=f"{p0}x{p1}", eff=e)
                    if mem_ref:
                        row["GiB/node"] = mem_ref[1] * (res / mem_ref[0]) ** 3 / nn
                    rows.append(row)
    return rows


def get_parser() -> argparse.ArgumentParser:
    """Command-line options (see the module docstring).

    Returns
    -------
    argparse.ArgumentParser
    """
    p = argparse.ArgumentParser(prog="python -m rbc.plan",
                                description="Balanced node counts and meshes for a resolution.")
    p.add_argument("--res", type=int, nargs="+", required=True, help="RES_Z value(s)")
    p.add_argument("--nodes", type=int, nargs=2, default=(1, 12), metavar=("MIN", "MAX"))
    p.add_argument("--cpt", type=int, nargs="+", default=[1, 2, 3], help="cpus-per-task values")
    p.add_argument("--cores", type=int, default=192, help="cores per node (Trillium: 192)")
    p.add_argument("--l3", type=int, default=8,
                   help="cores per L3 cache for packed placement (Trillium: 8; 0 = cores // C)")
    p.add_argument("--gamma", type=float, default=2.0, help="aspect ratio Lx/Lz")
    p.add_argument("--mem-ref", type=float, nargs=2, metavar=("RES_Z", "GIB"),
                   help="measured total memory of a run, scaled as RES_Z**3")
    p.add_argument("--all", action="store_true", help="also list unbalanced configurations")
    p.add_argument("--top", type=int, default=1, help="meshes per configuration (best first)")
    p.add_argument("--emit", action="store_true",
                   help="print 'RES_Z NODES C MESHX MESHY' lines for trill_bench instead of a table")
    return p


def main(argv=None) -> None:
    """Print the planning table.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments (default: ``sys.argv[1:]``).
    """
    a = get_parser().parse_args(argv)
    rows = plan(a.res, tuple(a.nodes), tuple(a.cpt), a.cores, a.gamma, a.mem_ref, a.all, a.l3, a.top)
    if a.emit:
        for r in rows:
            print(r["RES_Z"], r["nodes"], r["C"], *r["mesh"].split("x"))
        return
    if not rows:
        print("no balanced configuration in this range; try --all, other --res or --cpt")
        return
    cols = list(rows[0])
    print("  ".join(f"{c:>9}" for c in cols))
    for r in rows:
        print("  ".join(f"{r[c]:>9.3f}" if isinstance(r[c], float) else f"{r[c]:>9}" for c in cols))


if __name__ == "__main__":
    main()
