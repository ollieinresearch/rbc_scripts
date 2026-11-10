"""Two-stage benchmark driven by bench.conf: choose the resolution, then the node count.

Run from a benchmark folder holding ``params.conf`` (physics, output settings and
paths, as for a run) and ``bench.conf`` (what to benchmark; template in
``shell/bench.conf``), on a login node::

    python -m rbc.bench            # do the next step (asks before submitting anything)
    python -m rbc.bench status     # progress of the submitted jobs
    python -m rbc.bench plan       # show what the next step would submit, submit nothing
    python -m rbc.bench decide     # rewrite the decision files from the results so far

Run it repeatedly; each call does what is due:

1. plan stage 1 (resolutions at or above the floor) and submit one job per
   configuration, each on exactly the nodes its mesh fills (``shell/bench_one.sh``);
2. while jobs run, report their progress;
3. when stage 1 is done, write ``bench/res/decision_resolution.txt``, ``.pdf`` and ``.png``,
   plan stage 2 (node counts for the chosen resolution, from its measured memory) and
   submit it;
4. when stage 2 is done, write ``bench/nodes/decision_nodes.txt``, ``.pdf`` and ``.png``, which
   end with the ``RES_Z``/``RES_HOR``/``MESHX``/``MESHY``/``--nodes`` lines for production;
   if the best count is at an edge of those measured, offer counts beyond it (edge rule).

A resolution is a pair of Chebyshev modes RES_Z and horizontal modes RES_HOR (the
same in x and y); names put RES_Z first (``288_512``). Internally the pair is stored as
``(RES_HOR, RES_Z)``. By default RES_HOR follows RES_Z at the ratio in ``params.conf``;
setting ``RES_HOR_FLOOR`` in bench.conf lets the two vary independently, e.g. to
refine the boundary layers without refining horizontally.

Each benchmark job is a short production job: the production outputs are built before
the first step and it steps for 10 minutes (``BENCH_SECONDS``), the first minute
(``BENCH_DISCARD``) left out of the timing. The memory at its end, divided by
``BENCH_MEM_FRACTION`` (0.80: the share of a 24 h job's peak reached after 10 to 15
minutes), is the projected peak of a production job.

Every result with that memory measurement is also added to a history file
(``BENCH_HISTORY``, default ``~/rbc_bench_history.csv``) with the cluster it ran on.
New benchmarks on the same cluster use it: the node count chosen before at the nearest
resolution, scaled by size, is where stage 1 runs and what stage 2 sweeps around; memory
only sets the fewest nodes allowed. If the best count of stage 2 is the largest (or
smallest) measured, later runs offer counts beyond it until the best is bracketed or
``NODE_MAX`` is reached.

The rules behind each choice are written into the decision files; ``docs/benchmarking.md``
derives them.
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import re
import socket
import statistics
import subprocess
from datetime import datetime
from pathlib import Path

from rbc.plan import efficiency, meshes, ranks_per_node

NOISE = 0.03
STALL = 2.0               # a timed block slower than this multiple of the median block is a stall
SEGMENTS = 4              # the timed window is split into this many parts for the spread
OVERHEAD_S = 120          # per job besides the stepping: rehearsal, two solver builds, final writes
SLACK_MIN = 5             # minutes added to BENCH_SECONDS for the wall-time request
EDGE_STEPS = (1.25, 1.5)  # multiples of the edge count tried when the best count is at an edge
HISTORY_FIELDS = ("cluster", "recorded", "job_id", "nz", "nx", "nodes", "c", "mesh", "snapshots",
                  "s_per_iter", "noise", "mem_end_gib", "mem_fraction", "peak_gib", "var_peak_gib",
                  "var_per_cell", "fits", "mem_stop", "chosen", "folder")
TIER_PRIOR = {"A": 1.00, "B": 1.03, "C": 1.08}
TIER_NAME = {"A": "3-smooth", "B": "5-smooth", "C": "one factor of 7"}
ROYGBIV = ("#d62728", "#ff7f0e", "#e3b505", "#2ca02c", "#1f77b4", "#3f007d", "#9467bd")

DEFAULTS = dict(RES_FLOOR=None, RES_ABOVE=0.125, RES_HOR_FLOOR=None, RES_HOR_ABOVE=None,
                RES_BAND=0.15, RES_MAX=4, ORIENT="larger", SQUARE=1,
                BENCH_SECONDS=600.0, BENCH_DISCARD=60.0, BENCH_REHEARSAL=10, BENCH_TIME="",
                BENCH_MEM_FRACTION=0.80, BENCH_HISTORY="", ACCOUNT="rrg-goluskin", MEM_FIXED=90.0,
                MEM_PER_CELL=1.69e-4, MEM_GROWTH=0.015, JOB_HOURS=24.0, MEM_USABLE=747.0, MEM_HEADROOM=0.85, NODE_SPAN=2.0, NODE_POINTS=6,
                NODE_MIN=None, NODE_MAX=32, TIME_WEIGHT=1.0, DT_REF=None, DT_REF_RES=None,
                DT_REF_RES_HOR=None, MEM_FACTOR=1.0)
PARAMS = ("RES_Z", "RES_HOR", "LX", "RA", "PR", "SIM_TIME", "RESTART_FROM", "SCRIPTS_PATH", "USER",
          "SNAPSHOTS", "CC_CLUSTER")
INTS = ("RES_FLOOR", "RES_HOR_FLOOR", "RES_MAX", "BENCH_REHEARSAL", "NODE_POINTS",
        "NODE_MIN", "NODE_MAX", "DT_REF_RES", "DT_REF_RES_HOR", "SQUARE")


# ---------------------------------------------------------------- configuration

def source(folder: Path, files=("params.conf", "bench.conf")) -> dict:
    """Environment after sourcing configuration files in bash, as the jobs do.

    Parameters
    ----------
    folder : pathlib.Path
        Folder holding the files.
    files : tuple of str, optional
        Files sourced in order (missing ones skipped).

    Returns
    -------
    dict
        Variable name -> value.
    """
    cmd = "set -a; " + " ".join(f'[ -f "{f}" ] && source "{f}" >/dev/null;' for f in files) + " env -0"
    out = subprocess.run(["bash", "-c", cmd], cwd=folder, capture_output=True, check=True).stdout
    return dict(item.split("=", 1) for item in out.decode().split("\0") if "=" in item)


def load_config(folder: Path) -> dict:
    """``bench.conf`` settings (with defaults) plus the ``params.conf`` values used here.

    Parameters
    ----------
    folder : pathlib.Path
        Benchmark folder.

    Returns
    -------
    dict
        Also holds ``RATIO`` (RES_HOR / RES_Z when they are linked), ``LINKED`` and
        ``FLOOR`` (the floor resolution pair).

    Raises
    ------
    SystemExit
        Without ``params.conf``, or without a floor resolution.
    """
    if not (folder / "params.conf").exists():
        raise SystemExit(f"no params.conf in {folder}; copy one from a run (see shell/bench.conf)")
    env = source(folder)
    cfg = {}
    for key, default in DEFAULTS.items():
        raw = env.get(key, "")
        if raw == "":
            cfg[key] = default
        elif isinstance(default, str):
            cfg[key] = raw
        else:
            cfg[key] = float(raw)
    for key in PARAMS:
        cfg[key] = env.get(key) or None
    if cfg["RES_FLOOR"] is None:
        if not cfg["RES_Z"]:
            raise SystemExit("set RES_FLOOR in bench.conf (or RES_Z in params.conf)")
        cfg["RES_FLOOR"] = cfg["RES_Z"]
    for key in INTS:
        if cfg[key] is not None:
            cfg[key] = int(float(cfg[key]))
    cfg["SIM_TIME"] = float(cfg["SIM_TIME"]) if cfg["SIM_TIME"] else None
    cfg["LX"] = float(cfg["LX"]) if cfg["LX"] else 2.0
    if cfg["RES_HOR"] and cfg["RES_Z"] and int(cfg["RES_Z"]):
        cfg["RATIO"] = int(cfg["RES_HOR"]) / int(cfg["RES_Z"])
    else:
        cfg["RATIO"] = 2.0
    cfg["LINKED"] = cfg["RES_HOR_FLOOR"] is None
    if cfg["RES_HOR_ABOVE"] is None:
        cfg["RES_HOR_ABOVE"] = cfg["RES_ABOVE"]
    nz = cfg["RES_FLOOR"]
    cfg["FLOOR"] = (linked_hor(nz, cfg) if cfg["LINKED"] else cfg["RES_HOR_FLOOR"], nz)
    return cfg


def linked_hor(nz: int, cfg: dict) -> int:
    """Horizontal resolution that goes with ``nz`` when the two are linked.

    Parameters
    ----------
    nz : int
        RES_Z.
    cfg : dict
        Configuration (``RATIO``).

    Returns
    -------
    int
        Even, ``RATIO x nz`` rounded.
    """
    nx = round(cfg["RATIO"] * nz)
    return nx + nx % 2


# ---------------------------------------------------------------- arithmetic

def factorise(n: int) -> dict[int, int]:
    """Prime factorisation.

    Parameters
    ----------
    n : int
        Positive integer.

    Returns
    -------
    dict
        Prime -> exponent.
    """
    f, p = {}, 2
    while p * p <= n:
        while n % p == 0:
            f[p] = f.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def factor_str(n: int) -> str:
    """Factorisation as text, e.g. ``2^4 3^2``.

    Parameters
    ----------
    n : int
        Positive integer.

    Returns
    -------
    str
    """
    return " ".join(f"{p}^{e}" if e > 1 else f"{p}" for p, e in sorted(factorise(n).items()))


def tier(n: int) -> str | None:
    """Transform-size class of one axis.

    The grids (3/2 of the modes) contain the primes of the mode count plus a 3, so
    classifying the mode count classifies the transforms.

    Parameters
    ----------
    n : int
        Modes on the axis (must be even).

    Returns
    -------
    str or None
        ``"A"`` 3-smooth, ``"B"`` 5-smooth, ``"C"`` 7-smooth with a single 7, or None
        (odd, a prime of 11 or more, or 7^2).
    """
    if n % 2:
        return None
    f = factorise(n)
    if any(p > 7 for p in f) or f.get(7, 0) > 1:
        return None
    if 7 in f:
        return "C"
    return "B" if 5 in f else "A"


def pair_tier(res) -> str | None:
    """Worse tier of the two axes of a resolution.

    Parameters
    ----------
    res : tuple of int
        ``(RES_HOR, RES_Z)``.

    Returns
    -------
    str or None
    """
    t = [tier(n) for n in res]
    return None if None in t else max(t)


def label(res) -> str:
    """Short name, RES_Z first: ``"256"`` when RES_HOR = 2 RES_Z, otherwise ``"288_512"``.

    Parameters
    ----------
    res : tuple of int
        ``(RES_HOR, RES_Z)``.

    Returns
    -------
    str
    """
    nx, nz = res
    return f"{nz}" if nx == 2 * nz else f"{nz}_{nx}"


def grp(res) -> tuple[int, ...]:
    """Group counts that a mesh dimension may have to split.

    Parameters
    ----------
    res : tuple of int
        ``(RES_HOR, RES_Z)``.

    Returns
    -------
    tuple of int
        Fourier pairs, Fourier grid points, Chebyshev modes, Chebyshev grid points.
    """
    nx, nz = res
    return nx // 2, 3 * nx // 2, nz, 3 * nz // 2


def balance(mx: int, my: int, res) -> float:
    """Predicted load balance of a mesh (1 = every rank has the same work).

    Parameters
    ----------
    mx, my : int
        Mesh.
    res : tuple of int
        ``(RES_HOR, RES_Z)``.

    Returns
    -------
    float
    """
    return efficiency(mx, grp(res)) * efficiency(my, grp(res))


def best_mesh(res, nodes: int, c: int = 1):
    """Most balanced mesh for a resolution on a node count.

    Parameters
    ----------
    res : tuple of int
        ``(RES_HOR, RES_Z)``.
    nodes : int
        Node count.
    c : int, optional
        cpus-per-task.

    Returns
    -------
    tuple
        ``(meshx, meshy, balance)`` with ``meshx >= meshy``.
    """
    n = nodes * ranks_per_node(c)
    p0, p1, e = max(meshes(n, grp(res)), key=lambda m: (round(m[2], 4), -abs(m[0] - m[1])))
    return max(p0, p1), min(p0, p1), e


def square_mesh(res, nodes: int, c: int = 1):
    """Most nearly square mesh on a node count, whatever its balance.

    Parameters
    ----------
    res : tuple of int
        ``(RES_HOR, RES_Z)``.
    nodes : int
        Node count.
    c : int, optional
        cpus-per-task.

    Returns
    -------
    tuple
        ``(meshx, meshy, balance)`` with ``meshx >= meshy``.
    """
    n = nodes * ranks_per_node(c)
    p1 = max(d for d in range(1, math.isqrt(n) + 1) if n % d == 0)
    p0 = n // p1
    return p0, p1, balance(p0, p1, res)


def is_square(mx: int, my: int) -> bool:
    """Whether a mesh is the most nearly square factorisation of its rank count.

    Parameters
    ----------
    mx, my : int
        Mesh dimensions.

    Returns
    -------
    bool
    """
    n = mx * my
    return min(mx, my) == max(d for d in range(1, math.isqrt(n) + 1) if n % d == 0)


def work(res) -> float:
    """Work per step, proportional to grid points times a logarithm.

    Parameters
    ----------
    res : tuple of int
        ``(RES_HOR, RES_Z)``.

    Returns
    -------
    float
        ``RES_HOR^2 RES_Z ln(RES_HOR^2 RES_Z)``.
    """
    nx, nz = res
    return nx * nx * nz * math.log(nx * nx * nz)


def steps(res, lx: float) -> float:
    """Steps per unit simulated time, up to a constant.

    The advective CFL limit makes dt proportional to the smallest grid spacing; with
    the horizontal spacing Lx/RES_HOR and a mean vertical spacing 1/RES_Z, the number
    of steps is proportional to ``max(RES_HOR / Lx, RES_Z)``. For RES_HOR = Lx RES_Z
    this is RES_Z, as before.

    Parameters
    ----------
    res : tuple of int
        ``(RES_HOR, RES_Z)``.
    lx : float
        Box length (H = 1).

    Returns
    -------
    float
    """
    nx, nz = res
    return max(nx / lx, nz)


def cost_model(res, ref, lx: float) -> float:
    """Work to reach a fixed simulation time, relative to resolution ``ref``.

    Parameters
    ----------
    res, ref : tuple of int
        Resolutions.
    lx : float
        Box length.

    Returns
    -------
    float
    """
    return work(res) * steps(res, lx) / (work(ref) * steps(ref, lx))


def cells(res) -> float:
    """Size in units of RES_Z^3 at RES_HOR = 2 RES_Z: ``(RES_HOR/2)^2 RES_Z``.

    Parameters
    ----------
    res : tuple of int
        ``(RES_HOR, RES_Z)``.

    Returns
    -------
    float
    """
    nx, nz = res
    return (nx / 2) ** 2 * nz


# ---------------------------------------------------------------- memory

def fixed_per_node(cfg: dict, c: int = 1) -> float:
    """Resolution-independent memory per node (libraries, MPI, Python), GiB.

    Parameters
    ----------
    cfg : dict
        Configuration (``MEM_FIXED`` is for 192 ranks per node).
    c : int, optional
        cpus-per-task.

    Returns
    -------
    float
    """
    return cfg["MEM_FIXED"] * ranks_per_node(c) / 192


def growth(cfg: dict) -> float:
    """Growth of the resolution-dependent memory over a job, for results without a 10 min run.

    Applies to the ``MEM_PER_CELL`` prior and to results of the earlier, short benchmark;
    results of the 10 minute benchmark use ``BENCH_MEM_FRACTION`` instead.

    Parameters
    ----------
    cfg : dict
        Configuration (``MEM_GROWTH`` per hour, ``JOB_HOURS``).

    Returns
    -------
    float
    """
    return 1 + cfg["MEM_GROWTH"] * cfg["JOB_HOURS"]


def budget(cfg: dict) -> float:
    """Memory per node to plan for, GiB.

    Parameters
    ----------
    cfg : dict
        Configuration.

    Returns
    -------
    float
    """
    return cfg["MEM_USABLE"] * cfg["MEM_HEADROOM"]


def min_nodes(variable: float, cfg: dict, c: int = 1) -> int:
    """Fewest nodes that hold ``variable`` GiB plus the fixed part within the budget.

    Parameters
    ----------
    variable : float
        Resolution-dependent memory over all nodes at the peak of a job, GiB.
    cfg : dict
        Configuration.
    c : int, optional
        cpus-per-task.

    Returns
    -------
    int
    """
    return max(1, math.ceil(variable / (budget(cfg) - fixed_per_node(cfg, c))))


def prior_variable(res, cfg: dict) -> float:
    """Resolution-dependent memory at a job's peak, before any measurement, GiB.

    Parameters
    ----------
    res : tuple of int
        ``(RES_HOR, RES_Z)``.
    cfg : dict
        Configuration (``MEM_PER_CELL`` GiB per unit of :func:`cells`).

    Returns
    -------
    float
    """
    return cfg["MEM_PER_CELL"] * cells(res) * growth(cfg)


def measured_variable(row: dict, cfg: dict) -> tuple[float, float, str]:
    """Resolution-dependent memory of a benchmarked configuration, now and at a job's peak, GiB.

    Per node, the memory measured at the end of the benchmark is the larger of summed
    resident memory and the node's ``MemTotal - MemAvailable``, on the fullest node.
    Over all nodes, the resolution-dependent part is that minus the fixed part, times
    the node count. Projected to a production job's peak:

    * 10 minute benchmark (production outputs from the first step): the per-node
      memory divided by ``BENCH_MEM_FRACTION``, the share of a job's peak reached
      after 10 to 15 minutes (about 0.80 on 24 h jobs, from production memory logs);
    * earlier short benchmark (outputs built at the end): times ``MEM_FACTOR``
      (production over that benchmark) and :func:`growth`;
    * no memory measured: :func:`prior_variable`.

    Parameters
    ----------
    row : dict
        From :func:`load`.
    cfg : dict
        Configuration.

    Returns
    -------
    peak : float
        GiB over all nodes at the peak of a job.
    now : float
        GiB over all nodes at the end of the benchmark (scaled by ``MEM_FACTOR`` for
        the short benchmark).
    how : str
        Where the number came from.
    """
    n, fixed = row["nodes"], fixed_per_node(cfg, row["c"])
    if row["mem_io"] is None:
        v = prior_variable(row["res"], cfg)
        return v, v / growth(cfg), "prior (no memory measured with outputs in this result)"
    if row["long"]:
        f = cfg["BENCH_MEM_FRACTION"]
        return ((row["mem_io"] / f - fixed) * n, (row["mem_io"] - fixed) * n,
                f"{row['mem_io']:.0f} GiB/node after {row['stepping_s'] / 60:.0f} min on {n} nodes, "
                f"/ BENCH_MEM_FRACTION {f:g}")
    now = (row["mem_io"] - fixed) * n * cfg["MEM_FACTOR"]
    scaled = f", x MEM_FACTOR {cfg['MEM_FACTOR']:g}" if cfg["MEM_FACTOR"] != 1 else ""
    return (now * growth(cfg), now, f"short benchmark: {row['mem_io']:.0f} GiB/node on {n} nodes{scaled}, "
                                    f"grown x{growth(cfg):.2f} over {cfg['JOB_HOURS']:g} h")


def memory_rule(rows, cfg: dict) -> str:
    """One sentence on how the projected peak memory was obtained for these rows.

    Parameters
    ----------
    rows : list of dict
        From :func:`load`.
    cfg : dict
        Configuration.

    Returns
    -------
    str
    """
    parts = []
    if any(r["long"] for r in rows):
        parts.append(f"memory at the end of the {cfg['BENCH_SECONDS'] / 60:.0f} min benchmark (outputs built "
                     f"from the start) / BENCH_MEM_FRACTION {cfg['BENCH_MEM_FRACTION']:g}")
    if any(not r["long"] and r["mem_io"] is not None for r in rows):
        parts.append(f"short-benchmark memory x MEM_FACTOR {cfg['MEM_FACTOR']:g}, grown x{growth(cfg):.2f} "
                     f"over {cfg['JOB_HOURS']:g} h")
    if any(r["mem_io"] is None for r in rows):
        parts.append(f"MEM_PER_CELL prior, grown x{growth(cfg):.2f}")
    return "projected peak per node = " + "; or ".join(parts) + \
        f"; budget {budget(cfg):.0f} of {cfg['MEM_USABLE']:.0f} GiB/node"


def job_peak(variable: float, nodes: int, cfg: dict, c: int = 1) -> float:
    """Projected peak memory per node of a production job, GiB.

    Parameters
    ----------
    variable : float
        Resolution-dependent memory over all nodes at the peak, GiB (from
        :func:`measured_variable`).
    nodes : int
        Node count.
    cfg : dict
        Configuration.
    c : int, optional
        cpus-per-task.

    Returns
    -------
    float
    """
    return fixed_per_node(cfg, c) + variable / nodes


# ---------------------------------------------------------------- history across benchmarks

def cluster_name(cfg: dict | None = None) -> str:
    """Name of this cluster: ``CC_CLUSTER``, else Slurm's ClusterName, else the host name.

    Parameters
    ----------
    cfg : dict, optional
        Configuration (``CC_CLUSTER`` from the environment of the config files).

    Returns
    -------
    str
    """
    name = (cfg or {}).get("CC_CLUSTER") or os.environ.get("CC_CLUSTER")
    if name:
        return name
    try:
        out = subprocess.run(["scontrol", "show", "config"], capture_output=True, text=True, timeout=20).stdout
        m = re.search(r"ClusterName\s*=\s*(\S+)", out)
        if m:
            return m.group(1)
    except (OSError, subprocess.SubprocessError):
        pass
    return re.sub(r"[-\d]*(login|nia|tri)?\d*$", "", socket.gethostname().split(".")[0]) or "unknown"


def history_path(cfg: dict) -> Path:
    """Where the benchmark history lives (``BENCH_HISTORY``, default ``~/rbc_bench_history.csv``).

    Parameters
    ----------
    cfg : dict
        Configuration.

    Returns
    -------
    pathlib.Path
    """
    return Path(cfg["BENCH_HISTORY"] or "~/rbc_bench_history.csv").expanduser()


def read_history(cfg: dict) -> list[dict]:
    """Rows of the history file (empty if there is none).

    Parameters
    ----------
    cfg : dict
        Configuration.

    Returns
    -------
    list of dict
    """
    path = history_path(cfg)
    if not path.exists():
        return []
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def record_history(rows, cfg: dict, folder: Path, chosen=None) -> int:
    """Add the 10 minute benchmark results to the history; mark the stage 2 choice.

    Only results of the 10 minute benchmark are recorded (the short benchmark's memory
    is not comparable). A result already there (same cluster, job id and
    configuration) is replaced, so re-running a decision does not duplicate it.

    Parameters
    ----------
    rows : list of dict
        From :func:`load` (or :func:`decide_nodes`).
    cfg : dict
        Configuration.
    folder : pathlib.Path
        Benchmark folder, recorded for reference.
    chosen : dict, optional
        The configuration chosen in stage 2; its resolution's other rows get chosen = 0.

    Returns
    -------
    int
        Rows added or updated.
    """
    cluster = cluster_name(cfg)
    key = lambda h: (h["cluster"], h["job_id"], h["nz"], h["nx"], h["nodes"], h["mesh"])
    hist = {key(h): h for h in read_history(cfg)}
    n = 0
    for r in rows:
        if not r.get("long"):
            continue
        peak, _, _ = measured_variable(r, cfg)
        h = dict(cluster=r.get("cluster") or cluster, recorded=datetime.now().isoformat(timespec="seconds"),
                 job_id=r["job_id"], nz=str(r["nz"]), nx=str(r["nx"]), nodes=str(r["nodes"]), c=str(r["c"]),
                 mesh=r["mesh"], snapshots=str(r.get("snapshots", "")), s_per_iter=f"{r['s']:.5f}",
                 noise=f"{r['noise']:.4f}" if r["noise"] is not None else "", mem_end_gib=f"{r['mem_io']:.1f}",
                 mem_fraction=f"{cfg['BENCH_MEM_FRACTION']:g}",
                 peak_gib=f"{job_peak(peak, r['nodes'], cfg, r['c']):.1f}", var_peak_gib=f"{peak:.1f}",
                 var_per_cell=f"{peak / cells(r['res']):.6g}",
                 fits=str(int(job_peak(peak, r["nodes"], cfg, r["c"]) <= budget(cfg) and not r["mem_stop"])),
                 mem_stop=str(int(r["mem_stop"])), chosen="", folder=str(folder))
        old = hist.get(key(h))
        if old and old.get("chosen"):
            h["chosen"] = old["chosen"]
        if chosen is not None and r["res"] == chosen["res"]:
            h["chosen"] = "1" if r["nodes"] == chosen["nodes"] and r["mesh"] == chosen["mesh"] else "0"
        hist[key(h)] = h
        n += 1
    if n:
        path = history_path(cfg)
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=HISTORY_FIELDS, extrasaction="ignore")
            w.writeheader()
            w.writerows(sorted(hist.values(), key=lambda h: (h["cluster"], int(h["nz"]), int(h["nx"]),
                                                              int(h["nodes"]))))
    return n


def _nearest(hist, res):
    """History rows at the resolution closest to ``res`` in size (log of :func:`cells`)."""
    if not hist:
        return []
    size = lambda h: cells((int(h["nx"]), int(h["nz"])))
    best = min(hist, key=lambda h: abs(math.log(size(h) / cells(res))))
    return [h for h in hist if (h["nx"], h["nz"]) == (best["nx"], best["nz"])]


def history_variable(res, cfg: dict, c: int = 1):
    """Peak resolution-dependent memory for ``res``, scaled from the history, GiB over all nodes.

    From the resolution closest in size measured before on this cluster with the same
    cpus-per-task (C = 2 used about twice the memory of C = 1) and, when known, the same
    snapshot setting. Of its results, the largest memory per unit of :func:`cells` is
    used; the resolution-dependent memory is assumed proportional to :func:`cells`
    (field and transform buffers are), while the fixed part per node is not scaled.

    Parameters
    ----------
    res : tuple of int
        ``(RES_HOR, RES_Z)``.
    cfg : dict
        Configuration.
    c : int, optional
        cpus-per-task.

    Returns
    -------
    tuple or None
        ``(GiB, explanation)``, or None without a usable history.
    """
    cluster, snaps = cluster_name(cfg), str(cfg.get("SNAPSHOTS") or "")
    hist = [h for h in read_history(cfg) if h["cluster"] == cluster and h["c"] == str(c) and h["mem_stop"] != "1"
            and h.get("var_per_cell")]
    if snaps:
        hist = [h for h in hist if h["snapshots"] in ("", snaps)] or hist
    near = _nearest(hist, res)
    if not near:
        return None
    h = max(near, key=lambda h: float(h["var_per_cell"]))
    old = (int(h["nx"]), int(h["nz"]))
    v = float(h["var_per_cell"]) * cells(res)
    return v, (f"history ({cluster}): {label(old)} on {h['nodes']} nodes, {float(h['var_peak_gib']):,.0f} GiB "
               f"at the peak, x{cells(res) / cells(old):.2f} for the size")


def history_optimum(res, cfg: dict, c: int = 1):
    """Node count chosen before at the nearest resolution, scaled to ``res`` by size.

    Work per step grows like :func:`cells` (times a logarithm), so a count that was
    optimal at one size is scaled by the ratio of :func:`cells` as a first guess.

    Parameters
    ----------
    res : tuple of int
        ``(RES_HOR, RES_Z)``.
    cfg : dict
        Configuration.
    c : int, optional
        cpus-per-task.

    Returns
    -------
    tuple or None
        ``(nodes, explanation)`` with ``nodes`` unrounded, or None if no choice is
        recorded on this cluster.
    """
    cluster = cluster_name(cfg)
    hist = [h for h in read_history(cfg) if h["cluster"] == cluster and h["c"] == str(c) and h["chosen"] == "1"]
    near = _nearest(hist, res)
    if not near:
        return None
    h = near[-1]
    old = (int(h["nx"]), int(h["nz"]))
    n = int(h["nodes"]) * cells(res) / cells(old)
    return n, (f"{h['nodes']} nodes chosen at {label(old)} ({cluster}), x{cells(res) / cells(old):.2f} "
               f"for the size = {n:.1f}")


def even(x: float, floor: int = 2) -> int:
    """Nearest even node count to ``x``, at least ``floor`` (rounded up to even).

    Even counts are preferred because a 128 x M mesh needs 192 N / 128 = 1.5 N ranks per
    row, an integer only for even N.

    Parameters
    ----------
    x : float
        Target count.
    floor : int, optional
        Smallest count allowed.

    Returns
    -------
    int
    """
    return max(2 * round(x / 2), floor + floor % 2)


# ---------------------------------------------------------------- timestep for projections

def estimate_dt(cfg: dict):
    """Mean production timestep at a known resolution, to scale with :func:`steps`.

    From ``DT_REF`` (at ``DT_REF_RES_HOR`` x ``DT_REF_RES``, default the floor) in
    bench.conf, otherwise the mean of the last analysis writes of the run in
    ``RESTART_FROM``.

    Parameters
    ----------
    cfg : dict
        Configuration.

    Returns
    -------
    tuple or None
        ``(dt, res, source)``, or None if neither is available.
    """
    if cfg["DT_REF"]:
        nz = cfg["DT_REF_RES"] or cfg["FLOOR"][1]
        nx = cfg["DT_REF_RES_HOR"] or (cfg["FLOOR"][0] if not cfg["DT_REF_RES"] else linked_hor(nz, cfg))
        return cfg["DT_REF"], (nx, nz), "DT_REF in bench.conf"
    src = cfg.get("RESTART_FROM")
    if not src or not (Path(src) / "analysis").is_dir():
        return None
    try:
        import h5py
        import numpy as np

        files = sorted((Path(src) / "analysis").glob("analysis_s*.h5"),
                       key=lambda p: int(re.findall(r"_s(\d+)", p.name)[-1]))[-3:]
        ts = []
        for f in files:
            with h5py.File(f, "r") as h:
                ts += list(h["scales/timestep"][:])
        env = source(Path(src), ("params.conf",))
        nz = int(env["RES_Z"])
        nx = int(env.get("RES_HOR") or 2 * nz)
        return float(np.mean(ts)), (nx, nz), f"mean of the last {len(ts)} timesteps of {src}"
    except Exception:
        return None


def projection(s_per_iter: float, res, cfg: dict):
    """Wall time per production job (``SIM_TIME`` time units), hours.

    Parameters
    ----------
    s_per_iter : float
        Seconds per step.
    res : tuple of int
        ``(RES_HOR, RES_Z)``.
    cfg : dict
        Configuration.

    Returns
    -------
    float or None
    """
    est = estimate_dt(cfg)
    if est is None or not cfg["SIM_TIME"]:
        return None
    dt = est[0] * steps(est[1], cfg["LX"]) / steps(res, cfg["LX"])
    return cfg["SIM_TIME"] / dt * s_per_iter / 3600


# ---------------------------------------------------------------- results

def load(path) -> list[dict]:
    """Benchmark rows from a CSV file or every ``bench.csv`` under a folder.

    Repeated configurations are merged: ``s_per_iter`` is averaged and the noise
    estimate is the larger of the within-run block spread and the between-run range.

    Parameters
    ----------
    path : str or pathlib.Path
        A CSV file or a folder.

    Returns
    -------
    list of dict
        One per configuration, with numeric fields converted; ``res`` is
        ``(RES_HOR, RES_Z)``.
    """
    path = Path(path)
    files = [path] if path.is_file() else sorted(path.rglob("bench.csv")) if path.exists() else []
    raw, f_of = [], {}
    for f in files:
        with open(f, newline="") as h:
            here = list(csv.DictReader(h))
        long = [r for r in here if r.get("discard_s") and str(r.get("outputs")) == "1"]
        runs = trace_runs(f.parent)
        if long and len(runs) == len(long):          # step time recomputed without stalls (block_stats)
            for r, blocks in zip(long, runs):
                st = block_stats(blocks, float(r["discard_s"]))
                r.update(s_per_iter=f"{st['s_iter']:.5f}", s_per_iter_spread=f"{st['spread']:.4f}",
                         iters=str(st["iters"]), stall_blocks=str(st["stall_blocks"]), stall_s=f"{st['stall_s']:.1f}")
        for r in here:
            raw.append(r)
            f_of[id(r)] = f.parent
    merged = {}
    for r in raw:
        if str(r.get("finite", "True")) != "True":
            continue
        key = (int(r["nx"]), int(r["nz"]), int(r["nodes"]), int(r["cpus_per_rank"]), r["mesh"])
        merged.setdefault(key, []).append(r)

    def num(r, k):
        v = r.get(k)
        return float(v) if v not in (None, "", "nan") else None

    rows = []
    for (nx, nz, nodes, c, mesh), rs in merged.items():
        s = [float(r["s_per_iter"]) for r in rs]
        mean = sum(s) / len(s)
        spreads = [float(r["s_per_iter_spread"]) for r in rs
                   if r.get("s_per_iter_spread") and math.isfinite(float(r["s_per_iter_spread"]))]
        noise = max(spreads + [(max(s) - min(s)) / mean if len(s) > 1 else 0.0])
        last = rs[-1]
        io = [v for v in (num(last, "node_rss_io_max_gib"), num(last, "node_used_io_max_gib")) if v]
        mx, my = (int(v) for v in mesh.split("x"))
        res = (nx, nz)
        long = bool(last.get("discard_s")) and str(last.get("outputs")) == "1"
        rows.append(dict(nx=nx, nz=nz, res=res, label=label(res), nodes=nodes, c=c, mesh=mesh, mx=mx, my=my,
                         s=mean, runs=len(s), noise=noise if (spreads or len(s) > 1) else None,
                         node_s=nodes * mean, mem_max=float(last["node_rss_max_gib"]),
                         mem_min=float(last["node_rss_min_gib"]), mem_io=max(io) if io else None,
                         eff=balance(mx, my, res), long=long, stepping_s=num(last, "stepping_s") or 0.0,
                         iters=int(float(last.get("iters") or 0)), mem_stop=str(last.get("mem_stop")) == "1",
                         truncated=str(last.get("truncated")) == "1", job_id=last.get("job_id", ""),
                         cluster=last.get("cluster", ""), snapshots=last.get("snapshots", ""),
                         mem_discard=num(last, "node_used_discard_max_gib"), source=str(f_of[id(last)]),
                         stall_s=sum(num(r, "stall_s") or 0.0 for r in rs)))
    return sorted(rows, key=lambda r: (cells(r["res"]), r["nodes"], r["node_s"]))


def block_stats(blocks, discard: float, stall: float = STALL) -> dict:
    """Step time of a benchmark from its timed blocks, leaving out stalls.

    Blocks starting before ``discard`` seconds are not timed. Of the rest, a block
    whose time per step is more than ``stall`` times the median block's is a stall
    and is left out: most of it was spent waiting (a file-system or network hiccup),
    not stepping. Timing delays only ever add time (Chen & Revels 2016), so only slow
    blocks are rejected. Routine output writes make a block only somewhat slower
    (up to about 1.3x the median in the runs so far) and stay in, since production
    pays for them. Stalls are left out because one or two in a 9 minute window
    decide the comparison between node counts by chance; their total is reported.

    Parameters
    ----------
    blocks : list of tuple
        ``(start s, end s, steps)`` per block, in order.
    discard : float
        Seconds at the start not timed.
    stall : float, optional
        Rejection threshold, as a multiple of the median block's time per step.

    Returns
    -------
    dict
        ``s_iter`` (s per step over the kept blocks), ``spread`` ((max - min) / mean
        of the means of :data:`SEGMENTS` equal parts of the kept window), ``iters``,
        ``measured`` (s), ``blocks`` (kept), ``discard_iters``, ``stall_blocks``,
        ``stall_s``.
    """
    timed = [b for b in blocks if b[0] >= discard] or blocks[len(blocks) // 2:]
    rate = [(b[1] - b[0]) / b[2] for b in timed]
    med = statistics.median(rate)
    kept = [b for b, r in zip(timed, rate) if r <= stall * med]
    dropped = [b for b, r in zip(timed, rate) if r > stall * med]
    measured = sum(b[1] - b[0] for b in kept)
    iters = sum(b[2] for b in kept)
    s_iter = measured / iters
    lo, hi = kept[0][0], kept[-1][1]
    seg = []
    for i in range(SEGMENTS):
        a, b = lo + (hi - lo) * i / SEGMENTS, lo + (hi - lo) * (i + 1) / SEGMENTS
        part = [blk for blk in kept if a <= (blk[0] + blk[1]) / 2 < b or (i == SEGMENTS - 1 and blk[1] == hi)]
        if part:
            seg.append(sum(blk[1] - blk[0] for blk in part) / sum(blk[2] for blk in part))
    return dict(s_iter=s_iter, spread=(max(seg) - min(seg)) / s_iter if len(seg) > 1 else float("nan"),
                iters=iters, measured=measured, blocks=len(kept),
                discard_iters=sum(b[2] for b in blocks) - sum(b[2] for b in timed),
                stall_blocks=len(dropped), stall_s=sum(b[1] - b[0] for b in dropped))


def trace_runs(folder: Path) -> list[list[tuple]]:
    """Timed blocks of each run in a configuration's ``bench_mem.csv``.

    Each line of the file is the end of one block (time from the first step,
    iteration); a run's last line, written after the final output writes at the same
    iteration, is not a block. A rerun into the same folder appends a new run, whose
    times start again from zero.

    Parameters
    ----------
    folder : pathlib.Path
        Configuration folder.

    Returns
    -------
    list of list of tuple
        Per run, ``(start s, end s, steps)`` for every block after the first (whose
        start is not recorded; it lies in the untimed start).
    """
    path = folder / "bench_mem.csv"
    if not path.exists():
        return []
    runs, prev = [], None
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            t, it = float(r["t_s"]), int(float(r["iteration"]))
            if prev is None or t < prev[0]:
                runs.append([])
            elif it > prev[1]:
                runs[-1].append((prev[0], t, it - prev[1]))
            prev = (t, it)
    return [r for r in runs if r]


def noise_of(*rows) -> float:
    """Relative noise to compare rows against (at least :data:`NOISE`).

    Parameters
    ----------
    *rows : dict
        Rows from :func:`load`.

    Returns
    -------
    float
    """
    return max([NOISE] + [r["noise"] for r in rows if r["noise"] is not None])


def line_of(row) -> tuple:
    """Configuration line of a result row.

    Parameters
    ----------
    row : dict
        From :func:`load`.

    Returns
    -------
    tuple
        ``(RES_Z, nodes, C, MESHX, MESHY, RES_HOR)``.
    """
    return row["nz"], row["nodes"], row["c"], row["mx"], row["my"], row["nx"]


# ---------------------------------------------------------------- stage 1

def axis_values(floor: int, above: float) -> list[int]:
    """Even sizes from ``floor`` to ``floor (1 + above)``.

    Parameters
    ----------
    floor : int
        Lowest size.
    above : float
        Fraction above.

    Returns
    -------
    list of int
    """
    hi = math.floor(floor * (1 + above))
    return list(range(floor + floor % 2, hi + 1, 2))


def plan_resolutions(cfg: dict):
    """Resolutions at or above the floor worth benchmarking, each with a node count and mesh.

    Candidates: with RES_HOR linked to RES_Z, every even RES_Z in
    ``[RES_FLOOR, RES_FLOOR (1 + RES_ABOVE)]``; with ``RES_HOR_FLOOR`` set, every pair
    of such a RES_Z and an even RES_HOR in ``[RES_HOR_FLOOR, RES_HOR_FLOOR (1 +
    RES_HOR_ABOVE)]``. A pair is classified by the worse :func:`tier` of its axes. Each
    acceptable one is run on the node count chosen before at the nearest resolution on
    this cluster, scaled by size (:func:`history_optimum`) and rounded to an even count,
    but not below its memory minimum. Without such a choice it gets the most balanced
    count between its memory minimum and 1.35 times that. The memory minimum comes from
    the history (:func:`history_variable`) when this cluster has one, otherwise from the
    ``MEM_PER_CELL`` prior. Its predicted
    cost is :func:`cost_model` times the tier prior over the mesh balance. Kept: the
    floor, the best-predicted resolution of each tier, and nothing predicted more than
    ``RES_BAND`` above the cheapest; at most ``RES_MAX``.

    Parameters
    ----------
    cfg : dict
        Configuration.

    Returns
    -------
    rows : list of dict
        Every resolution examined, with ``keep`` and ``why``.
    lines : list of tuple
        ``(RES_Z, nodes, 1, MESHX, MESHY, RES_HOR)`` to benchmark.
    """
    floor = cfg["FLOOR"]
    zs = axis_values(floor[1], cfg["RES_ABOVE"])
    pairs = [(linked_hor(z, cfg), z) for z in zs] if cfg["LINKED"] else \
        [(x, z) for z in zs for x in axis_values(floor[0], cfg["RES_HOR_ABOVE"])]
    rows = []
    for res in pairs:
        t = pair_tier(res)
        row = dict(res=res, tier=t, keep=False, why="")
        if t is None and res != floor:
            row["why"] = "rejected: a prime of 11 or more, or 7^2"
            rows.append(row)
            continue
        hv = history_variable(res, cfg)
        n0 = min_nodes(hv[0] if hv else prior_variable(res, cfg), cfg)
        row["mem_how"] = hv[1] if hv else f"MEM_PER_CELL prior, grown x{growth(cfg):.2f}"
        seed = history_optimum(res, cfg)
        if seed:
            nn = even(seed[0], n0)
            if cfg["NODE_MAX"]:
                nn = max(min(nn, cfg["NODE_MAX"]), n0)
            mx, my, e = best_mesh(res, nn)
            row["nodes_how"] = f"{seed[1]} -> {nn}"
        else:
            best = None
            for nn in range(n0, max(n0, math.ceil(1.35 * n0)) + 1):
                mx, my, e = best_mesh(res, nn)
                key = (round(e, 3), -nn)
                if best is None or key > best[0]:
                    best = (key, nn, mx, my, e)
            _, nn, mx, my, e = best
            row["nodes_how"] = "no choice in the history: memory minimum"
        row.update(nodes=nn, n0=n0, mx=mx, my=my, eff=e,
                   pred=cost_model(res, floor, cfg["LX"]) * TIER_PRIOR.get(t, 1.15) / e)
        rows.append(row)

    ok = [r for r in rows if "pred" in r]
    cheapest = min(r["pred"] for r in ok)
    picks = {floor}
    for t in "ABC":
        tier_rows = [r for r in ok if r["tier"] == t]
        if tier_rows:
            picks.add(min(tier_rows, key=lambda r: r["pred"])["res"])
    for r in ok:
        if r["res"] not in picks:
            r["why"] = "not the best of its tier"
        elif r["res"] != floor and r["pred"] > (1 + cfg["RES_BAND"]) * cheapest:
            r["why"] = f"predicted more than {cfg['RES_BAND']:.0%} above the cheapest: prediction decides"
            picks.discard(r["res"])
    ranked = sorted((r for r in ok if r["res"] in picks), key=lambda r: (r["res"] != floor, r["pred"]))
    for i, r in enumerate(ranked):
        if i < cfg["RES_MAX"]:
            r["keep"] = True
            r["why"] = "floor" if r["res"] == floor else f"best of tier {r['tier']}"
            if r["tier"] is None:
                r["why"] += " (warning: its factors make transforms slow)"
        else:
            r["why"] = f"over the limit of {cfg['RES_MAX']} resolutions"
    lines = []
    for r in sorted(ranked[:cfg["RES_MAX"]], key=lambda r: cells(r["res"])):
        (nx, nz), mx, my, nn = r["res"], r["mx"], r["my"], r["nodes"]
        both = cfg["ORIENT"] == "both" or (cfg["ORIENT"] == "elongated" and mx > 2 * my)
        lines += [(nz, nn, 1, mx, my, nx)] + ([(nz, nn, 1, my, mx, nx)] if both and mx != my else [])
        if cfg["SQUARE"]:
            sx, sy, _ = square_mesh(r["res"], nn)
            if (sx, sy) != (mx, my):
                lines.append((nz, nn, 1, sx, sy, nx))
    return rows, lines


def resolution_report(rows, lines, cfg) -> str:
    """Text explaining which resolutions will be benchmarked and why (``res.txt``).

    Parameters
    ----------
    rows, lines : list
        From :func:`plan_resolutions`.
    cfg : dict
        Configuration.

    Returns
    -------
    str
    """
    fx, fz = cfg["FLOOR"]
    if cfg["LINKED"]:
        span = (f"RES_Z from {fz} to {fz * (1 + cfg['RES_ABOVE']):.0f}, "
                f"RES_HOR = {cfg['RATIO']:g} x RES_Z")
    else:
        span = (f"RES_Z from {fz} to {fz * (1 + cfg['RES_ABOVE']):.0f}, "
                f"RES_HOR from {fx} to {fx * (1 + cfg['RES_HOR_ABOVE']):.0f}")
    out = [f"Stage 1 plan: {span}",
           f"Predicted cost = work per step (RES_HOR^2 RES_Z log) x steps (max(RES_HOR/Lx, RES_Z)), relative "
           f"to the floor {label(cfg['FLOOR'])}, x tier prior / mesh balance (priors, assumed: "
           + ", ".join(f"{t} {TIER_NAME[t]} {v:.2f}" for t, v in TIER_PRIOR.items()) + ").",
           "Nodes: the count chosen before at the nearest resolution on this cluster, scaled by size and "
           "rounded to even, not below the memory minimum; without one, the most balanced count from the "
           "memory minimum up to 1.35x it. Node counts from:"]
    out += [f"  {label(r['res'])}: {r['nodes_how']}" for r in rows if r.get("keep") and "nodes_how" in r]
    out += [f"Memory minimum (min): the fewest nodes that hold the projected peak (fixed {fixed_per_node(cfg):.0f} "
            f"GiB/node, budget {budget(cfg):.0f} GiB/node). Memory from:"]
    out += [f"  {h}" for h in sorted({r['mem_how'] for r in rows if 'mem_how' in r})]
    out += ["",
           f"{'RES_Z x RES_HOR':>16} {'factors (z | hor)':<24}{'tier':>5} {'min':>4} {'nodes':>6} {'mesh':>9} "
           f"{'bal':>6} {'pred':>6}  decision"]
    rejected = [r for r in rows if "pred" not in r]
    for r in rows:
        if "pred" in r:
            nx, nz = r["res"]
            out.append(f"{str(nz) + ' x ' + str(nx):>16} {factor_str(nz) + ' | ' + factor_str(nx):<24}"
                       f"{r['tier'] or '-':>5} {r['n0']:>4} {r['nodes']:>6} {str(r['mx']) + 'x' + str(r['my']):>9} "
                       f"{r['eff']:6.3f} {r['pred']:6.3f}  {'BENCHMARK' if r['keep'] else 'skip'}: {r['why']}")
    if rejected:
        out.append(f"Rejected (a prime of 11 or more, or 7^2): {len(rejected)} resolutions, e.g. "
                   + ", ".join(label(r["res"]) for r in rejected[:8]))
    out += ["", f"{len(lines)} jobs (orientations: {cfg['ORIENT']}; most square mesh added: "
                f"{'yes' if cfg['SQUARE'] else 'no'}):"]
    out += [f"  {label((nx, nz)):>9} on {b:>2} nodes, mesh {d}x{e}  balance {balance(d, e, (nx, nz)):.3f}"
            f"{'  (most square)' if is_square(d, e) else ''}" for nz, b, _, d, e, nx in lines]
    return "\n".join(out) + "\n"


def decide_resolution(rows, cfg):
    """Choose the resolution with the lowest node-hours per unit of simulated time.

    Cost per unit simulated time is ``node-seconds per step x steps per unit time``,
    with steps from :func:`steps`, reported relative to the floor. Only resolutions
    at or above the floor on both axes are eligible. If the runner-up is within the
    measured noise, the larger resolution is chosen: it is more accurate for the
    same cost.

    Parameters
    ----------
    rows : list of dict
        From :func:`load`.
    cfg : dict
        Configuration.

    Returns
    -------
    best : dict
        Chosen configuration.
    per_res : list of dict
        Best configuration of each resolution.
    text : str
        Contents of ``decision_resolution.txt``.
    """
    fx, fz = floor = cfg["FLOOR"]
    lx = cfg["LX"]
    by_res = {}
    for r in rows:
        if r["nz"] >= fz and r["nx"] >= fx and (r["res"] not in by_res or r["node_s"] < by_res[r["res"]]["node_s"]):
            by_res[r["res"]] = r
    if not by_res:
        raise SystemExit(f"no stage 1 result at or above the floor {label(floor)}")
    ref_res = floor if floor in by_res else min(by_res, key=cells)
    ref = by_res[ref_res]
    per_res = []
    for res, r in sorted(by_res.items(), key=lambda kv: cells(kv[0])):
        var, _, how = measured_variable(r, cfg)
        per_res.append(dict(r, cost=r["node_s"] * steps(res, lx) / (ref["node_s"] * steps(ref_res, lx)),
                            wall=r["s"] * steps(res, lx) / (ref["s"] * steps(ref_res, lx)),
                            work=(r["node_s"] / work(res)) / (ref["node_s"] / work(ref_res)),
                            var=var, how=how, n_min=min_nodes(var, cfg, r["c"]),
                            peak=job_peak(var, r["nodes"], cfg, r["c"])))
    ranked = sorted(per_res, key=lambda r: r["cost"])
    best, why = ranked[0], "lowest node-hours per unit of simulated time"
    if len(ranked) > 1:
        second = ranked[1]
        margin = second["cost"] / best["cost"] - 1
        tol = noise_of(best, second)
        if margin < tol and cells(second["res"]) > cells(best["res"]):
            best, why = second, (f"within noise ({margin:.1%} < {tol:.1%}) of {ranked[0]['label']}, "
                                 f"so the larger resolution is taken")

    out = [f"Resolution decision (floor RES_Z x RES_HOR = {fz} x {fx})", "",
           f"CHOICE: RES_Z = {best['nz']}, RES_HOR = {best['nx']} (best measured layout: {best['nodes']} nodes, "
           f"mesh {best['mesh']})",
           f"Reason: {why}.", "",
           f"{'RES_Z x RES_HOR':>16} {'nodes':>5} {'mesh':>8} {'bal':>6} {'s/step':>8} {'node-s':>7} "
           f"{'cost/t':>7} {'wall/t':>7} {'work':>6} {'noise':>6} {'GiB/node':>9} {'min nodes':>9}"]
    for r in per_res:
        noise = f"{r['noise']:.1%}" if r["noise"] is not None else "n/a"
        mem = f"{r['mem_io']:.0f}" if r["mem_io"] else "n/a"
        out.append(f"{str(r['nz']) + ' x ' + str(r['nx']):>16} {r['nodes']:>5} {r['mesh']:>8} {r['eff']:6.3f} "
                   f"{r['s']:8.4f} {r['node_s']:7.2f} {r['cost']:7.3f} {r['wall']:7.3f} {r['work']:6.3f} "
                   f"{noise:>6} {mem:>9} {r['n_min']:>9}{'  <-' if r['res'] == best['res'] else ''}")
    out += ["", "cost/t, wall/t = node-hours and wall time per unit simulated time, relative to the floor, with",
            "steps per unit time proportional to max(RES_HOR/Lx, RES_Z). work = node-seconds per RES_HOR^2 RES_Z log,",
            "relative to the floor (factorisation, balance and communication combined). bal = predicted mesh",
            "load balance. noise = spread of the timed segments. GiB/node = measured with production outputs",
            "(fullest node). min nodes = fewest nodes whose projected peak fits:",
            f"  {memory_rule(per_res, cfg)}.", ""]
    if len({r["long"] for r in per_res}) > 1:
        out += ["WARNING: these results mix the 10 minute benchmark with the earlier short one; their step",
                "times are not comparable (the 10 minute one includes output writes). Rerun the short ones.", ""]
    stopped = [r for r in per_res if r["mem_stop"]]
    if stopped:
        out += ["Stopped early because a node passed 95% of its memory (more nodes needed): "
                + ", ".join(f"{r['label']} on {r['nodes']} nodes" for r in stopped), ""]
    same = [r for r in rows if r["res"] == best["res"] and r is not by_res[best["res"]]]
    if same:
        out.append(f"Other configurations at {best['label']}:")
        out += [f"  {r['mesh']:>8} on {r['nodes']} nodes: {r['node_s'] / best['node_s']:.2f}x the node-seconds "
                f"per step of {best['mesh']}" for r in same]
        out.append("")
    below = sorted({r["label"] for r in rows if r["res"] not in by_res})
    if below:
        out += [f"Below the floor (not eligible): {', '.join(below)}", ""]
    hours = projection(best["s"], best["res"], cfg)
    if hours:
        out.append(f"At {best['nodes']} nodes, a job of SIM_TIME = {cfg['SIM_TIME']:g} would take about "
                   f"{hours:.1f} h ({estimate_dt(cfg)[2]}, scaled with the steps per unit time).")
    return best, per_res, "\n".join(out) + "\n"


# ---------------------------------------------------------------- stage 2

def node_floor(rows, res, cfg: dict):
    """Fewest nodes a production job at ``res`` fits on, from the stage 1 measurement.

    Memory only bounds the node counts from below; it does not decide where the
    sweep lies.

    Parameters
    ----------
    rows : list of dict
        Stage 1 results.
    res : tuple of int
        ``(RES_HOR, RES_Z)``.
    cfg : dict
        Configuration (``NODE_MIN`` raises the floor).

    Returns
    -------
    n_min : int
        Floor.
    var : float
        Resolution-dependent memory at a job's peak, GiB over all nodes.
    how : str
        Where the memory figure comes from.
    ref : dict
        The stage 1 row used (fewest node-seconds per step at ``res``).

    Raises
    ------
    SystemExit
        Without a stage 1 result at ``res``, or if the floor is above ``NODE_MAX``.
    """
    at = [r for r in rows if r["res"] == res]
    if not at:
        raise SystemExit(f"no stage 1 result at {label(res)}")
    ref = min(at, key=lambda r: r["node_s"])
    var, _, how = measured_variable(ref, cfg)
    n_min = min_nodes(var, cfg)
    if cfg["NODE_MIN"]:
        n_min = max(n_min, cfg["NODE_MIN"])
    if cfg["NODE_MAX"] and n_min > cfg["NODE_MAX"]:
        raise SystemExit(f"{label(res)} needs at least {n_min} nodes for memory, above NODE_MAX = {cfg['NODE_MAX']}")
    return n_min, var, how, ref


def mesh_rule(ref, res, cfg: dict):
    """How stage 2 picks the mesh for a node count, following the shape that won stage 1.

    The most square mesh if a square mesh beat a different, better-balanced one in
    stage 1; otherwise the most balanced. The orientation of the stage 1 winner is kept.

    Parameters
    ----------
    ref : dict
        Stage 1 row at ``res`` with the fewest node-seconds per step.
    res : tuple of int
        ``(RES_HOR, RES_Z)``.
    cfg : dict
        Configuration (``SQUARE``).

    Returns
    -------
    mesh : callable
        ``mesh(nodes) -> (meshx, meshy, balance)``.
    square : bool
        Whether the most square mesh is used.
    """
    flip = ref["mx"] < ref["my"]
    balanced = best_mesh(res, ref["nodes"])[:2]
    square = bool(cfg["SQUARE"]) and is_square(ref["mx"], ref["my"]) and \
        (max(ref["mx"], ref["my"]), min(ref["mx"], ref["my"])) != balanced
    pick = square_mesh if square else best_mesh

    def mesh(nn):
        mx, my, e = pick(res, nn)
        return (my, mx, e) if flip else (mx, my, e)
    return mesh, square


def plan_nodes(rows, res, cfg: dict):
    """Node counts to benchmark for one resolution, spanning ``NODE_SPAN`` around a seed.

    The seed is the count chosen before at the nearest resolution on this cluster,
    scaled by size (:func:`history_optimum`); without one, the stage 1 count. The sweep
    runs from seed / sqrt(NODE_SPAN) to seed x sqrt(NODE_SPAN); if the memory floor
    (:func:`node_floor`) cuts the bottom off, the top moves up so the sweep still spans
    a factor NODE_SPAN. ``NODE_MAX`` caps it. Meshes follow :func:`mesh_rule`; with the
    most balanced meshes, counts with balance below min(0.93, best - 0.03) are dropped.
    At most ``NODE_POINTS`` evenly spaced counts, plus the seed rounded to even. If the
    best count turns out to be at an edge of the sweep, :func:`edge_counts` adds more.

    Parameters
    ----------
    rows : list of dict
        Stage 1 results.
    res : tuple of int
        Chosen ``(RES_HOR, RES_Z)``.
    cfg : dict
        Configuration.

    Returns
    -------
    lines : list of tuple
        Configuration lines.
    text : str
        Explanation (``nodes.txt``).
    """
    n_min, var, how, ref = node_floor(rows, res, cfg)
    seed = history_optimum(res, cfg)
    centre, centre_how = seed if seed else (float(ref["nodes"]), f"the stage 1 count, {ref['nodes']}")
    root = math.sqrt(cfg["NODE_SPAN"])
    lo = max(n_min, math.floor(centre / root))
    hi = max(math.ceil(centre * root), math.ceil(lo * cfg["NODE_SPAN"]))
    if cfg["NODE_MAX"]:
        hi = min(hi, cfg["NODE_MAX"])
    mesh, square = mesh_rule(ref, res, cfg)
    options = [(nn, *mesh(nn)) for nn in range(lo, hi + 1)]
    top = max(o[3] for o in options)
    cut = 0.0 if square else min(top - 0.03, 0.93)
    kept = [o for o in options if o[3] >= cut]
    pts = int(cfg["NODE_POINTS"])
    if len(kept) > pts:
        idx = sorted({round(i * (len(kept) - 1) / (pts - 1)) for i in range(pts)})
        kept = [kept[i] for i in idx]
    target = min(max(even(centre, lo), lo), hi)
    near = min(options, key=lambda o: (abs(o[0] - target), -o[3]))
    if near not in kept:
        kept = sorted(kept + [near])
    nx, nz = res
    lines = [(nz, nn, 1, mx, my, nx) for nn, mx, my, _ in kept]
    out = [f"Stage 2 plan: node counts for RES_Z x RES_HOR = {nz} x {nx}", "",
           f"Seed: {centre_how}. Sweep {lo}-{hi}: a factor NODE_SPAN = {cfg['NODE_SPAN']:g} around the seed"
           + (f", capped at NODE_MAX = {cfg['NODE_MAX']}" if cfg["NODE_MAX"] and hi == cfg["NODE_MAX"] else "") + ".",
           f"Memory floor ({how}): {var:,.0f} GiB depends on the resolution at a job's peak; per node, budget "
           f"{budget(cfg):.0f} GiB minus {fixed_per_node(cfg):.0f} GiB fixed -> at least {n_min} nodes"
           + (f" (NODE_MIN = {cfg['NODE_MIN']})" if cfg["NODE_MIN"] else "")
           + (". The floor cut the sweep's bottom, so its top moved up." if lo == n_min and n_min > centre / root
              else "."),
           f"Meshes: {'most square' if square else 'most balanced'} for each node count, as {ref['mesh']} "
           f"won stage 1. " + ("Kept: evenly spaced counts" if square else f"Kept: balance >= {cut:.3f}")
           + f", plus {near[0]} (the seed rounded to even).",
           "If the best count is the largest (or smallest) of the sweep, later runs add counts beyond it.", "",
           f"{'nodes':>5} {'ranks':>6} {'mesh':>9} {'bal':>6}  kept"]
    for nn, mx, my, e in options:
        out.append(f"{nn:>5} {nn * ranks_per_node(1):>6} {str(mx) + 'x' + str(my):>9} {e:6.3f}  "
                   f"{'yes' if (nn, mx, my, e) in kept else ''}")
    return lines, "\n".join(out) + "\n"


def edge_counts(pts, best, cfg: dict):
    """Node counts to add when the best count is at an edge of those measured.

    A lowest J at the largest count measured is not yet a minimum: a larger count may
    be lower still. The minimum is bracketed once the best count has a measured count
    on each side (Press et al., Numerical Recipes, section 10.1). So when the best is
    the largest, :data:`EDGE_STEPS` multiples of it (rounded to even, capped at
    ``NODE_MAX``) are added; when it is the smallest and the memory floor allows, the
    same fractions below it. Nothing is added once the best is bracketed.

    Parameters
    ----------
    pts : list of dict
        From :func:`decide_nodes`, one per node count.
    best : dict
        Chosen point.
    cfg : dict
        Configuration.

    Returns
    -------
    counts : list of int
        Node counts not measured yet (may be empty).
    note : str
        One-line explanation.
    """
    counts = sorted(p["nodes"] for p in pts)
    n = best["nodes"]
    cap = cfg["NODE_MAX"] or float("inf")
    below = [p for p in pts if p["nodes"] < n]
    if any(p["fits"] for p in below) and n < counts[-1]:
        return [], f"bracketed: measured counts on both sides of {n} have a higher J."
    if n == counts[-1]:
        if n >= cap:
            return [], f"the best count, {n}, is the largest measured and NODE_MAX = {cfg['NODE_MAX']}; not extended."
        new = [m for m in sorted({min(even(n * f), cap) for f in EDGE_STEPS} - set(counts)) if m > n]
        return new, (f"the best count, {n}, is the largest measured, so larger counts are added: "
                     + ", ".join(map(str, new)) + "." if new else f"the best count, {n}, is the largest measured.")
    floor = max(min_nodes(max(p["var"] for p in pts), cfg, best["c"]), cfg["NODE_MIN"] or 1)
    if below or n <= floor:
        return [], f"the best count, {n}, is the smallest with memory for a full job; not extended."
    new = [m for m in sorted({max(even(n / f), floor) for f in EDGE_STEPS} - set(counts)) if m < n]
    return new, (f"the best count, {n}, is the smallest measured, so smaller counts are added: "
                 + ", ".join(map(str, new)) + "." if new else f"the best count, {n}, is the smallest measured.")


def edge_lines(counts, stage1_rows, res, cfg: dict):
    """Configuration lines for the counts from :func:`edge_counts`, meshes as in :func:`plan_nodes`.

    Parameters
    ----------
    counts : list of int
        Node counts.
    stage1_rows : list of dict
        Stage 1 results.
    res : tuple of int
        ``(RES_HOR, RES_Z)``.
    cfg : dict
        Configuration.

    Returns
    -------
    list of tuple
    """
    _, _, _, ref = node_floor(stage1_rows, res, cfg)
    mesh, _ = mesh_rule(ref, res, cfg)
    nx, nz = res
    return [(nz, nn, 1, *mesh(nn)[:2], nx) for nn in counts]


def decide_nodes(rows, res, cfg: dict):
    """Choose the node count that minimises node-hours x wall-time^k.

    For a node count N with step time T, node-hours per step are N T and wall time is
    T, so the objective is ``J = N T^(1+k)`` (``k = TIME_WEIGHT``). With ``r = N / N0``,
    speed-up ``S = T0 / T`` and relative efficiency ``E = S / r`` (``N0`` the smallest
    count with memory for a full job), ``J / J0 = r^-k E^-(1+k)``: a larger count beats
    the baseline exactly when ``E > r^(-k/(1+k))``, i.e. ``E > r^-1/2`` for k = 1 (0.82
    at 1.5x the nodes, 0.71 at 2x). Counts whose projected end-of-job memory exceeds
    the budget are not eligible. The count with the lowest J is chosen; the measured
    times are taken as exact (the noise is reported, but does not change the choice).

    Also reported: marginal efficiency against the next smaller count,
    ``(T_prev / T - 1) / (N / N_prev - 1)``, the share of ideal speed-up the added
    nodes delivered; and jumps, where a step to the next count is faster than ideal by
    more than the noise. The code cannot scale better than ideal, so a jump means the
    smaller counts were slowed by something that went away (memory pressure: a node
    without swap evicts and re-reads program code as it fills; or the network). The
    decision text lists the measured memory on both sides of each jump.

    Parameters
    ----------
    rows : list of dict
        Stage 2 results (and stage 1 results at the same resolution).
    res : tuple of int
        ``(RES_HOR, RES_Z)``.
    cfg : dict
        Configuration.

    Returns
    -------
    best : dict
        Chosen configuration.
    pts : list of dict
        One per node count.
    text : str
        Contents of ``decision_nodes.txt``.
    """
    k = cfg["TIME_WEIGHT"]
    by_nodes = {}
    for r in rows:
        if r["res"] == res and (r["nodes"] not in by_nodes or r["s"] < by_nodes[r["nodes"]]["s"]):
            by_nodes[r["nodes"]] = r
    if len(by_nodes) < 2:
        raise SystemExit(f"need results at two or more node counts for {label(res)}")
    pts = [dict(r) for _, r in sorted(by_nodes.items())]
    for p in pts:
        p["var"], p["var_now"], _ = measured_variable(p, cfg)
        p["peak"] = job_peak(p["var"], p["nodes"], cfg, p["c"])
        p["fits"] = p["peak"] <= budget(cfg) and not p["mem_stop"]
        p["bench_gib"] = max(p["mem_max"], p["mem_io"] or 0)
    base = next((p for p in pts if p["fits"]), pts[0])
    prev = None
    for p in pts:
        p["r"] = p["nodes"] / base["nodes"]
        p["S"] = base["s"] / p["s"]
        p["E"] = p["S"] / p["r"]
        p["thr"] = p["r"] ** (-k / (1 + k))
        p["Em"] = None if prev is None else (prev["s"] / p["s"] - 1) / (p["nodes"] / prev["nodes"] - 1)
        p["J"] = p["nodes"] * p["s"] ** (1 + k) / (base["nodes"] * base["s"] ** (1 + k))
        p["ok"] = p["E"] >= p["thr"]
        p["jump"] = prev is not None and (prev["s"] / p["s"]) > (p["nodes"] / prev["nodes"]) * (1 + 2 * noise_of(p, prev))
        prev = p
    eligible = [p for p in pts if p["fits"]] or pts
    top = min(eligible, key=lambda p: p["J"])
    best = top
    for p in pts:
        p["Jrel"] = p["J"] / top["J"]
    for p in pts:
        if p["mem_stop"]:
            p["why"] = "stopped: a node passed 95% of its memory"
        elif not p["fits"]:
            p["why"] = "too little memory for a full job"
        elif p is best:
            p["why"] = "CHOSEN"
        else:
            p["why"] = f"J {p['J'] / best['J'] - 1:+.0%} vs choice"

    nx, nz = res
    out = [f"Node-count decision for RES_Z x RES_HOR = {nz} x {nx}", "",
           f"CHOICE: {best['nodes']} nodes, mesh {best['mesh']}", ""]
    if best is base:
        out.append(f"Reason: no larger node count lowered node-hours x wall-time^{k:g}.")
    else:
        out.append(f"Reason: lowest node-hours x wall-time^{k:g}. Against {base['nodes']} nodes it runs "
                   f"{best['S']:.2f}x as fast ({1 - 1 / best['S']:.0%} less wall time) for "
                   f"{1 / best['E'] - 1:+.0%} node-hours (efficiency {best['E']:.2f}; break-even {best['thr']:.2f}).")
    out += ["", f"Rule: minimise J = nodes x T^{1 + k:g}; against the baseline this means E > r^({-k / (1 + k):.2f}). "
                "Counts without memory for a full job are excluded. Times are taken as measured: the noise",
                f"column is for information. Timed blocks slower than {STALL:g}x the median block (stalls) are left",
                "out of s/step; stall s = their total.", "",
            f"{'nodes':>5} {'mesh':>8} {'bal':>6} {'s/step':>8} {'node-s':>7} {'speedup':>8} {'E':>6} "
            f"{'E min':>6} {'marginal':>9} {'J/Jmin':>7} {'noise':>6} {'stall s':>8} {'GiB bench':>10} "
            f"{'GiB end':>8}  result"]
    for p in pts:
        em = f"{p['Em']:.2f}" if p["Em"] is not None else "-"
        noise = f"{p['noise']:.1%}" if p["noise"] is not None else "n/a"
        out.append(f"{p['nodes']:>5} {p['mesh']:>8} {p['eff']:6.3f} {p['s']:8.4f} {p['node_s']:7.2f} {p['S']:8.2f} "
                   f"{p['E']:6.2f} {p['thr']:6.2f} {em:>9} {p['Jrel']:7.3f} {noise:>6} {p['stall_s']:8.0f} "
                   f"{p['bench_gib']:10.0f} "
                   f"{p['peak']:8.0f}  {p['why']}")
    out += ["", f"speedup and E (efficiency) are relative to {base['nodes']} nodes, the smallest count with memory for",
            "a full job; E min = break-even against it; marginal = share of ideal speed-up delivered by the nodes",
            "added since the previous row (above 1: faster than ideal); J/Jmin = objective relative to the lowest;",
            "GiB bench = memory per node at the end of the benchmark (fullest node); GiB end = projected peak",
            f"per node of a production job: {memory_rule(pts, cfg)}."]
    if len({p["long"] for p in pts}) > 1:
        out += ["", "WARNING: these results mix the 10 minute benchmark with the earlier short one; their step",
                "times are not comparable (the 10 minute one includes output writes). Rerun the short ones."]
    short = [p for p in pts if p["truncated"]]
    if short:
        out += ["", "Shortened by the Slurm time limit (fewer timed steps; raise BENCH_TIME): "
                + ", ".join(f"{p['nodes']} nodes ({p['stepping_s'] / 60:.0f} min)" for p in short)]
    jumps = [(q, p) for q, p in zip(pts, pts[1:]) if p["jump"]]
    for q, p in jumps:
        out += ["", f"WARNING: {q['nodes']} -> {p['nodes']} nodes ran {q['s'] / p['s']:.2f}x as fast with "
                    f"{p['nodes'] / q['nodes']:.2f}x the nodes, faster than ideal beyond the noise.",
                "  The code cannot scale better than ideal (cache effects give a few percent, not this), so the",
                f"  counts up to {q['nodes']} were slowed by something that went away. Measured peak per node in the",
                f"  benchmark: {q['bench_gib']:.0f} GiB ({q['bench_gib'] / cfg['MEM_USABLE']:.0%} of MEM_USABLE) at "
                f"{q['nodes']} nodes, {p['bench_gib']:.0f} GiB ({p['bench_gib'] / cfg['MEM_USABLE']:.0%}) at {p['nodes']}.",
                "  - Near full (above ~90%): memory pressure. Without swap, a nearly full node evicts program code",
                "    and libraries and reads them again from the file system. Then even the benchmark had too few",
                "    nodes, and MEM_USABLE (or MEM_HEADROOM) should be lowered to below the slow counts' usage.",
                "  - Not near full: look in the Slurm output of the slow jobs for btl_tcp or 'Connection reset'",
                "    (MPI fell back to TCP) and rerun them; also check MEM_USABLE is right for this cluster",
                "    (scontrol show node <node> | grep RealMemory).",
                "  Until explained, the times at the slow counts are not trustworthy; the choice among the faster",
                "  counts stands."]
    if best["nodes"] * 192 < 1024 and "tri" in cluster_name(cfg).lower():
        out += ["", f"Note: {best['nodes']} nodes is {best['nodes'] * 192} cores, below Trillium's stated intent "
                    "of jobs of 1024 cores or more."]
    _, edge = edge_counts(pts, best, cfg)
    out += ["", f"Edge rule: {edge}"]
    others = [r for r in rows if r["res"] == res and r is not by_nodes[r["nodes"]]]
    if others:
        out += ["", "Other meshes measured at the same node counts (the faster mesh of each count is used above):"]
        for r in sorted(others, key=lambda r: (r["nodes"], r["mesh"])):
            q = by_nodes[r["nodes"]]
            out.append(f"  {r['nodes']:>3} nodes: {r['mesh']:>8} {r['s']:.4f} s/step = {r['s'] / q['s']:.3f}x "
                       f"{q['mesh']}; balance alone predicts {q['eff'] / r['eff']:.3f}x")
    est = estimate_dt(cfg)
    if est and cfg["SIM_TIME"]:
        out += ["", f"Per job of SIM_TIME = {cfg['SIM_TIME']:g} ({est[2]}, scaled with the steps per unit time):"]
        for p in pts:
            h = projection(p["s"], res, cfg)
            out.append(f"  {p['nodes']:>3} nodes: {h:6.1f} h wall, {h * p['nodes']:7.0f} node-hours")
    out += ["", "For production:",
            f"  params.conf:          export RES_Z={nz} RES_HOR={nx} MESHX={best['mx']} MESHY={best['my']}",
            f"  run_sim_trillium.sh:  #SBATCH --nodes={best['nodes']}"]
    return best, pts, "\n".join(out) + "\n"


# ---------------------------------------------------------------- jobs

def write_conf(lines, path: Path, header: str) -> None:
    """Write configuration lines (``RES_Z NODES C MESHX MESHY RES_HOR``).

    Parameters
    ----------
    lines : list of tuple
        Configurations.
    path : pathlib.Path
        Output file.
    header : str
        Comment for the top.
    """
    body = "".join(" ".join(map(str, line)) + "\n" for line in lines)
    path.write_text(f"# {header}\n# RES_Z NODES C MESHX MESHY RES_HOR\n" + body)


def tag(line) -> str:
    """Folder name of a configuration (``256_8n_c1_64x24``, or ``288_512_8n_c1_48x32``).

    Must match the folder name ``shell/bench_one.sh`` builds.

    Parameters
    ----------
    line : tuple
        ``(RES_Z, nodes, C, MESHX, MESHY, RES_HOR)``.

    Returns
    -------
    str
    """
    nz, nodes, c, mx, my, nx = line
    return f"{label((nx, nz))}_{nodes}n_c{c}_{mx}x{my}"


def node_hours(lines, cfg) -> float:
    """Upper estimate of the node-hours a list of jobs will use.

    Parameters
    ----------
    lines : list of tuple
        Configurations.
    cfg : dict
        Configuration.

    Returns
    -------
    float
    """
    return sum(line[1] for line in lines) * (cfg["BENCH_SECONDS"] + OVERHEAD_S) / 3600


def bench_time(cfg: dict) -> str:
    """Slurm ``--time`` for one benchmark job: ``BENCH_TIME``, or BENCH_SECONDS + :data:`SLACK_MIN` minutes.

    The slack covers the rehearsal launch, the two solver builds and the final writes.
    If it is too short, the job still records its result: stepping stops early enough
    before the allocation ends (``truncated`` in bench.csv).

    Parameters
    ----------
    cfg : dict
        Configuration.

    Returns
    -------
    str
        ``HH:MM:SS``.
    """
    if cfg["BENCH_TIME"]:
        return cfg["BENCH_TIME"]
    m = math.ceil(cfg["BENCH_SECONDS"] / 60) + SLACK_MIN
    return f"{m // 60:02d}:{m % 60:02d}:00"


def submit(lines, outdir: Path, cfg: dict, folder: Path) -> None:
    """Submit one job per configuration and record the job ids in ``outdir/jobs.tsv``.

    Parameters
    ----------
    lines : list of tuple
        Configurations.
    outdir : pathlib.Path
        Results folder.
    cfg : dict
        Configuration.
    folder : pathlib.Path
        Benchmark folder (the jobs' working directory).

    Raises
    ------
    SystemExit
        If a mesh does not fill its nodes or ``bench_one.sh`` is missing.
    """
    script = Path(cfg["SCRIPTS_PATH"] or "") / "shell" / "bench_one.sh"
    if not script.exists():
        raise SystemExit(f"{script} not found; is SCRIPTS_PATH set in params.conf?")
    outdir.mkdir(parents=True, exist_ok=True)
    with open(outdir / "jobs.tsv", "a") as log:
        for line in lines:
            nz, nodes, c, mx, my, nx = line
            rpn = ranks_per_node(c)
            if mx * my != nodes * rpn:
                raise SystemExit(f"{line}: {mx}x{my} = {mx * my} ranks, {nodes} nodes give {nodes * rpn}")
            (outdir / tag(line)).mkdir(exist_ok=True)
            args = ["sbatch", "--parsable", f"--nodes={nodes}", f"--ntasks-per-node={rpn}",
                    f"--cpus-per-task={c}", f"--time={bench_time(cfg)}", f"--account={cfg['ACCOUNT']}",
                    f"--job-name=bench_{tag(line)}", f"--output={outdir / tag(line)}/slurm_%j.out",
                    str(script), str(nz), str(nodes), str(c), str(mx), str(my), str(outdir.resolve()), str(nx)]
            jid = subprocess.run(args, cwd=folder, capture_output=True, text=True, check=True).stdout.strip()
            jid = jid.split(";")[0]
            log.write(f"{tag(line)}\t{jid}\n")
            print(f"  submitted {tag(line):<32} job {jid}")


def job_status(outdir: Path) -> list[dict]:
    """State of every submitted job of a stage.

    Parameters
    ----------
    outdir : pathlib.Path
        Results folder with ``jobs.tsv``.

    Returns
    -------
    list of dict
        ``tag``, ``job``, ``state`` (queue state, or ``done``/``failed``), ``result``.
    """
    jobs = [line.split("\t") for line in (outdir / "jobs.tsv").read_text().splitlines() if line.strip()]
    try:
        q = subprocess.run(["squeue", "-h", "-u", os.environ.get("USER", ""), "-o", "%i %T"],
                           capture_output=True, text=True, check=True).stdout.split("\n")
        queue = dict(line.split() for line in q if line.strip())
    except (OSError, subprocess.CalledProcessError):
        queue = {}
    out = []
    for t, j in jobs:
        has = (outdir / t / "bench.csv").exists()
        state = queue.get(j) or ("done" if has else "failed")
        out.append(dict(tag=t, job=j, state=state, result=has))
    return out


def print_status(name: str, st: list[dict], outdir: Path) -> bool:
    """Print a stage's jobs; True when none are queued or running.

    Parameters
    ----------
    name : str
        Stage name.
    st : list of dict
        From :func:`job_status`.
    outdir : pathlib.Path
        Results folder.

    Returns
    -------
    bool
    """
    print(f"{name}:")
    for s in st:
        extra = "" if s["state"] != "failed" else f"  (no result; see {outdir / s['tag']}/slurm_{s['job']}.out)"
        print(f"  {s['tag']:<32} job {s['job']:<10} {s['state']}{extra}")
    finished = all(s["state"] in ("done", "failed") for s in st)
    if not finished:
        print("  still running: run this again later (squeue -u $USER shows the queue)")
    return finished


def confirm(lines, cfg, yes: bool) -> bool:
    """Ask before submitting.

    Parameters
    ----------
    lines : list of tuple
        Configurations.
    cfg : dict
        Configuration.
    yes : bool
        Skip the question.

    Returns
    -------
    bool
    """
    nh = node_hours(lines, cfg)
    print(f"\n{len(lines)} jobs, each about {(cfg['BENCH_SECONDS'] + OVERHEAD_S) / 60:.0f} min (requesting "
          f"{bench_time(cfg)}); at most about {nh:.1f} node-hours ({nh * 192:,.0f} core-hours).")
    if yes:
        return True
    return input("Submit them? [y/N] ").strip().lower() in ("y", "yes")


# ---------------------------------------------------------------- plots

def colours(keys, expense) -> dict:
    """ROYGBIV colours by expense: red for the cheapest, violet for the dearest.

    Parameters
    ----------
    keys : iterable
        Items to colour.
    expense : callable
        Key -> expense, for ordering.

    Returns
    -------
    dict
        Key -> colour.
    """
    vals = sorted(set(keys), key=expense)
    n = len(vals)
    return {v: ROYGBIV[round(i * (len(ROYGBIV) - 1) / max(n - 1, 1))] for i, v in enumerate(vals)}


def plot_resolution(per_res, rows, best, path: Path, lx: float) -> None:
    """Stage 1 figure: cost and wall time per unit simulated time, work efficiency, memory.

    Parameters
    ----------
    per_res : list of dict
        From :func:`decide_resolution`.
    rows : list of dict
        All stage 1 rows.
    best : dict
        Chosen configuration.
    path : pathlib.Path
        Output name (``.pdf`` and ``.png`` are written).
    lx : float
        Box length (for the expense ordering of colours).
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"pdf.fonttype": 42})
    expense = lambda res: work(res) * steps(res, lx)
    col = colours((r["res"] for r in rows), expense)
    fig, ax = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    names = [r["label"] for r in per_res]
    for a, key, ylabel in ((ax[0, 0], "cost", "node-hours per unit sim time (rel. to floor)"),
                           (ax[0, 1], "wall", "wall time per unit sim time (rel. to floor)"),
                           (ax[1, 0], "work", "node-s per unit of work (rel. to floor)")):
        vals = [r[key] for r in per_res]
        err = [r[key] * (r["noise"] or 0) for r in per_res]
        bars = a.bar(names, vals, yerr=err, color=[col[r["res"]] for r in per_res], capsize=4)
        for b, r in zip(bars, per_res):
            b.set_edgecolor("black" if r["res"] == best["res"] else "none")
            b.set_linewidth(2)
            a.annotate(f"{r['mesh']}\n{r['nodes']} nodes", (b.get_x() + b.get_width() / 2, b.get_height()),
                       ha="center", va="bottom", fontsize=8, xytext=(0, 3), textcoords="offset points")
        a.axhline(1, color="grey", lw=0.8, ls="--")
        a.set_ylabel(ylabel)
        a.set_xlabel("resolution (RES_Z, or RES_Z_RES_HOR)")
    a = ax[1, 1]
    res_list = sorted({r["res"] for r in rows}, key=expense)
    for i, res in enumerate(res_list):
        group = [r for r in rows if r["res"] == res]
        for j, r in enumerate(group):
            x = i + (j - (len(group) - 1) / 2) * 0.8 / max(len(group), 1)
            a.vlines(x, r["mem_min"], r["mem_max"], color=col.get(res, "grey"), lw=4)
            top = r["mem_max"]
            if r["mem_io"]:
                a.plot(x, r["mem_io"], "D", color=col.get(res, "grey"), mec="black")
                top = max(top, r["mem_io"])
            a.annotate(f"{r['mesh']}\n{r['nodes']}n", (x, top), fontsize=7, ha="center",
                       va="bottom", xytext=(0, 4), textcoords="offset points")
    a.set_xticks(range(len(res_list)), [label(v) for v in res_list])
    a.set_xlim(-0.6, len(res_list) - 0.4)
    a.margins(y=0.15)
    a.set_xlabel("resolution")
    a.set_ylabel("peak GiB per node: bar = timing phase, ◆ = with outputs")
    for a in (ax[0, 0], ax[0, 1], ax[1, 0]):
        a.margins(y=0.18)
    fig.suptitle(f"Stage 1: resolution (choice {best['label']}, outlined)")
    save(fig, path)


def save(fig, path: Path) -> None:
    """Write a figure as a vector PDF (for papers and zooming) and a 200 dpi PNG.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
    path : pathlib.Path
        Output name; the suffix is replaced.
    """
    import matplotlib.pyplot as plt

    fig.savefig(path.with_suffix(".pdf"))
    fig.savefig(path.with_suffix(".png"), dpi=200)
    plt.close(fig)


def plot_nodes(pts, best, path: Path, cfg: dict) -> None:
    """Stage 2 figure, four panels, every comparison made against the choice.

    (a) Step time against nodes (log-log), with perfect scaling and the break-even
        curve through the choice: a point below the break-even curve has lower J.
    (b) Memory per node, projected to the end of a production job, split into the
        fixed part, the resolution-dependent part at the start and its growth.
    (c) J relative to its minimum (lower is better), with
        the marginal efficiency of each step up in node count.
    (d) Cost against wall time, each relative to the best feasible value, with curves
        of equal J: the time-to-solution / cost Pareto plot of HPC scaling studies.

    Parameters
    ----------
    pts : list of dict
        From :func:`decide_nodes`.
    best : dict
        Chosen configuration.
    path : pathlib.Path
        Output name (``.pdf`` and ``.png`` are written).
    cfg : dict
        Configuration.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter

    plt.rcParams.update({"pdf.fonttype": 42, "font.size": 9})
    blue, red, grey, green, light = "#1f77b4", "#d62728", "#7f7f7f", "#2ca02c", "#aec7e8"
    k = cfg["TIME_WEIGHT"]
    n = np.array([p["nodes"] for p in pts], float)
    s = np.array([p["s"] for p in pts])
    noise = np.array([p["noise"] or 0 for p in pts])
    nb, sb = best["nodes"], best["s"]
    plain = FuncFormatter(lambda v, _: f"{v:g}")
    wt = "wall-time" if k == 1 else f"wall-time^{k:g}"
    idx = np.arange(len(pts))

    def colour(p):
        return "black" if p is best else (blue if p["fits"] else red)

    def ticks(axis, lo, hi):
        """Plain-number ticks on a log axis, denser when it spans less than a factor of 3."""
        dense = (1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.7, 2, 2.5, 3, 4, 5, 6, 7, 8, 9)
        axis.set_major_locator(LogLocator(subs=dense if hi / lo < 3 else (1, 1.5, 2, 3, 5, 7)))
        axis.set_major_formatter(plain)
        axis.set_minor_formatter(NullFormatter())

    def logy(a, lo, hi, room=1.15):
        """Log y axis showing every value in [lo, hi] (``room``: space above) with plain labels."""
        a.set_yscale("log")
        a.set_ylim(lo / 1.15, hi * room)
        ticks(a.yaxis, lo / 1.15, hi * room)

    tight = not all(p["fits"] for p in pts)
    no_mem = "too little memory for a job"

    fig, ax = plt.subplots(2, 2, figsize=(12.5, 9.5), constrained_layout=True)

    # (a) step time
    a = ax[0, 0]
    nn = np.geomspace(n.min(), n.max(), 200)
    ideal, even = sb * nb / nn, sb * (nb / nn) ** (1 / (1 + k))
    a.plot(nn, ideal, "--", color=grey, lw=1)
    a.plot(nn, even, "-", color=green, lw=1.2)
    for p, x, y, e in zip(pts, n, s, noise):
        a.errorbar(x, y, yerr=y * e, fmt="o", ms=7 if p is best else 5, capsize=3, color=colour(p),
                   mfc=colour(p) if p["fits"] else "white")
    for q, p in zip(pts, pts[1:]):
        if p["jump"]:
            a.annotate(f"{q['s'] / p['s']:.2f}x faster with {p['nodes'] / q['nodes']:.2f}x nodes:\n"
                       "faster than ideal, see the text file", xy=(p["nodes"], p["s"]),
                       xytext=((q["nodes"] * p["nodes"]) ** 0.5, (q["s"] * p["s"]) ** 0.5),
                       color=red, fontsize=7.5, ha="right", va="center",
                       bbox=dict(fc="white", ec="none", alpha=0.85, pad=1),
                       arrowprops=dict(arrowstyle="->", color=red, lw=0.8))
    a.set_xscale("log")
    a.set_xticks(n)
    a.xaxis.set_major_formatter(plain)
    a.xaxis.set_minor_formatter(NullFormatter())
    a.set_xlim(n.min() / 1.08, n.max() * 1.08)
    logy(a, min((s * (1 - noise)).min(), ideal.min(), even.min()),
         max((s * (1 + noise)).max(), ideal.max(), even.max()), room=2.3)
    a.set_xlabel("nodes")
    a.set_ylabel("wall time per step (s)")
    a.set_title("(a) Step time (log-log)")
    a.legend(handles=[Line2D([], [], ls="--", color=grey, label="perfect scaling through the choice (T ∝ 1/N)"),
                      Line2D([], [], color=green, label=f"equal J to the choice (T ∝ N^-{1 / (1 + k):.2g}): "
                                                         "points below it are better"),
                      Line2D([], [], ls="", marker="o", color="black", label=f"choice ({nb} nodes)"),
                      Line2D([], [], ls="", marker="o", color=blue, label="measured (bars: noise)"),
                      ] + [Line2D([], [], ls="", marker="o", color=red, mfc="white", label=no_mem)] * tight,
             fontsize=7.5, loc="upper right")

    # (b) memory per node
    a = ax[0, 1]
    fixed = np.array([fixed_per_node(cfg, p["c"]) for p in pts])
    start = np.array([p["var_now"] / p["nodes"] for p in pts])
    grow = np.array([(p["var"] - p["var_now"]) / p["nodes"] for p in pts])
    a.bar(idx, fixed, color="#c7c7c7", edgecolor="white")
    a.bar(idx, start, bottom=fixed, color=light, edgecolor="white")
    a.bar(idx, grow, bottom=fixed + start, color="white", edgecolor=blue, hatch="///", lw=0.6)
    a.plot(idx, [p["bench_gib"] for p in pts], "D", color="black", ms=5, zorder=5)
    a.axhline(budget(cfg), color=red, ls="--", lw=1)
    a.axhline(cfg["MEM_USABLE"], color="black", ls=":", lw=1)
    for x, p in zip(idx, pts):
        a.annotate(f"{p['peak']:.0f}" + ("" if p["fits"] else "\ntoo full"), (x, p["peak"]), ha="center",
                   va="bottom", fontsize=7.5, color=red if not p["fits"] else "black",
                   xytext=(0, 2), textcoords="offset points", bbox=dict(fc="white", ec="none", pad=0.5))
    a.set_xticks(idx, [str(p["nodes"]) for p in pts])
    for t, p in zip(a.get_xticklabels(), pts):
        t.set_fontweight("bold" if p is best else "normal")
    a.set_ylim(0, max(cfg["MEM_USABLE"], max(p["peak"] for p in pts), max(p["bench_gib"] for p in pts)) * 1.35)
    a.set_xlabel("nodes")
    a.set_ylabel("GiB per node")
    a.set_title(f"(b) Memory per node at the peak of a {cfg['JOB_HOURS']:g} h job (projected)\n"
                "the same total spread over more nodes, so less per node", fontsize=9)
    a.legend(handles=[Patch(color="#c7c7c7", label="fixed (libraries, MPI)"),
                      Patch(color=light, label="fields and outputs, end of the benchmark"),
                      Patch(facecolor="white", edgecolor=blue, hatch="///",
                            label=f"projected rise to the peak (÷ {cfg['BENCH_MEM_FRACTION']:g})"
                            if all(p["long"] for p in pts) else f"projected rise over {cfg['JOB_HOURS']:g} h"),
                      Line2D([], [], ls="", marker="D", color="black", label="measured (fullest node)"),
                      Line2D([], [], ls="--", color=red, label=f"budget ({cfg['MEM_HEADROOM']:.0%} of node)"),
                      Line2D([], [], ls=":", color="black", label=f"node memory ({cfg['MEM_USABLE']:.0f} GiB)")],
             fontsize=7.5, loc="upper right", ncol=2)

    # (c) objective
    a = ax[1, 0]
    jr = np.array([p["Jrel"] for p in pts])
    for x, p in zip(idx, pts):
        a.bar(x, p["Jrel"], color=colour(p) if p["fits"] else "white", edgecolor=colour(p),
              hatch=None if p["fits"] else "xx", lw=1, width=0.7)
        dn, dt = p["nodes"] * p["s"] / (nb * sb) - 1, p["s"] / sb - 1
        txt = "choice" if p is best else f"{dn:+.0%} node-h\n{dt:+.0%} wall"
        a.annotate(txt, (x, p["Jrel"]), ha="center", va="bottom", fontsize=7, xytext=(0, 2),
                   textcoords="offset points")
    logy(a, 1, jr.max(), room=1.9)
    a.set_xticks(idx, [str(p["nodes"]) if p["Em"] is None else f"{p['nodes']}\n{p['Em']:.0%}" for p in pts])
    for t, p in zip(a.get_xticklabels(), pts):
        t.set_color(red if p["jump"] else "black")
        t.set_fontweight("bold" if p is best else "normal")
    a.set_xlim(-0.6, len(pts) - 0.4)
    a.set_xlabel("nodes\n(below: share of perfect speed-up gained from the previous count; red: above 100%)")
    a.set_ylabel(f"J / lowest J   (J = node-hours x {wt})")
    a.set_title("(c) Objective, lower is better; labels: change against the choice", fontsize=9)
    if tight:
        a.legend(handles=[Patch(facecolor="white", edgecolor=red, hatch="xx", label=no_mem)],
                 fontsize=7.5, loc="upper right")

    # (d) cost against speed
    a = ax[1, 1]
    feas = [p for p in pts if p["fits"]] or pts
    cmin = min(p["nodes"] * p["s"] for p in feas)
    tmin = min(p["s"] for p in feas)
    x = np.array([p["nodes"] * p["s"] / cmin for p in pts])
    y = s / tmin
    front = []
    for p, xi, yi in sorted(zip(pts, x, y), key=lambda t: (t[1], t[2])):
        if p["fits"] and (not front or yi < front[-1][2]):
            front.append((p, xi, yi))
    a.plot([f[1] for f in front], [f[2] for f in front], "-", color=blue, lw=1, zorder=1)
    on_front = {id(f[0]) for f in front}
    for p, xi, yi in zip(pts, x, y):
        filled = p["fits"] and id(p) in on_front
        a.plot(xi, yi, "o" if p["fits"] else "X", ms=8 if p is best else 6, color=colour(p),
               mfc=colour(p) if filled else "white", zorder=4)
        a.annotate(f"{p['nodes']}", (xi, yi), xytext=(5, 3), textcoords="offset points", fontsize=8,
                   fontweight="bold" if p is best else "normal", zorder=5)
    a.set_xscale("log")
    a.set_yscale("log")
    xl, yl = (x.min() / 1.12, x.max() * 1.25), (y.min() / 1.12, y.max() * 1.9)
    a.set_xlim(*xl)
    a.set_ylim(*yl)
    ticks(a.xaxis, *xl)
    ticks(a.yaxis, *yl)
    xb, yb = best["nodes"] * sb / cmin, sb / tmin
    xx = np.geomspace(*xl, 200)
    for m in (1, 1.25, 1.5, 2, 3):
        if k > 0:
            yy = yb * (m * xb / xx) ** (1 / k)
            keep = (yy >= yl[0]) & (yy <= yl[1])
            if keep.sum() < 2:
                continue
            a.plot(xx[keep], yy[keep], "-", color=green, lw=0.8 if m > 1 else 1.4, alpha=0.8, zorder=0)
            if m > 1:
                mid = np.flatnonzero(keep)[keep.sum() // 2]
                a.annotate(f"J {m:g}x choice", (xx[mid], yy[mid]), fontsize=7, color=green, ha="center",
                           va="center", bbox=dict(fc="white", ec="none", pad=0.3), zorder=2)
        elif xl[0] < m * xb < xl[1]:
            a.axvline(m * xb, color=green, lw=0.8)
    a.set_autoscale_on(False)
    a.legend(handles=[Line2D([], [], color=green, lw=1.4, label="equal J to the choice; thin: 1.25x to 3x"),
                      Line2D([], [], ls="", marker="o", color=blue, label="Pareto front (filled)"),
                      Line2D([], [], ls="", marker="o", color=blue, mfc="white",
                             label="dominated: another count is cheaper and faster"),
                      ] + [Line2D([], [], ls="", marker="X", color=red, mfc="white", label=no_mem)] * tight,
             fontsize=7.5, loc="upper left")
    a.set_xlabel("node-hours per step ÷ cheapest     (← cheaper)")
    a.set_ylabel("wall time per step ÷ fastest     (↓ faster)")
    a.set_title("(d) Cost against speed, each relative to the best feasible count; lower left is better",
                fontsize=9)

    nx, nz = best["res"]
    fig.suptitle(f"Stage 2: node count for RES_Z x RES_HOR = {nz} x {nx}.  "
                 f"Choice: {nb} nodes, mesh {best['mesh']}  (minimises J = node-hours x {wt})",
                 fontsize=11)
    save(fig, path)


# ---------------------------------------------------------------- the workflow

def chosen_resolution(res_dir: Path):
    """Resolution chosen in stage 1, from its decision file.

    Parameters
    ----------
    res_dir : pathlib.Path
        Stage 1 results folder.

    Returns
    -------
    tuple of int or None
        ``(RES_HOR, RES_Z)``; decision files without RES_HOR mean 2 RES_Z.
    """
    f = res_dir / "decision_resolution.txt"
    if not f.exists():
        return None
    m = re.search(r"CHOICE: RES_Z = (\d+)(?:, RES_HOR = (\d+))?", f.read_text())
    if not m:
        return None
    nz = int(m.group(1))
    return int(m.group(2)) if m.group(2) else 2 * nz, nz


def do_decide_res(res_dir: Path, cfg: dict):
    """Write the stage 1 decision; return the chosen resolution.

    Parameters
    ----------
    res_dir : pathlib.Path
        Stage 1 results folder.
    cfg : dict
        Configuration.

    Returns
    -------
    tuple of int
        ``(RES_HOR, RES_Z)``.
    """
    rows = load(res_dir)
    best, per_res, text = decide_resolution(rows, cfg)
    n = record_history(rows, cfg, res_dir.parent.parent)
    if n:
        text += f"\n{n} results recorded in {history_path(cfg)} (cluster {cluster_name(cfg)}).\n"
    (res_dir / "decision_resolution.txt").write_text(text)
    plot_resolution(per_res, rows, best, res_dir / "decision_resolution", cfg["LX"])
    print(text + f"(written to {res_dir}/decision_resolution.txt, .pdf and .png)\n")
    return best["res"]


def do_decide_nodes(nodes_dir: Path, res, cfg: dict) -> None:
    """Write the stage 2 decision, using the stage 1 results at the same resolution too.

    Parameters
    ----------
    nodes_dir : pathlib.Path
        Stage 2 results folder.
    res : tuple of int
        ``(RES_HOR, RES_Z)``.
    cfg : dict
        Configuration.

    Returns
    -------
    best : dict
        Chosen configuration.
    pts : list of dict
        One per node count (see :func:`decide_nodes`).
    """
    reused = [r for r in load(nodes_dir.parent / "res") if r["res"] == res]
    rows = load(nodes_dir) + reused
    best, pts, text = decide_nodes(rows, res, cfg)
    n = record_history(rows, cfg, nodes_dir.parent.parent, chosen=best)
    if n:
        text += f"\n{n} results recorded in {history_path(cfg)} (cluster {cluster_name(cfg)}), the choice marked.\n"
    (nodes_dir / "decision_nodes.txt").write_text(text)
    plot_nodes(pts, best, nodes_dir / "decision_nodes", cfg)
    print(text + f"(written to {nodes_dir}/decision_nodes.txt, .pdf and .png)")
    return best, pts


def step(folder: Path, yes: bool = False, dry: bool = False) -> None:
    """Do whatever is due next (see the module docstring).

    Parameters
    ----------
    folder : pathlib.Path
        Benchmark folder.
    yes : bool, optional
        Submit without asking.
    dry : bool, optional
        Show the next plan but submit nothing.
    """
    cfg = load_config(folder)
    res_dir, nodes_dir = folder / "bench" / "res", folder / "bench" / "nodes"

    if not (res_dir / "jobs.tsv").exists():
        rows, lines = plan_resolutions(cfg)
        text = resolution_report(rows, lines, cfg)
        print(text)
        if dry:
            return
        (folder / "res.txt").write_text(text)
        write_conf(lines, folder / "res.conf", "stage 1 plan (see res.txt)")
        if confirm(lines, cfg, yes):
            submit(lines, res_dir, cfg, folder)
            print("\nStage 1 submitted. Run 'python -m rbc.bench' again to check on it.")
        return

    if not print_status("Stage 1 (resolution)", job_status(res_dir), res_dir):
        return
    if not (nodes_dir / "jobs.tsv").exists():
        done = {line.split("\t")[0] for line in (res_dir / "jobs.tsv").read_text().splitlines() if line.strip()}
        _, planned = plan_resolutions(cfg)
        missing = [line for line in planned if tag(line) not in done]
        if missing:
            print("\nStage 1 plan has configurations not yet run (settings changed since it was submitted):")
            print("\n".join(f"  {tag(line)}" for line in missing))
            if not dry and confirm(missing, cfg, yes):
                submit(missing, res_dir, cfg, folder)
                print("\nSubmitted. Run 'python -m rbc.bench' again to check on them.")
            return
        print()
        res = do_decide_res(res_dir, cfg)
        lines, text = plan_nodes(load(res_dir), res, cfg)
        have = {tag(line_of(r)) for r in load(res_dir)}
        reuse = [line for line in lines if tag(line) in have]
        lines = [line for line in lines if tag(line) not in have]
        if reuse:
            text += "\nAlready measured in stage 1, reused: " + ", ".join(tag(line) for line in reuse) + "\n"
        print(text)
        if dry:
            return
        (folder / "nodes.txt").write_text(text)
        write_conf(lines, folder / "nodes.conf", "stage 2 plan (see nodes.txt)")
        if not lines:
            nodes_dir.mkdir(parents=True, exist_ok=True)
            (nodes_dir / "jobs.tsv").touch()
            print("Nothing new to run. Run 'python -m rbc.bench' again for the decision.")
        elif confirm(lines, cfg, yes):
            submit(lines, nodes_dir, cfg, folder)
            print("\nStage 2 submitted. Run 'python -m rbc.bench' again to check on it.")
        return

    print()
    if not print_status("Stage 2 (nodes)", job_status(nodes_dir), nodes_dir):
        return
    print()
    res = chosen_resolution(res_dir)
    best, pts = do_decide_nodes(nodes_dir, res, cfg)
    counts, note = edge_counts(pts, best, cfg)
    submitted = {line.split("\t")[0] for d in (res_dir, nodes_dir)
                 for line in (d / "jobs.tsv").read_text().splitlines() if line.strip()}
    lines = [line for line in edge_lines(counts, load(res_dir), res, cfg) if tag(line) not in submitted]
    if counts and not lines:
        print("\nEdge rule: those counts were submitted before without a result (see their slurm_*.out); "
              "not resubmitted.")
    if not lines or dry:
        return
    print("\nEdge rule, jobs to add:\n" + "\n".join(f"  {tag(line)}" for line in lines))
    if confirm(lines, cfg, yes):
        submit(lines, nodes_dir, cfg, folder)
        print("\nSubmitted. Run 'python -m rbc.bench' again for the decision with them.")


def main(argv=None) -> None:
    """Command line.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments (default: ``sys.argv[1:]``).
    """
    p = argparse.ArgumentParser(prog="python -m rbc.bench", description=__doc__.split("\n")[0])
    p.add_argument("cmd", nargs="?", default="run", choices=("run", "status", "plan", "decide"),
                   help="run (default): do the next step; status; plan: show the next plan only; "
                        "decide: rewrite the decision files")
    p.add_argument("--dir", type=Path, default=Path.cwd(), help="benchmark folder (default: here)")
    p.add_argument("--yes", action="store_true", help="submit without asking")
    a = p.parse_args(argv)
    folder = a.dir.resolve()
    res_dir, nodes_dir = folder / "bench" / "res", folder / "bench" / "nodes"
    if a.cmd == "run":
        step(folder, yes=a.yes)
    elif a.cmd == "plan":
        step(folder, dry=True)
    elif a.cmd == "status":
        for name, d in (("Stage 1 (resolution)", res_dir), ("Stage 2 (nodes)", nodes_dir)):
            if (d / "jobs.tsv").exists():
                print_status(name, job_status(d), d)
            else:
                print(f"{name}: not submitted")
    elif a.cmd == "decide":
        cfg = load_config(folder)
        if not any(res_dir.rglob("bench.csv")):
            raise SystemExit("no stage 1 results yet")
        res = do_decide_res(res_dir, cfg)
        if (nodes_dir / "jobs.tsv").exists():
            do_decide_nodes(nodes_dir, chosen_resolution(res_dir) or res, cfg)


if __name__ == "__main__":
    main()
