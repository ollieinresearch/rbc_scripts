"""Two-stage benchmark driven by bench.conf: choose the resolution, then the node count.

Run from a benchmark folder holding ``params.conf`` (physics, output settings and
paths, as for a run) and ``bench.conf`` (what to benchmark; template in
``shell/bench.conf``), on a login node::

    python -m rbc.bench            # do the next step (asks before submitting anything)
    python -m rbc.bench status     # progress of the submitted jobs
    python -m rbc.bench plan       # show what the next step would submit, submit nothing
    python -m rbc.bench decide     # rewrite the decision files from the results so far

Run it repeatedly; each call does what is due:

1. plan stage 1 (resolutions at or above ``RES_FLOOR``) and submit one job per
   configuration, each on exactly the nodes its mesh fills (``shell/bench_one.sh``);
2. while jobs run, report their progress;
3. when stage 1 is done, write ``bench/res/decision_resolution.txt`` and ``.png``,
   plan stage 2 (node counts for the chosen resolution, from its measured memory) and
   submit it;
4. when stage 2 is done, write ``bench/nodes/decision_nodes.txt`` and ``.png``, which
   end with the ``MESHX``/``MESHY``/``--nodes`` lines for production.

The rules behind each choice are written into the decision files; ``docs/benchmarking.md``
derives them.
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import re
import subprocess
from pathlib import Path

from rbc.plan import efficiency, groups, meshes, ranks_per_node

NOISE = 0.03
OVERHEAD_S = 240
TIER_PRIOR = {"A": 1.00, "B": 1.03, "C": 1.08}
TIER_NAME = {"A": "3-smooth", "B": "5-smooth", "C": "one factor of 7"}
ROYGBIV = ("#d62728", "#ff7f0e", "#e3b505", "#2ca02c", "#1f77b4", "#3f007d", "#9467bd")

DEFAULTS = dict(RES_FLOOR=None, RES_ABOVE=0.125, RES_BAND=0.15, RES_MAX=4, ORIENT="both",
                BENCH_SECONDS=240.0, BENCH_IO=30, BENCH_REHEARSAL=10, BENCH_TIME="00:30:00",
                ACCOUNT="rrg-goluskin", MEM_FIXED=90.0, MEM_PER_CELL=1.69e-4, MEM_GROWTH=0.015,
                JOB_HOURS=12.0, MEM_USABLE=747.0, MEM_HEADROOM=0.85, NODE_SPAN=2.0, NODE_POINTS=6,
                NODE_MAX=None, TIME_WEIGHT=1.0, DT_REF=None, DT_REF_RES=None)
PARAMS = ("RES_Z", "RA", "PR", "SIM_TIME", "RESTART_FROM", "SCRIPTS_PATH", "USER")


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
    for key in ("RES_FLOOR", "RES_MAX", "BENCH_IO", "BENCH_REHEARSAL", "NODE_POINTS", "NODE_MAX", "DT_REF_RES"):
        if cfg[key] is not None:
            cfg[key] = int(float(cfg[key]))
    cfg["SIM_TIME"] = float(cfg["SIM_TIME"]) if cfg["SIM_TIME"] else None
    return cfg


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


def tier(res: int) -> str | None:
    """Transform-size class of a resolution.

    The Fourier grid (3 RES_Z for gamma = 2) and Chebyshev grid (1.5 RES_Z) contain the
    primes of RES_Z plus a 3, so classifying RES_Z classifies every transform.

    Parameters
    ----------
    res : int
        ``RES_Z`` (must be even so that 1.5 RES_Z is an integer).

    Returns
    -------
    str or None
        ``"A"`` 3-smooth, ``"B"`` 5-smooth, ``"C"`` 7-smooth with a single 7, or None
        (odd, a prime of 11 or more, or 7^2).
    """
    if res % 2:
        return None
    f = factorise(res)
    if any(p > 7 for p in f) or f.get(7, 0) > 1:
        return None
    if 7 in f:
        return "C"
    return "B" if 5 in f else "A"


def best_mesh(res: int, nodes: int, c: int = 1):
    """Most balanced mesh for a resolution on a node count.

    Parameters
    ----------
    res : int
        ``RES_Z``.
    nodes : int
        Node count.
    c : int, optional
        cpus-per-task.

    Returns
    -------
    tuple
        ``(meshx, meshy, efficiency)`` with ``meshx >= meshy``.
    """
    n = nodes * ranks_per_node(c)
    p0, p1, e = max(meshes(n, groups(res)), key=lambda m: (round(m[2], 4), -abs(m[0] - m[1])))
    return max(p0, p1), min(p0, p1), e


def cost_model(res: int, ref: int) -> float:
    """Work to reach a fixed simulation time, relative to resolution ``ref``.

    Work per step scales as RES_Z^3 log RES_Z at fixed aspect ratio; the advective CFL
    limit makes dt scale as 1/RES_Z, so the number of steps scales as RES_Z.

    Parameters
    ----------
    res, ref : int
        Resolutions.

    Returns
    -------
    float
        ``(res/ref)^4 log(res) / log(ref)``.
    """
    return (res / ref) ** 4 * math.log(res) / math.log(ref)


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
    """Factor by which the resolution-dependent memory grows over one production job.

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
    """Fewest nodes that hold ``variable`` GiB (grown over a job) plus the fixed part.

    Parameters
    ----------
    variable : float
        Resolution-dependent memory over all nodes at the start of a job, GiB.
    cfg : dict
        Configuration.
    c : int, optional
        cpus-per-task.

    Returns
    -------
    int
    """
    return max(1, math.ceil(variable * growth(cfg) / (budget(cfg) - fixed_per_node(cfg, c))))


def prior_variable(res: int, cfg: dict) -> float:
    """Resolution-dependent production memory before any measurement, GiB.

    Parameters
    ----------
    res : int
        RES_Z.
    cfg : dict
        Configuration (``MEM_PER_CELL`` GiB per RES_Z^3).

    Returns
    -------
    float
    """
    return cfg["MEM_PER_CELL"] * res ** 3


def measured_variable(row: dict, cfg: dict) -> tuple[float, str]:
    """Resolution-dependent production memory from a benchmark row, GiB.

    Uses the peak per node with production outputs (the larger of summed resident
    memory and the node's MemTotal - MemAvailable), times the node count, minus the
    fixed part. Rows without the outputs phase fall back to :func:`prior_variable`.

    Parameters
    ----------
    row : dict
        From :func:`load`.
    cfg : dict
        Configuration.

    Returns
    -------
    variable : float
        GiB over all nodes.
    how : str
        Where the number came from.
    """
    if row["mem_io"] is None:
        return prior_variable(row["nz"], cfg), "prior (no outputs phase in this result)"
    v = row["mem_io"] * row["nodes"] - fixed_per_node(cfg, row["c"]) * row["nodes"]
    return v, f"measured with outputs: {row['mem_io']:.0f} GiB/node on {row['nodes']} nodes"


def job_peak(variable: float, nodes: int, cfg: dict, c: int = 1) -> float:
    """Projected peak memory per node at the end of a production job, GiB.

    Parameters
    ----------
    variable : float
        Resolution-dependent memory over all nodes at the start, GiB.
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
    return fixed_per_node(cfg, c) + variable * growth(cfg) / nodes


# ---------------------------------------------------------------- timestep for projections

def estimate_dt(cfg: dict):
    """Mean production timestep at a known resolution, to scale as 1/RES_Z.

    From ``DT_REF`` (at ``DT_REF_RES``, default ``RES_FLOOR``) in bench.conf, otherwise
    the mean of the last analysis writes of the run in ``RESTART_FROM``.

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
        return cfg["DT_REF"], cfg["DT_REF_RES"] or cfg["RES_FLOOR"], "DT_REF in bench.conf"
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
        res = int(source(Path(src), ("params.conf",))["RES_Z"])
        return float(np.mean(ts)), res, f"mean of the last {len(ts)} timesteps of {src}"
    except Exception:
        return None


def projection(s_per_iter: float, res: int, cfg: dict):
    """Wall time per production job (``SIM_TIME`` time units), hours.

    Parameters
    ----------
    s_per_iter : float
        Seconds per step.
    res : int
        RES_Z.
    cfg : dict
        Configuration.

    Returns
    -------
    float or None
    """
    est = estimate_dt(cfg)
    if est is None or not cfg["SIM_TIME"]:
        return None
    dt = est[0] * est[1] / res
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
        One per configuration, with numeric fields converted.
    """
    path = Path(path)
    files = [path] if path.is_file() else sorted(path.rglob("bench.csv"))
    raw = []
    for f in files:
        with open(f, newline="") as h:
            raw += list(csv.DictReader(h))
    merged = {}
    for r in raw:
        if str(r.get("finite", "True")) != "True":
            continue
        key = (int(r["nz"]), int(r["nodes"]), int(r["cpus_per_rank"]), r["mesh"])
        merged.setdefault(key, []).append(r)

    def num(r, k):
        v = r.get(k)
        return float(v) if v not in (None, "", "nan") else None

    rows = []
    for (nz, nodes, c, mesh), rs in merged.items():
        s = [float(r["s_per_iter"]) for r in rs]
        mean = sum(s) / len(s)
        spreads = [float(r["s_per_iter_spread"]) for r in rs if r.get("s_per_iter_spread")]
        noise = max(spreads + [(max(s) - min(s)) / mean if len(s) > 1 else 0.0])
        io = [v for v in (num(rs[-1], "node_rss_io_max_gib"), num(rs[-1], "node_used_io_max_gib")) if v]
        mx, my = (int(v) for v in mesh.split("x"))
        rows.append(dict(nz=nz, nodes=nodes, c=c, mesh=mesh, mx=mx, my=my, s=mean, runs=len(s),
                         noise=noise if (spreads or len(s) > 1) else None, node_s=nodes * mean,
                         mem_max=float(rs[-1]["node_rss_max_gib"]), mem_min=float(rs[-1]["node_rss_min_gib"]),
                         mem_io=max(io) if io else None,
                         eff=efficiency(mx, groups(nz)) * efficiency(my, groups(nz))))
    return sorted(rows, key=lambda r: (r["nz"], r["nodes"], r["node_s"]))


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


# ---------------------------------------------------------------- stage 1

def plan_resolutions(cfg: dict):
    """Resolutions at or above the floor worth benchmarking, each with a node count and mesh.

    Every even RES_Z in ``[RES_FLOOR, RES_FLOOR (1 + RES_ABOVE)]`` is classified by
    :func:`tier`. Each acceptable one gets the most balanced node count between its
    memory minimum (from ``MEM_PER_CELL``, before anything is measured) and 1.35 times
    that. Its predicted cost is :func:`cost_model` times the tier prior over the mesh
    balance. Kept: the floor, the best-predicted resolution of each tier, and nothing
    predicted more than ``RES_BAND`` above the cheapest; at most ``RES_MAX``.

    Parameters
    ----------
    cfg : dict
        Configuration.

    Returns
    -------
    rows : list of dict
        Every resolution examined, with ``keep`` and ``why``.
    lines : list of tuple
        ``(res, nodes, 1, meshx, meshy)`` to benchmark.
    """
    floor = cfg["RES_FLOOR"]
    hi = math.floor(floor * (1 + cfg["RES_ABOVE"]))
    rows = []
    for res in range(floor + floor % 2, hi + 1, 2):
        t = tier(res)
        row = dict(res=res, factors=factor_str(res), tier=t, keep=False, why="")
        if t is None and res != floor:
            row["why"] = "rejected: prime of 11 or more, or 7^2"
            rows.append(row)
            continue
        n0 = min_nodes(prior_variable(res, cfg), cfg)
        best = None
        for nn in range(n0, max(n0, math.ceil(1.35 * n0)) + 1):
            mx, my, e = best_mesh(res, nn)
            key = (round(e, 3), -nn)
            if best is None or key > best[0]:
                best = (key, nn, mx, my, e)
        _, nn, mx, my, e = best
        row.update(nodes=nn, n0=n0, mx=mx, my=my, eff=e,
                   pred=cost_model(res, floor) * TIER_PRIOR.get(t, 1.15) / e)
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
    for r in sorted(ranked[:cfg["RES_MAX"]], key=lambda r: r["res"]):
        mx, my = r["mx"], r["my"]
        both = cfg["ORIENT"] == "both" or (cfg["ORIENT"] == "elongated" and mx > 2 * my)
        lines += [(r["res"], r["nodes"], 1, mx, my)] + ([(r["res"], r["nodes"], 1, my, mx)] if both and mx != my else [])
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
    floor = cfg["RES_FLOOR"]
    out = [f"Stage 1 plan: resolutions from RES_FLOOR = {floor} to {floor * (1 + cfg['RES_ABOVE']):.0f}",
           f"Predicted cost = (N/{floor})^4 log N/log {floor} x tier prior / mesh balance "
           f"(priors, assumed: " + ", ".join(f"{t} {TIER_NAME[t]} {v:.2f}" for t, v in TIER_PRIOR.items()) + ").",
           f"Nodes: memory minimum from MEM_PER_CELL = {cfg['MEM_PER_CELL']:g} GiB per RES_Z^3, grown "
           f"x{growth(cfg):.2f} over {cfg['JOB_HOURS']:g} h, fixed {fixed_per_node(cfg):.0f} GiB/node, "
           f"budget {budget(cfg):.0f} GiB/node; then the most balanced count up to 1.35x that.", "",
           f"{'RES_Z':>6} {'factors':<14}{'tier':>5} {'min':>4} {'nodes':>6} {'mesh':>9} {'bal':>6} "
           f"{'pred':>6}  decision"]
    rejected = [r for r in rows if "pred" not in r]
    for r in rows:
        if "pred" not in r:
            continue
        else:
            out.append(f"{r['res']:>6} {r['factors']:<14}{r['tier'] or '-':>5} {r['n0']:>4} {r['nodes']:>6} "
                       f"{str(r['mx']) + 'x' + str(r['my']):>9} {r['eff']:6.3f} {r['pred']:6.3f}  "
                       f"{'BENCHMARK' if r['keep'] else 'skip'}: {r['why']}")
    if rejected:
        out.append("Rejected (a prime of 11 or more, or 7^2): " + ", ".join(str(r["res"]) for r in rejected))
    out += ["", f"{len(lines)} jobs (orientations: {cfg['ORIENT']}):"]
    out += [f"  RES_Z {a:>4} on {b:>2} nodes, mesh {d}x{e}" for a, b, _, d, e in lines]
    return "\n".join(out) + "\n"


def decide_resolution(rows, cfg):
    """Choose the resolution with the lowest node-hours per unit of simulated time.

    Cost per unit simulated time is ``node-seconds per step / dt`` with dt proportional
    to 1/RES_Z, i.e. ``nodes x s_per_iter x RES_Z`` up to a constant, reported relative
    to the floor. Only resolutions at or above ``RES_FLOOR`` are eligible. If the
    runner-up is within the measured noise, the higher resolution is chosen: it is more
    accurate for the same cost.

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
    floor = cfg["RES_FLOOR"]
    by_res = {}
    for r in rows:
        if r["nz"] >= floor and (r["nz"] not in by_res or r["node_s"] < by_res[r["nz"]]["node_s"]):
            by_res[r["nz"]] = r
    if not by_res:
        raise SystemExit(f"no stage 1 result at or above RES_FLOOR = {floor}")
    ref_res = floor if floor in by_res else min(by_res)
    ref = by_res[ref_res]
    per_res = []
    for nz, r in sorted(by_res.items()):
        var, how = measured_variable(r, cfg)
        per_res.append(dict(r, cost=r["node_s"] * nz / (ref["node_s"] * ref_res),
                            wall=r["s"] * nz / (ref["s"] * ref_res),
                            work=(r["node_s"] / (nz ** 3 * math.log(nz))) / (ref["node_s"] / (ref_res ** 3 * math.log(ref_res))),
                            var=var, how=how, n_min=min_nodes(var, cfg, r["c"]),
                            peak=job_peak(var, r["nodes"], cfg, r["c"])))
    ranked = sorted(per_res, key=lambda r: r["cost"])
    best, why = ranked[0], "lowest node-hours per unit of simulated time"
    if len(ranked) > 1:
        second = ranked[1]
        margin = second["cost"] / best["cost"] - 1
        tol = noise_of(best, second)
        if margin < tol and second["nz"] > best["nz"]:
            best, why = second, (f"within noise ({margin:.1%} < {tol:.1%}) of RES_Z = {ranked[0]['nz']}, "
                                 f"so the higher resolution is taken")

    out = [f"Resolution decision (RES_FLOOR = {floor})", "",
           f"CHOICE: RES_Z = {best['nz']} (best measured layout: {best['nodes']} nodes, mesh {best['mesh']})",
           f"Reason: {why}.", "",
           f"{'RES_Z':>6} {'factors':<12} {'nodes':>5} {'mesh':>8} {'bal':>6} {'s/step':>8} {'node-s':>7} "
           f"{'cost/t':>7} {'wall/t':>7} {'work':>6} {'noise':>6} {'GiB/node':>9} {'min nodes':>9}"]
    for r in per_res:
        noise = f"{r['noise']:.1%}" if r["noise"] is not None else "n/a"
        mem = f"{r['mem_io']:.0f}" if r["mem_io"] else "n/a"
        out.append(f"{r['nz']:>6} {factor_str(r['nz']):<12} {r['nodes']:>5} {r['mesh']:>8} {r['eff']:6.3f} "
                   f"{r['s']:8.4f} {r['node_s']:7.2f} {r['cost']:7.3f} {r['wall']:7.3f} {r['work']:6.3f} "
                   f"{noise:>6} {mem:>9} {r['n_min']:>9}{'  <-' if r['nz'] == best['nz'] else ''}")
    out += ["", "cost/t, wall/t = node-hours and wall time per unit simulated time, relative to the floor,",
            "assuming dt proportional to 1/RES_Z. work = node-seconds per RES_Z^3 ln RES_Z, relative to the floor",
            "(factorisation, balance and communication combined). bal = predicted mesh load balance.",
            "noise = spread of the timed blocks. GiB/node = peak with production outputs.",
            f"min nodes = fewest nodes holding that memory grown x{growth(cfg):.2f} over a {cfg['JOB_HOURS']:g} h job "
            f"within {budget(cfg):.0f} GiB/node.", ""]
    same = [r for r in rows if r["nz"] == best["nz"] and r is not by_res[best["nz"]]]
    if same:
        out.append(f"Other configurations at RES_Z = {best['nz']}:")
        out += [f"  {r['mesh']:>8} on {r['nodes']} nodes: {r['node_s'] / best['node_s']:.2f}x the node-seconds "
                f"per step of {best['mesh']}" for r in same]
        out.append("")
    below = sorted({r["nz"] for r in rows if r["nz"] < floor})
    if below:
        out += [f"Below the floor (not eligible): {', '.join(map(str, below))}", ""]
    hours = projection(best["s"], best["nz"], cfg)
    if hours:
        out.append(f"At {best['nodes']} nodes, a job of SIM_TIME = {cfg['SIM_TIME']:g} would take about "
                   f"{hours:.1f} h ({estimate_dt(cfg)[2]}, scaled as 1/RES_Z).")
    return best, per_res, "\n".join(out) + "\n"


# ---------------------------------------------------------------- stage 2

def plan_nodes(rows, res: int, cfg: dict):
    """Node counts to benchmark for one resolution, from its memory minimum up to ``NODE_SPAN`` x.

    The minimum uses the stage 1 measurement with production outputs. Node counts
    whose best mesh balance is at least min(0.93, best in range - 0.03) are kept, at
    most ``NODE_POINTS``, always including the smallest and largest. Meshes use the
    orientation that won in stage 1.

    Parameters
    ----------
    rows : list of dict
        Stage 1 results.
    res : int
        Chosen RES_Z.
    cfg : dict
        Configuration.

    Returns
    -------
    lines : list of tuple
        Configuration lines.
    text : str
        Explanation (``nodes.txt``).
    """
    at = [r for r in rows if r["nz"] == res]
    if not at:
        raise SystemExit(f"no stage 1 result at RES_Z = {res}")
    ref = min(at, key=lambda r: r["node_s"])
    var, how = measured_variable(ref, cfg)
    n_min = min_nodes(var, cfg)
    n_max = math.ceil(cfg["NODE_SPAN"] * n_min)
    if cfg["NODE_MAX"]:
        n_max = min(n_max, cfg["NODE_MAX"])
    flip = ref["mx"] < ref["my"]
    options = []
    for nn in range(n_min, n_max + 1):
        mx, my, e = best_mesh(res, nn)
        options.append((nn, *((my, mx) if flip else (mx, my)), e))
    top = max(o[3] for o in options)
    cut = min(top - 0.03, 0.93)
    kept = [o for o in options if o[3] >= cut]
    pts = int(cfg["NODE_POINTS"])
    if len(kept) > pts:
        idx = sorted({round(i * (len(kept) - 1) / (pts - 1)) for i in range(pts)})
        kept = [kept[i] for i in idx]
    lines = [(res, nn, 1, mx, my) for nn, mx, my, _ in kept]
    out = [f"Stage 2 plan: node counts for RES_Z = {res}", "",
           f"Memory ({how}): {var:,.0f} GiB depends on the resolution; grown x{growth(cfg):.2f} over "
           f"{cfg['JOB_HOURS']:g} h -> {var * growth(cfg):,.0f} GiB.",
           f"Per node: budget {budget(cfg):.0f} GiB minus {fixed_per_node(cfg):.0f} GiB fixed -> "
           f"minimum {n_min} nodes. Sweep {n_min}-{n_max}.",
           f"Kept: balance >= {cut:.3f}. Orientation: {'smaller' if flip else 'larger'} dimension first "
           f"(the faster one in stage 1, {ref['mesh']}).", "",
           f"{'nodes':>5} {'ranks':>6} {'mesh':>9} {'bal':>6}  kept"]
    for nn, mx, my, e in options:
        out.append(f"{nn:>5} {nn * ranks_per_node(1):>6} {str(mx) + 'x' + str(my):>9} {e:6.3f}  "
                   f"{'yes' if (nn, mx, my, e) in kept else ''}")
    return lines, "\n".join(out) + "\n"


def decide_nodes(rows, res: int, cfg: dict):
    """Choose the node count that minimises node-hours x wall-time^k.

    For a node count N with step time T, node-hours per step are N T and wall time is
    T, so the objective is ``J = N T^(1+k)`` (``k = TIME_WEIGHT``). With ``r = N / N0``,
    speed-up ``S = T0 / T`` and relative efficiency ``E = S / r`` (``N0`` the smallest
    count tested), ``J / J0 = r^-k E^-(1+k)``: a larger count beats the baseline exactly
    when ``E > r^(-k/(1+k))``, i.e. ``E > r^-1/2`` for k = 1 (0.82 at 1.5x the nodes,
    0.71 at 2x). Counts whose projected end-of-job memory exceeds the budget are not
    eligible. Counts whose J is within the measurement noise of the minimum are ties,
    and the fewest nodes among them is chosen.

    Also reported: marginal efficiency against the next smaller count,
    ``(T_prev / T - 1) / (N / N_prev - 1)``, the share of ideal speed-up the added
    nodes delivered.

    Parameters
    ----------
    rows : list of dict
        Stage 2 results.
    res : int
        RES_Z.
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
        if r["nz"] == res and (r["nodes"] not in by_nodes or r["s"] < by_nodes[r["nodes"]]["s"]):
            by_nodes[r["nodes"]] = r
    if len(by_nodes) < 2:
        raise SystemExit(f"need results at two or more node counts for RES_Z = {res}")
    pts = [dict(r) for _, r in sorted(by_nodes.items())]
    base = pts[0]
    prev = None
    for p in pts:
        p["r"] = p["nodes"] / base["nodes"]
        p["S"] = base["s"] / p["s"]
        p["E"] = p["S"] / p["r"]
        p["thr"] = p["r"] ** (-k / (1 + k))
        p["Em"] = None if prev is None else (prev["s"] / p["s"] - 1) / (p["nodes"] / prev["nodes"] - 1)
        p["J"] = p["nodes"] * p["s"] ** (1 + k) / (base["nodes"] * base["s"] ** (1 + k))
        var, _ = measured_variable(p, cfg)
        p["peak"] = job_peak(var, p["nodes"], cfg, p["c"])
        p["fits"] = p["peak"] <= budget(cfg)
        p["ok"] = p["E"] >= p["thr"]
        prev = p
    eligible = [p for p in pts if p["fits"]] or pts
    top = min(eligible, key=lambda p: p["J"])
    ties = [p for p in eligible if p["J"] <= top["J"] * (1 + (1 + k) * noise_of(p, top))]
    best = min(ties, key=lambda p: p["nodes"])
    for p in pts:
        if not p["fits"]:
            p["why"] = "too little memory for a full job"
        elif p is best:
            p["why"] = "CHOSEN" + ("" if p is top else f" (within noise of {top['nodes']} nodes; fewer nodes)")
        elif p in ties:
            p["why"] = "tie (within noise)"
        else:
            p["why"] = f"J {p['J'] / best['J'] - 1:+.0%} vs choice"

    out = [f"Node-count decision for RES_Z = {res}", "",
           f"CHOICE: {best['nodes']} nodes, mesh {best['mesh']}", ""]
    if best is base:
        out.append(f"Reason: no larger node count lowered node-hours x wall-time^{k:g} beyond the noise.")
    else:
        out.append(f"Reason: lowest node-hours x wall-time^{k:g}. Against {base['nodes']} nodes it runs "
                   f"{best['S']:.2f}x as fast ({1 - 1 / best['S']:.0%} less wall time) for "
                   f"{1 / best['E'] - 1:+.0%} node-hours (efficiency {best['E']:.2f}; break-even {best['thr']:.2f}).")
    out += ["", f"Rule: minimise J = nodes x T^{1 + k:g}; against the baseline this means E > r^({-k / (1 + k):.2f}). "
                "Ties within noise go to fewer nodes; counts without memory for a full job are excluded.", "",
            f"{'nodes':>5} {'mesh':>8} {'bal':>6} {'s/step':>8} {'node-s':>7} {'speedup':>8} {'E':>6} "
            f"{'E min':>6} {'marginal':>9} {'J':>6} {'noise':>6} {'GiB end':>8}  result"]
    for p in pts:
        em = f"{p['Em']:.2f}" if p["Em"] is not None else "-"
        noise = f"{p['noise']:.1%}" if p["noise"] is not None else "n/a"
        out.append(f"{p['nodes']:>5} {p['mesh']:>8} {p['eff']:6.3f} {p['s']:8.4f} {p['node_s']:7.2f} {p['S']:8.2f} "
                   f"{p['E']:6.2f} {p['thr']:6.2f} {em:>9} {p['J']:6.3f} {noise:>6} {p['peak']:8.0f}  {p['why']}")
    out += ["", "E = efficiency relative to the smallest count; E min = break-even against it; marginal = share of",
            "ideal speed-up delivered by the nodes added since the previous row; J relative to the smallest count;",
            f"GiB end = projected peak per node after a {cfg['JOB_HOURS']:g} h job (budget {budget(cfg):.0f})."]
    if best["nodes"] * 192 < 1024:
        out += ["", f"Note: {best['nodes']} nodes is {best['nodes'] * 192} cores, below Trillium's stated intent "
                    "of jobs of 1024 cores or more."]
    if best is pts[-1]:
        out += ["", "Note: the largest count tested was chosen; a larger one might be better still."]
    est = estimate_dt(cfg)
    if est and cfg["SIM_TIME"]:
        out += ["", f"Per job of SIM_TIME = {cfg['SIM_TIME']:g} ({est[2]}, scaled as 1/RES_Z):"]
        for p in pts:
            h = projection(p["s"], res, cfg)
            out.append(f"  {p['nodes']:>3} nodes: {h:6.1f} h wall, {h * p['nodes']:7.0f} node-hours")
    out += ["", "For production:", f"  params.conf:          export RES_Z={res} MESHX={best['mx']} MESHY={best['my']}",
            f"  run_sim_trillium.sh:  #SBATCH --nodes={best['nodes']}"]
    return best, pts, "\n".join(out) + "\n"


# ---------------------------------------------------------------- jobs

def write_conf(lines, path: Path, header: str) -> None:
    """Write configuration lines (``RES_Z NODES C MESHX MESHY``).

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
    path.write_text(f"# {header}\n# RES_Z NODES C MESHX MESHY\n" + body)


def tag(line) -> str:
    """Folder name of a configuration.

    Parameters
    ----------
    line : tuple
        ``(res, nodes, c, meshx, meshy)``.

    Returns
    -------
    str
    """
    res, nodes, c, mx, my = line
    return f"{res}_{nodes}n_c{c}_{mx}x{my}"


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
            res, nodes, c, mx, my = line
            rpn = ranks_per_node(c)
            if mx * my != nodes * rpn:
                raise SystemExit(f"{line}: {mx}x{my} = {mx * my} ranks, {nodes} nodes give {nodes * rpn}")
            (outdir / tag(line)).mkdir(exist_ok=True)
            args = ["sbatch", "--parsable", f"--nodes={nodes}", f"--ntasks-per-node={rpn}",
                    f"--cpus-per-task={c}", f"--time={cfg['BENCH_TIME']}", f"--account={cfg['ACCOUNT']}",
                    f"--job-name=bench_{tag(line)}", f"--output={outdir / tag(line)}/slurm_%j.out",
                    str(script), str(res), str(nodes), str(c), str(mx), str(my), str(outdir.resolve())]
            jid = subprocess.run(args, cwd=folder, capture_output=True, text=True, check=True).stdout.strip()
            jid = jid.split(";")[0]
            log.write(f"{tag(line)}\t{jid}\n")
            print(f"  submitted {tag(line):<28} job {jid}")


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
        print(f"  {s['tag']:<28} job {s['job']:<10} {s['state']}{extra}")
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
    print(f"\n{len(lines)} jobs, each about {(cfg['BENCH_SECONDS'] + OVERHEAD_S) / 60:.0f} min; at most about "
          f"{nh:.1f} node-hours ({nh * 192:,.0f} core-hours).")
    if yes:
        return True
    return input("Submit them? [y/N] ").strip().lower() in ("y", "yes")


# ---------------------------------------------------------------- plots

def colours(values) -> dict:
    """ROYGBIV colours by expense: red for the cheapest value, violet for the dearest.

    Parameters
    ----------
    values : iterable
        Sortable values, cheapest first when sorted.

    Returns
    -------
    dict
        Value -> colour.
    """
    vals = sorted(set(values))
    n = len(vals)
    return {v: ROYGBIV[round(i * (len(ROYGBIV) - 1) / max(n - 1, 1))] for i, v in enumerate(vals)}


def plot_resolution(per_res, rows, best, path: Path) -> None:
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
        Output PNG.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    col = colours(r["nz"] for r in per_res)
    fig, ax = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    names = [str(r["nz"]) for r in per_res]
    for a, key, label in ((ax[0, 0], "cost", "node-hours per unit sim time (rel. to floor)"),
                          (ax[0, 1], "wall", "wall time per unit sim time (rel. to floor)"),
                          (ax[1, 0], "work", "node-s per RES_Z³ ln RES_Z (rel. to floor)")):
        vals = [r[key] for r in per_res]
        err = [r[key] * (r["noise"] or 0) for r in per_res]
        bars = a.bar(names, vals, yerr=err, color=[col[r["nz"]] for r in per_res], capsize=4)
        for b, r in zip(bars, per_res):
            b.set_edgecolor("black" if r["nz"] == best["nz"] else "none")
            b.set_linewidth(2)
            a.annotate(f"{r['mesh']}\n{r['nodes']} nodes", (b.get_x() + b.get_width() / 2, b.get_height()),
                       ha="center", va="bottom", fontsize=8, xytext=(0, 3), textcoords="offset points")
        a.axhline(1, color="grey", lw=0.8, ls="--")
        a.set_ylabel(label)
        a.set_xlabel("RES_Z")
    a = ax[1, 1]
    res_list = sorted({r["nz"] for r in rows})
    for i, nz in enumerate(res_list):
        group = [r for r in rows if r["nz"] == nz]
        for j, r in enumerate(group):
            x = i + (j - (len(group) - 1) / 2) * 0.8 / max(len(group), 1)
            a.vlines(x, r["mem_min"], r["mem_max"], color=col.get(nz, "grey"), lw=4)
            top = r["mem_max"]
            if r["mem_io"]:
                a.plot(x, r["mem_io"], "D", color=col.get(nz, "grey"), mec="black")
                top = max(top, r["mem_io"])
            a.annotate(f"{r['mesh']}\n{r['nodes']}n", (x, top), fontsize=7, ha="center",
                       va="bottom", xytext=(0, 4), textcoords="offset points")
    a.set_xticks(range(len(res_list)), [str(v) for v in res_list])
    a.set_xlim(-0.6, len(res_list) - 0.4)
    a.margins(y=0.15)
    a.set_xlabel("RES_Z")
    a.set_ylabel("peak GiB per node: bar = timing phase, ◆ = with outputs")
    fig.suptitle(f"Stage 1: resolution (choice {best['nz']}, outlined)")
    fig.savefig(path, dpi=130)
    plt.close(fig)


def plot_nodes(pts, best, path: Path, cfg: dict) -> None:
    """Stage 2 figure: step time, speed-up, efficiency against the threshold, cost vs time.

    Parameters
    ----------
    pts : list of dict
        From :func:`decide_nodes`.
    best : dict
        Chosen configuration.
    path : pathlib.Path
        Output PNG.
    cfg : dict
        Configuration (``TIME_WEIGHT`` for the threshold curve, memory budget).
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    k = cfg["TIME_WEIGHT"]
    n = np.array([p["nodes"] for p in pts], float)
    s = np.array([p["s"] for p in pts])
    noise = np.array([p["noise"] or 0 for p in pts])
    n0 = n[0]
    fig, ax = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    mk = dict(color="#1f77b4", marker="o")

    a = ax[0, 0]
    a.errorbar(n, s, yerr=s * noise, capsize=3, **mk)
    a.plot(n, s[0] * n0 / n, "--", color="grey", label="ideal")
    a.set_ylabel("s per step")
    a.legend()

    a = ax[0, 1]
    a.plot(n, [p["S"] for p in pts], **mk, label="measured")
    a.plot(n, n / n0, "--", color="grey", label="ideal")
    a.set_ylabel(f"speed-up over {int(n0)} nodes")
    a2 = a.twinx()
    a2.plot(n, [p["peak"] for p in pts], "D:", color="#9467bd", label="GiB/node at job end")
    a2.axhline(budget(cfg), color="#9467bd", lw=0.8, ls="--")
    a2.set_ylabel("projected peak GiB per node", color="#9467bd")
    a.legend(loc="upper left")
    a2.legend(loc="lower right", fontsize=8)

    a = ax[1, 0]
    rr = np.linspace(1, n[-1] / n0, 100)
    a.plot(n, [p["E"] for p in pts], **mk, label="relative efficiency E")
    a.plot(n0 * rr, rr ** (-k / (1 + k)), "--", color="#d62728", label=f"break-even r^(-{k / (1 + k):.2f})")
    em = [(p["nodes"], p["Em"]) for p in pts if p["Em"] is not None]
    if em:
        a.plot(*zip(*em), "s:", color="#2ca02c", label="marginal efficiency")
    a.set_ylim(0, 1.15)
    a.set_ylabel("efficiency")
    a.legend(fontsize=8)

    a = ax[1, 1]
    wall = s / s[0]
    cost = n * s / (n0 * s[0])
    a.plot(cost, wall, "-", color="grey", lw=0.8)
    for p, x, y in zip(pts, cost, wall):
        a.plot(x, y, "o", color="black" if p is best else ("#1f77b4" if p["fits"] else "#d62728"))
        a.annotate(f"{p['nodes']}", (x, y), xytext=(4, 4), textcoords="offset points", fontsize=8)
    a.set_xlabel(f"node-hours (rel. to {int(n0)} nodes)")
    a.set_ylabel(f"wall time (rel. to {int(n0)} nodes)")
    for a in (ax[0, 0], ax[0, 1], ax[1, 0]):
        a.set_xlabel("nodes")
        a.axvline(best["nodes"], color="black", lw=0.8, ls=":")
    fig.suptitle(f"Stage 2: nodes for RES_Z = {pts[0]['nz']} (choice {best['nodes']}, black)")
    fig.savefig(path, dpi=130)
    plt.close(fig)


# ---------------------------------------------------------------- the workflow

def chosen_resolution(res_dir: Path) -> int | None:
    """RES_Z chosen in stage 1, from its decision file.

    Parameters
    ----------
    res_dir : pathlib.Path
        Stage 1 results folder.

    Returns
    -------
    int or None
    """
    f = res_dir / "decision_resolution.txt"
    if not f.exists():
        return None
    m = re.search(r"CHOICE: RES_Z = (\d+)", f.read_text())
    return int(m.group(1)) if m else None


def do_decide_res(res_dir: Path, cfg: dict) -> int:
    """Write the stage 1 decision; return the chosen RES_Z.

    Parameters
    ----------
    res_dir : pathlib.Path
        Stage 1 results folder.
    cfg : dict
        Configuration.

    Returns
    -------
    int
    """
    rows = load(res_dir)
    best, per_res, text = decide_resolution(rows, cfg)
    (res_dir / "decision_resolution.txt").write_text(text)
    plot_resolution(per_res, rows, best, res_dir / "decision_resolution.png")
    print(text + f"(written to {res_dir}/decision_resolution.txt and .png)\n")
    return best["nz"]


def do_decide_nodes(nodes_dir: Path, res: int, cfg: dict) -> None:
    """Write the stage 2 decision.

    Parameters
    ----------
    nodes_dir : pathlib.Path
        Stage 2 results folder.
    res : int
        RES_Z.
    cfg : dict
        Configuration.
    """
    best, pts, text = decide_nodes(load(nodes_dir), res, cfg)
    (nodes_dir / "decision_nodes.txt").write_text(text)
    plot_nodes(pts, best, nodes_dir / "decision_nodes.png", cfg)
    print(text + f"(written to {nodes_dir}/decision_nodes.txt and .png)")


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
        print()
        res = do_decide_res(res_dir, cfg)
        lines, text = plan_nodes(load(res_dir), res, cfg)
        print(text)
        if dry:
            return
        (folder / "nodes.txt").write_text(text)
        write_conf(lines, folder / "nodes.conf", "stage 2 plan (see nodes.txt)")
        if confirm(lines, cfg, yes):
            submit(lines, nodes_dir, cfg, folder)
            print("\nStage 2 submitted. Run 'python -m rbc.bench' again to check on it.")
        return

    print()
    if not print_status("Stage 2 (nodes)", job_status(nodes_dir), nodes_dir):
        return
    print()
    do_decide_nodes(nodes_dir, chosen_resolution(res_dir), cfg)


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
        if any(nodes_dir.rglob("bench.csv")):
            do_decide_nodes(nodes_dir, chosen_resolution(res_dir) or res, cfg)


if __name__ == "__main__":
    main()
