"""Two-stage benchmark: choose the resolution, then the node count, from measured runs.

Stage 1 (resolution)::

    python -m rbc.bench candidates --target 256 --nodes 6 -o res.conf     # writes res.txt
    python -m rbc.bench submit res.conf --dir bench/res --submit
    python -m rbc.bench decide-res bench/res --target 256                  # decision + plot

Stage 2 (nodes, for the chosen resolution)::

    python -m rbc.bench sweep bench/res --res 240 -o nodes.conf            # writes nodes.txt
    python -m rbc.bench submit nodes.conf --dir bench/nodes --submit
    python -m rbc.bench decide-nodes bench/nodes --res 240                 # decision + plots

Every configuration runs as its own job (``shell/bench_one.sh``) on exactly the nodes it
needs, so no allocated node sits idle. Each job writes ``<dir>/<tag>/bench.csv``; the
decide steps read every such file. Configuration lines are ``RES_Z NODES C MESHX MESHY``,
the same format ``trill_bench`` reads, so a list can also be run in one allocation.

The rules behind each choice are written into the decision files; ``docs/benchmarking.md``
derives them.
"""
from __future__ import annotations

import argparse
import csv
import math
import shlex
import subprocess
from pathlib import Path

from rbc.plan import efficiency, groups, meshes, ranks_per_node

NOISE = 0.03
TIER_PRIOR = {"A": 1.00, "B": 1.03, "C": 1.08}
TIER_NAME = {"A": "3-smooth", "B": "5-smooth", "C": "one factor of 7"}
ROYGBIV = ("#d62728", "#ff7f0e", "#e3b505", "#2ca02c", "#1f77b4", "#3f007d", "#9467bd")


# ---------------------------------------------------------------- arithmetic helpers

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
        (odd, or a prime of 11 or more, or 7^2).
    """
    if res % 2:
        return None
    f = factorise(res)
    if any(p > 7 for p in f) or f.get(7, 0) > 1:
        return None
    if 7 in f:
        return "C"
    return "B" if 5 in f else "A"


def best_mesh(res: int, nodes: int, c: int = 1, cores: int = 192, l3: int = 8):
    """Most balanced mesh for a resolution on a node count.

    Parameters
    ----------
    res : int
        ``RES_Z``.
    nodes : int
        Node count.
    c : int, optional
        cpus-per-task.
    cores, l3 : int, optional
        Node layout (see :func:`rbc.plan.ranks_per_node`).

    Returns
    -------
    tuple
        ``(meshx, meshy, efficiency)`` with ``meshx >= meshy``.
    """
    n = nodes * ranks_per_node(c, cores, l3)
    p0, p1, e = max(meshes(n, groups(res)), key=lambda m: (round(m[2], 4), -abs(m[0] - m[1])))
    return max(p0, p1), min(p0, p1), e


def orientations(mx: int, my: int) -> list[tuple[int, int]]:
    """Mesh orientations worth benchmarking.

    The larger dimension first was faster in all three measured pairs (by 17% for
    128x9 against 9x128, 4% for 40x24, 0.4% for 48x36), so that orientation is always
    run. The other one is added only for meshes elongated by more than 2:1, where the
    difference has been large enough to matter.

    Parameters
    ----------
    mx, my : int
        Mesh with ``mx >= my``.

    Returns
    -------
    list of tuple
    """
    return [(mx, my), (my, mx)] if mx > 2 * my else [(mx, my)]


def cost_model(res: int, target: int) -> float:
    """Work to reach a fixed simulation time, relative to ``target``.

    Work per step scales as the grid points times a logarithm (transforms), RES_Z^3 log
    RES_Z at fixed aspect ratio; the advective CFL limit makes dt scale as 1/RES_Z, so
    the number of steps scales as RES_Z.

    Parameters
    ----------
    res, target : int
        Resolutions.

    Returns
    -------
    float
        ``(res/target)^4 log(res) / log(target)``.
    """
    return (res / target) ** 4 * math.log(res) / math.log(target)


# ---------------------------------------------------------------- reading results

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
    rows = []
    for (nz, nodes, c, mesh), rs in merged.items():
        s = [float(r["s_per_iter"]) for r in rs]
        mean = sum(s) / len(s)
        spreads = [float(r["s_per_iter_spread"]) for r in rs if r.get("s_per_iter_spread")]
        noise = max(spreads + [(max(s) - min(s)) / mean if len(s) > 1 else 0.0])
        mx, my = (int(v) for v in mesh.split("x"))
        rows.append(dict(nz=nz, nodes=nodes, c=c, mesh=mesh, mx=mx, my=my, s=mean, runs=len(s),
                         noise=noise if (spreads or len(s) > 1) else None,
                         node_s=nodes * mean,
                         mem_max=float(rs[-1]["node_rss_max_gib"]),
                         mem_min=float(rs[-1]["node_rss_min_gib"]),
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


# ---------------------------------------------------------------- stage 1: candidates

def candidates(target: int, ref_nodes: int, below: float = 0.07, above: float = 0.125,
               band: float = 0.15, max_res: int = 4, cores: int = 192, l3: int = 8):
    """Resolutions near ``target`` worth benchmarking, with a node count and mesh each.

    Every even RES_Z in ``[target (1 - below), target (1 + above)]`` is classified by
    :func:`tier`. Each acceptable one gets the most balanced node count between
    ``0.9 s`` and ``1.35 s``, with ``s = ref_nodes (res / target)^3`` (the same work and
    memory per node as the target; never much fewer nodes, which might not fit).
    Its predicted cost is :func:`cost_model` times the tier prior over the mesh
    efficiency. Kept: the target, the best-predicted resolution of each tier, and
    nothing predicted more than ``band`` above the cheapest; at most ``max_res``.

    Parameters
    ----------
    target : int
        Resolution you would otherwise use.
    ref_nodes : int
        Node count you would use at ``target``.
    below, above : float, optional
        Window as fractions of ``target``.
    band : float, optional
        Predicted-cost band; outside it the prediction decides without a benchmark.
    max_res : int, optional
        Most resolutions to benchmark.
    cores, l3 : int, optional
        Node layout.

    Returns
    -------
    all_rows : list of dict
        Every resolution examined, with ``keep`` and ``why``.
    lines : list of tuple
        ``(res, nodes, 1, meshx, meshy)`` to benchmark.
    """
    lo, hi = math.ceil(target * (1 - below)), math.floor(target * (1 + above))
    rows = []
    for res in range(lo + lo % 2, hi + 1, 2):
        t = tier(res)
        row = dict(res=res, factors=factor_str(res), tier=t, keep=False, why="")
        if t is None:
            row["why"] = "rejected: prime of 11 or more, or 7^2"
            rows.append(row)
            continue
        scaled = ref_nodes * (res / target) ** 3
        best = None
        for nn in range(max(1, math.ceil(0.9 * scaled)), max(1, math.ceil(1.35 * scaled)) + 1):
            mx, my, e = best_mesh(res, nn, 1, cores, l3)
            key = (round(e, 3), -abs(nn - scaled))
            if best is None or key > best[0]:
                best = (key, nn, mx, my, e)
        _, nn, mx, my, e = best
        row.update(nodes=nn, mx=mx, my=my, eff=e,
                   pred=cost_model(res, target) * TIER_PRIOR[t] / e)
        rows.append(row)

    ok = [r for r in rows if r["tier"]]
    if not ok:
        raise SystemExit("no acceptable resolution in the window; widen --below/--above")
    cheapest = min(r["pred"] for r in ok)
    picks = {r["res"] for r in ok if r["res"] == target}
    for t in "ABC":
        tier_rows = [r for r in ok if r["tier"] == t]
        if tier_rows:
            picks.add(min(tier_rows, key=lambda r: r["pred"])["res"])
    for r in ok:
        if r["res"] not in picks:
            r["why"] = "not the best of its tier"
        elif r["res"] != target and r["pred"] > (1 + band) * cheapest:
            r["why"] = f"predicted more than {band:.0%} above the cheapest: prediction decides"
            picks.discard(r["res"])
    ranked = sorted((r for r in ok if r["res"] in picks), key=lambda r: (r["res"] != target, r["pred"]))
    for i, r in enumerate(ranked):
        if i < max_res:
            r["keep"], r["why"] = True, "target" if r["res"] == target else f"best of tier {r['tier']}"
        else:
            r["why"] = f"over the limit of {max_res} resolutions"
    lines = [(r["res"], r["nodes"], 1, a, b)
             for r in sorted(ranked[:max_res], key=lambda r: r["res"]) for a, b in orientations(r["mx"], r["my"])]
    return rows, lines


def write_conf(lines, path: Path, header: str) -> None:
    """Write configuration lines for :func:`submit` or ``trill_bench``.

    Parameters
    ----------
    lines : list of tuple
        ``(res, nodes, c, meshx, meshy)``.
    path : pathlib.Path
        Output file (ends with a newline so the last line is read).
    header : str
        Comment written at the top.
    """
    body = "".join(f"{' '.join(map(str, line))}\n" for line in lines)
    path.write_text("".join(f"# {h}\n" for h in header.splitlines()) + "# RES_Z NODES C MESHX MESHY\n" + body)


def candidates_report(rows, lines, target, ref_nodes, below, above, band) -> str:
    """Text explaining which resolutions will be benchmarked and why.

    Parameters
    ----------
    rows, lines : list
        From :func:`candidates`.
    target, ref_nodes : int
        Inputs to :func:`candidates`.
    below, above, band : float
        Inputs to :func:`candidates`.

    Returns
    -------
    str
    """
    out = [f"Stage 1 candidates around RES_Z = {target} (reference: {ref_nodes} nodes)",
           f"Window {target * (1 - below):.0f}-{target * (1 + above):.0f}; even sizes only. "
           f"Predicted cost = (N/{target})^4 log N/log {target} x tier prior / mesh efficiency.",
           "Tier priors (assumed, used only to pre-select): " +
           ", ".join(f"{t} {TIER_NAME[t]} {v:.2f}" for t, v in TIER_PRIOR.items()), ""]
    out.append(f"{'RES_Z':>6} {'factors':<14}{'tier':>5} {'nodes':>6} {'mesh':>8} {'eff':>6} "
               f"{'pred':>6}  decision")
    for r in rows:
        if r["tier"] is None:
            out.append(f"{r['res']:>6} {r['factors']:<14}{'-':>5} {'':>6} {'':>8} {'':>6} {'':>6}  {r['why']}")
        else:
            out.append(f"{r['res']:>6} {r['factors']:<14}{r['tier']:>5} {r['nodes']:>6} "
                       f"{r['mx']}x{r['my']:<4} {r['eff']:6.3f} {r['pred']:6.3f}  "
                       f"{'BENCHMARK' if r['keep'] else 'skip'}: {r['why']}")
    out += ["", f"{len(lines)} configurations (larger mesh dimension first; both orientations only "
                f"when elongated more than 2:1).",
            f"Below {target}: only choose these if a lower resolution is acceptable physically "
            f"(use --below 0 otherwise)."]
    return "\n".join(out) + "\n"


# ---------------------------------------------------------------- submitting

def submit(conf: Path, outdir: Path, script: str, time: str, account: str | None,
           cores: int = 192, l3: int = 8, run: bool = False) -> list[str]:
    """One ``sbatch`` per configuration line, each on exactly the nodes it needs.

    Parameters
    ----------
    conf : pathlib.Path
        Configuration lines.
    outdir : pathlib.Path
        Results folder; each run writes ``outdir/<tag>/bench.csv``.
    script : str
        Path to ``bench_one.sh``.
    time : str
        Wall-time request per job.
    account : str or None
        Slurm account (``--account``).
    cores, l3 : int, optional
        Node layout.
    run : bool, optional
        Submit; otherwise only return the commands.

    Returns
    -------
    list of str
        The commands.

    Raises
    ------
    SystemExit
        If a line's mesh does not fill its nodes.
    """
    cmds = []
    for line in conf.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        res, nodes, c, mx, my = (int(v) for v in line.split())
        rpn = ranks_per_node(c, cores, l3)
        if mx * my != nodes * rpn:
            raise SystemExit(f"'{line}': {mx}x{my} = {mx * my} ranks, {nodes} nodes give {nodes * rpn}")
        tag = f"{res}_{nodes}n_c{c}_{mx}x{my}"
        (outdir / tag).mkdir(parents=True, exist_ok=True)
        args = ["sbatch", f"--nodes={nodes}", f"--ntasks-per-node={rpn}", f"--cpus-per-task={c}",
                f"--time={time}", f"--job-name=bench_{tag}", f"--output={outdir / tag}/slurm_%j.out"]
        if account:
            args.append(f"--account={account}")
        args += [script, str(res), str(nodes), str(c), str(mx), str(my), str(outdir.resolve())]
        cmds.append(shlex.join(args))
        if run:
            subprocess.run(args, check=True)
    return cmds


# ---------------------------------------------------------------- stage 1: decision

def decide_resolution(rows, target: int, floor: int | None = None, dt_ref: float | None = None,
                      sim_time: float | None = None):
    """Choose the resolution with the lowest node-hours per unit of simulated time.

    Cost per unit simulated time is ``node_s_per_step / dt`` with dt proportional to
    1/RES_Z, i.e. ``nodes * s_per_iter * RES_Z`` up to a constant; it is reported
    relative to the target's best configuration. Only resolutions at or above
    ``floor`` are eligible. If the runner-up is within the measurement noise, the
    higher resolution is chosen: it is more accurate for the same cost.

    Parameters
    ----------
    rows : list of dict
        From :func:`load`.
    target : int
        Reference resolution.
    floor : int, optional
        Lowest acceptable resolution (default: every benchmarked one).
    dt_ref : float, optional
        Mean timestep at ``target`` (scaled as 1/RES_Z for the others).
    sim_time : float, optional
        Simulation time still to run, for absolute projections.

    Returns
    -------
    best : dict
        Chosen configuration (with derived fields).
    per_res : list of dict
        Best configuration of each resolution.
    text : str
        Contents of ``decision_resolution.txt``.
    """
    by_res = {}
    for r in rows:
        if r["nz"] not in by_res or r["node_s"] < by_res[r["nz"]]["node_s"]:
            by_res[r["nz"]] = r
    if target not in by_res:
        raise SystemExit(f"no result at the target RES_Z = {target}")
    ref = by_res[target]["node_s"] * target
    work_ref = by_res[target]["node_s"] / (target ** 3 * math.log(target))
    per_res = []
    for nz, r in sorted(by_res.items()):
        r = dict(r, cost=r["node_s"] * nz / ref, wall=r["s"] * nz / (by_res[target]["s"] * target),
                 work=r["node_s"] / (nz ** 3 * math.log(nz)) / work_ref, tier=tier(nz))
        per_res.append(r)
    eligible = [r for r in per_res if floor is None or r["nz"] >= floor]
    eligible.sort(key=lambda r: r["cost"])
    best, why = eligible[0], "lowest node-hours per unit of simulated time"
    if len(eligible) > 1:
        second = eligible[1]
        margin = second["cost"] / best["cost"] - 1
        tol = noise_of(best, second)
        if margin < tol and second["nz"] > best["nz"]:
            best, why = second, (f"within noise ({margin:.1%} < {tol:.1%}) of RES_Z = {eligible[0]['nz']}, "
                                 f"so the higher resolution is taken")

    out = [f"Resolution decision (target RES_Z = {target}"
           + (f", floor {floor}" if floor else "") + ")", ""]
    out.append(f"CHOICE: RES_Z = {best['nz']} on {best['nodes']} nodes, mesh {best['mesh']} "
               f"(C = {best['c']})")
    out.append(f"Reason: {why}.")
    out.append("")
    out.append("Per resolution (best configuration of each):")
    out.append(f"{'RES_Z':>6} {'factors':<12} {'nodes':>5} {'mesh':>8} {'eff':>6} {'s/step':>8} "
               f"{'node-s':>7} {'cost/t':>7} {'wall/t':>7} {'work':>6} {'noise':>6}")
    for r in per_res:
        mark = "  <-" if r is best or r["nz"] == best["nz"] else ""
        noise = f"{r['noise']:.1%}" if r["noise"] is not None else "n/a"
        out.append(f"{r['nz']:>6} {factor_str(r['nz']):<12} {r['nodes']:>5} {r['mesh']:>8} {r['eff']:6.3f} "
                   f"{r['s']:8.4f} {r['node_s']:7.2f} {r['cost']:7.3f} {r['wall']:7.3f} {r['work']:6.3f} "
                   f"{noise:>6}{mark}")
    out += ["", "Columns: cost/t = node-hours per unit simulated time and wall/t = wall time per unit",
            "simulated time, both relative to the target and assuming dt proportional to 1/RES_Z;",
            "work = node-seconds per (RES_Z^3 ln RES_Z), relative to the target: how efficiently each",
            "configuration does a unit of work (factorisation, mesh balance and communication combined).",
            "eff = predicted mesh load balance; noise = spread of the timed blocks (n/a: older runs).", ""]
    same = [r for r in rows if r["nz"] == best["nz"] and r is not by_res[best["nz"]]]
    if same:
        out.append(f"Other configurations at RES_Z = {best['nz']}:")
        for r in same:
            out.append(f"  {r['mesh']:>8} on {r['nodes']} nodes, C = {r['c']}: "
                       f"{r['node_s'] / best['node_s']:.2f}x the node-seconds per step of the choice")
        out.append("")
    below = [r for r in per_res if floor is not None and r["nz"] < floor]
    if below:
        out.append("Not eligible (below the floor): " + ", ".join(str(r["nz"]) for r in below))
    if dt_ref and sim_time:
        dt = dt_ref * target / best["nz"]
        steps = sim_time / dt
        wall_h = steps * best["s"] / 3600
        out.append(f"Projection for {sim_time:g} time units (dt ~ {dt:.3g}): {steps:,.0f} steps, "
                   f"{wall_h:,.1f} h wall, {wall_h * best['nodes']:,.0f} node-hours "
                   f"({wall_h * best['nodes'] * 192:,.0f} core-hours), excluding output and restarts.")
    out += ["", "Next: python -m rbc.bench sweep <this folder> --res "
            f"{best['nz']}  (node count for this resolution)"]
    return best, per_res, "\n".join(out) + "\n"


# ---------------------------------------------------------------- stage 2: node sweep

def sweep(rows, res: int, mem_factor: float = 3.5, usable: float = 747.0, headroom: float = 0.85,
          span: float = 2.0, max_points: int = 6, max_nodes: int | None = None,
          cores: int = 192, l3: int = 8):
    """Node counts to benchmark for one resolution, from the memory minimum up to ``span`` x.

    The minimum comes from the stage 1 run: total benchmark memory (fullest node times
    nodes, an upper bound) times ``mem_factor`` (production / benchmark memory; 3.5 was
    measured at RES_Z = 256), divided by ``usable * headroom`` per node. Node counts
    whose best mesh balance is at least min(0.93, best in range - 0.03) are kept, at
    most ``max_points`` of them, always including the smallest and largest.

    Parameters
    ----------
    rows : list of dict
        Stage 1 results (from :func:`load`).
    res : int
        Chosen RES_Z.
    mem_factor : float, optional
        Production memory over benchmark memory.
    usable : float, optional
        GiB a job can use per node (sacct MaxRSS at the out-of-memory kill: 747).
    headroom : float, optional
        Fraction of ``usable`` to plan for.
    span : float, optional
        Largest node count as a multiple of the minimum.
    max_points : int, optional
        Most node counts to benchmark.
    max_nodes : int, optional
        Hard upper limit.
    cores, l3 : int, optional
        Node layout.

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
    bench_total = ref["mem_max"] * ref["nodes"]
    need = bench_total * mem_factor
    n_min = max(1, math.ceil(need / (usable * headroom)))
    n_max = math.ceil(span * n_min)
    if max_nodes:
        n_max = min(n_max, max_nodes)
    options = [(nn, *best_mesh(res, nn, 1, cores, l3)) for nn in range(n_min, n_max + 1)]
    top = max(o[3] for o in options)
    cut = min(top - 0.03, 0.93)
    kept = [o for o in options if o[3] >= cut]
    if len(kept) > max_points:
        idx = sorted({round(i * (len(kept) - 1) / (max_points - 1)) for i in range(max_points)})
        kept = [kept[i] for i in idx]
    lines = [(res, nn, 1, mx, my) for nn, mx, my, _ in kept]
    out = [f"Stage 2 node sweep for RES_Z = {res}", "",
           f"Memory: stage 1 used at most {ref['mem_max']:.0f} GiB/node on {ref['nodes']} nodes "
           f"({ref['mesh']}) -> {bench_total:,.0f} GiB in total.",
           f"Production estimate = {bench_total:,.0f} x {mem_factor} = {need:,.0f} GiB; per node we plan for "
           f"{usable:.0f} x {headroom} = {usable * headroom:.0f} GiB -> minimum {n_min} nodes.",
           f"Sweep {n_min}-{n_max} nodes; kept counts whose best mesh balance is at least {cut:.3f} "
           f"(the smaller of 0.93 and the best in range, {top:.3f}, minus 0.03).", "",
           f"{'nodes':>5} {'ranks':>6} {'mesh':>8} {'eff':>6}  kept"]
    for nn, mx, my, e in options:
        out.append(f"{nn:>5} {nn * ranks_per_node(1, cores, l3):>6} {mx}x{my:<4} {e:6.3f}  "
                   f"{'yes' if (nn, mx, my, e) in kept else ''}")
    return lines, "\n".join(out) + "\n"


def decide_nodes(rows, res: int, k: float = 1.0, dt: float | None = None,
                 sim_time: float | None = None, cores: int = 192):
    """Choose the node count that minimises node-hours x wall-time^k.

    For a node count N with step time T, node-hours per step are N T and wall time is
    T, so the objective is ``J = N T^(1+k)``. With ``r = N / N0``, speed-up
    ``S = T0 / T`` and relative efficiency ``E = S / r`` (``N0`` the smallest count
    tested), ``J / J0 = r^-k E^-(1+k)``, so a larger count beats the baseline exactly
    when ``E > r^(-k/(1+k))``: ``E > r^-1/2`` for ``k = 1``, i.e. 0.82 for 1.5x the
    nodes and 0.71 for 2x. Counts whose J is within the measurement noise of the
    minimum are treated as ties and the fewest nodes among them is chosen.

    Also reported (not used to decide): marginal efficiency against the next smaller
    count, ``(T_prev / T - 1) / (N / N_prev - 1)``, the share of ideal speed-up the
    added nodes delivered; minimising J with k = 1 keeps it above about 0.45.

    Parameters
    ----------
    rows : list of dict
        Stage 2 results (from :func:`load`).
    res : int
        RES_Z.
    k : float, optional
        Weight of wall time against node-hours (0: cheapest; large: fastest).
    dt : float, optional
        Mean timestep, for absolute projections.
    sim_time : float, optional
        Simulation time to run, for absolute projections.
    cores : int, optional
        Cores per node (for the job-size note).

    Returns
    -------
    best : dict
        Chosen configuration.
    pts : list of dict
        One per node count, with ``r``, ``S``, ``E``, ``thr``, ``Em``, ``J``, ``ok``.
    text : str
        Contents of ``decision_nodes.txt``.
    """
    by_nodes = {}
    for r in rows:
        if r["nz"] == res and (r["nodes"] not in by_nodes or r["s"] < by_nodes[r["nodes"]]["s"]):
            by_nodes[r["nodes"]] = r
    if len(by_nodes) < 2:
        raise SystemExit(f"need at least two node counts at RES_Z = {res}")
    pts = [dict(r) for _, r in sorted(by_nodes.items())]
    base = pts[0]
    expo = -k / (1 + k)
    prev = None
    for p in pts:
        p["r"] = p["nodes"] / base["nodes"]
        p["S"] = base["s"] / p["s"]
        p["E"] = p["S"] / p["r"]
        p["thr"] = p["r"] ** expo
        p["Em"] = None if prev is None else (prev["s"] / p["s"] - 1) / (p["nodes"] / prev["nodes"] - 1)
        p["J"] = p["nodes"] * p["s"] ** (1 + k) / (base["nodes"] * base["s"] ** (1 + k))
        prev = p
    jmin = min(p["J"] for p in pts)
    top = min(pts, key=lambda p: p["J"])
    ties = [p for p in pts if p["J"] <= jmin * (1 + (1 + k) * noise_of(p, top))]
    best = min(ties, key=lambda p: p["nodes"])
    for p in pts:
        p["ok"] = p["E"] >= p["thr"]
        if p is best:
            p["why"] = "CHOSEN" + (" (lowest J)" if p is top else f" (J within noise of {top['nodes']} nodes; fewer nodes)")
        elif p in ties:
            p["why"] = "tie (within noise)"
        else:
            p["why"] = f"J {p['J'] / best['J'] - 1:+.0%} vs choice"

    out = [f"Node-count decision for RES_Z = {res}", "",
           f"CHOICE: {best['nodes']} nodes, mesh {best['mesh']} (MESHX={best['mx']} MESHY={best['my']}, "
           f"N={best['mx'] * best['my']})", ""]
    if best is base:
        out.append(f"Reason: no larger node count lowered node-hours x wall-time^{k:g} beyond the noise.")
    else:
        out.append(f"Reason: lowest node-hours x wall-time^{k:g}. Against {base['nodes']} nodes it runs "
                   f"{best['S']:.2f}x as fast ({1 - 1 / best['S']:.0%} less wall time) for "
                   f"{1 / best['E'] - 1:+.0%} node-hours (efficiency {best['E']:.2f}; needed > {best['thr']:.2f}).")
    out += ["", f"Rule: minimise J = nodes x T^{1 + k:g} (node-hours x wall-time^{k:g}); against the baseline this "
                f"means E > r^({expo:.2f}). Ties within noise go to fewer nodes.", "",
            f"{'nodes':>5} {'mesh':>8} {'bal':>6} {'s/step':>8} {'node-s':>7} {'speedup':>8} "
            f"{'E':>6} {'E min':>6} {'marginal':>9} {'J':>6} {'noise':>6}  result"]
    for p in pts:
        em = f"{p['Em']:.2f}" if p["Em"] is not None else "-"
        noise = f"{p['noise']:.1%}" if p["noise"] is not None else "n/a"
        out.append(f"{p['nodes']:>5} {p['mesh']:>8} {p['eff']:6.3f} {p['s']:8.4f} {p['node_s']:7.2f} "
                   f"{p['S']:8.2f} {p['E']:6.2f} {p['thr']:6.2f} {em:>9} {p['J']:6.3f} {noise:>6}  {p['why']}")
    out += ["", "E = efficiency relative to the smallest count; E min = r^(-k/(1+k)), the break-even against it;",
            "marginal = share of ideal speed-up delivered by the nodes added since the previous row;",
            "bal = predicted mesh load balance; J relative to the smallest count."]
    if best["nodes"] * cores < 1024:
        out += ["", f"Note: {best['nodes']} nodes is {best['nodes'] * cores} cores, below Trillium's stated "
                    "intent of jobs of 1024 cores or more."]
    if best is pts[-1]:
        out += ["", "Note: the largest count tested was chosen; a larger one might be better still."]
    if dt and sim_time:
        out += ["", f"Projection for {sim_time:g} time units at dt ~ {dt:.3g} ({sim_time / dt:,.0f} steps), "
                    "excluding output and restarts:"]
        for p in pts:
            wall = sim_time / dt * p["s"] / 3600
            out.append(f"  {p['nodes']:>3} nodes: {wall:8.1f} h wall, {wall * p['nodes']:8.0f} node-hours")
    out += ["", f"params.conf: export MESHX={best['mx']} MESHY={best['my']}; "
                f"run_sim_trillium.sh: #SBATCH --nodes={best['nodes']}"]
    return best, pts, "\n".join(out) + "\n"


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
    for a, key, label in ((ax[0, 0], "cost", "node-hours per unit sim time (rel. to target)"),
                          (ax[0, 1], "wall", "wall time per unit sim time (rel. to target)"),
                          (ax[1, 0], "work", "node-s per RES_Z³ ln RES_Z (rel. to target)")):
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
            a.annotate(f"{r['mesh']}\n{r['nodes']}n", (x, r["mem_max"]), fontsize=7, ha="center",
                       va="bottom", xytext=(0, 2), textcoords="offset points")
    a.set_xticks(range(len(res_list)), [str(v) for v in res_list])
    a.set_xlim(-0.6, len(res_list) - 0.4)
    a.margins(y=0.15)
    a.set_xlabel("RES_Z")
    a.set_ylabel("benchmark peak memory per node, GiB (min-max)")
    fig.suptitle(f"Stage 1: resolution (choice {best['nz']}, outlined)")
    fig.savefig(path, dpi=130)
    plt.close(fig)


def plot_nodes(pts, best, path: Path, k: float = 1.0) -> None:
    """Stage 2 figure: step time, speed-up, efficiency against the threshold, cost vs time.

    Parameters
    ----------
    pts : list of dict
        From :func:`decide_nodes`.
    best : dict
        Chosen configuration.
    path : pathlib.Path
        Output PNG.
    k : float, optional
        The rule's exponent (for the threshold curve).
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    n = np.array([p["nodes"] for p in pts], float)
    s = np.array([p["s"] for p in pts])
    noise = np.array([p["noise"] or 0 for p in pts])
    n0 = n[0]
    fig, ax = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    ok = [p["ok"] for p in pts]
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
    a.legend()

    a = ax[1, 0]
    rr = np.linspace(1, n[-1] / n0, 100)
    a.plot(n, [p["E"] for p in pts], **mk, label="relative efficiency E")
    a.plot(n0 * rr, rr ** (-k / (1 + k)), "--", color="#d62728", label=f"threshold r^(-{k / (1 + k):.2f})")
    em = [(p["nodes"], p["Em"]) for p in pts if p["Em"] is not None]
    if em:
        a.plot(*zip(*em), "s:", color="#2ca02c", label="marginal efficiency")
    for p, good in zip(pts, ok):
        if not good:
            a.plot(p["nodes"], p["E"], "x", color="#d62728", ms=10)
    a.set_ylim(0, 1.15)
    a.set_ylabel("efficiency")
    a.legend(fontsize=8)

    a = ax[1, 1]
    wall = s / s[0]
    cost = n * s / (n0 * s[0])
    a.plot(cost, wall, "-", color="grey", lw=0.8)
    for p, x, y in zip(pts, cost, wall):
        a.plot(x, y, "o", color="black" if p is best else ("#1f77b4" if p["ok"] else "#d62728"))
        a.annotate(f"{p['nodes']}", (x, y), xytext=(4, 4), textcoords="offset points", fontsize=8)
    a.set_xlabel(f"node-hours (rel. to {int(n0)} nodes)")
    a.set_ylabel(f"wall time (rel. to {int(n0)} nodes)")
    for a in ax.flat[:3]:
        a.set_xlabel("nodes")
        a.axvline(best["nodes"], color="black", lw=0.8, ls=":")
    fig.suptitle(f"Stage 2: nodes for RES_Z = {pts[0]['nz']} (choice {best['nodes']})")
    fig.savefig(path, dpi=130)
    plt.close(fig)


# ---------------------------------------------------------------- command line

def get_parser() -> argparse.ArgumentParser:
    """Command-line options (see the module docstring).

    Returns
    -------
    argparse.ArgumentParser
    """
    p = argparse.ArgumentParser(prog="python -m rbc.bench", description=__doc__.split("\n")[0])
    sub = p.add_subparsers(dest="cmd", required=True)

    c = sub.add_parser("candidates", help="stage 1: resolutions to benchmark")
    c.add_argument("--target", type=int, required=True, help="RES_Z you would otherwise use")
    c.add_argument("--nodes", type=int, required=True, help="node count you would use at the target")
    c.add_argument("--below", type=float, default=0.07, help="window below the target (fraction)")
    c.add_argument("--above", type=float, default=0.125, help="window above the target (fraction)")
    c.add_argument("--band", type=float, default=0.15, help="predicted-cost band worth measuring")
    c.add_argument("--max-res", type=int, default=4, help="most resolutions to benchmark")
    c.add_argument("-o", "--out", type=Path, default=Path("res.conf"))

    s = sub.add_parser("submit", help="sbatch one job per configuration line")
    s.add_argument("conf", type=Path)
    s.add_argument("--dir", type=Path, required=True, help="results folder")
    s.add_argument("--script", default="$SCRIPTS_PATH/shell/bench_one.sh")
    s.add_argument("--time", default="00:20:00")
    s.add_argument("--account", default="rrg-goluskin")
    s.add_argument("--submit", action="store_true", help="submit (default: print the commands)")

    d = sub.add_parser("decide-res", help="stage 1 decision")
    d.add_argument("results", type=Path, help="results folder or CSV")
    d.add_argument("--target", type=int, required=True)
    d.add_argument("--floor", type=int, help="lowest acceptable RES_Z")
    d.add_argument("--dt-ref", type=float, help="mean dt at the target, for projections")
    d.add_argument("--sim-time", type=float, help="simulation time to run, for projections")

    w = sub.add_parser("sweep", help="stage 2: node counts to benchmark")
    w.add_argument("results", type=Path, help="stage 1 results folder or CSV")
    w.add_argument("--res", type=int, required=True)
    w.add_argument("--mem-factor", type=float, default=3.5, help="production / benchmark memory")
    w.add_argument("--usable", type=float, default=747.0, help="usable GiB per node")
    w.add_argument("--headroom", type=float, default=0.85)
    w.add_argument("--span", type=float, default=2.0, help="largest count / smallest")
    w.add_argument("--max-points", type=int, default=6)
    w.add_argument("--max-nodes", type=int)
    w.add_argument("-o", "--out", type=Path, default=Path("nodes.conf"))

    n = sub.add_parser("decide-nodes", help="stage 2 decision")
    n.add_argument("results", type=Path, help="stage 2 results folder or CSV")
    n.add_argument("--res", type=int, required=True)
    n.add_argument("--k", type=float, default=1.0, help="weight of wall time (1: node-hours x wall time)")
    n.add_argument("--dt", type=float, help="mean dt, for projections")
    n.add_argument("--sim-time", type=float, help="simulation time to run, for projections")
    return p


def main(argv=None) -> None:
    """Run one subcommand.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments (default: ``sys.argv[1:]``).
    """
    a = get_parser().parse_args(argv)
    if a.cmd == "candidates":
        rows, lines = candidates(a.target, a.nodes, a.below, a.above, a.band, a.max_res)
        text = candidates_report(rows, lines, a.target, a.nodes, a.below, a.above, a.band)
        write_conf(lines, a.out, f"stage 1 around RES_Z = {a.target} (see {a.out.with_suffix('.txt').name})")
        a.out.with_suffix(".txt").write_text(text)
        print(text + f"\nwrote {a.out}")
    elif a.cmd == "submit":
        script = a.script.replace("$SCRIPTS_PATH", __import__("os").environ.get("SCRIPTS_PATH", "$SCRIPTS_PATH"))
        for cmd in submit(a.conf, a.dir, script, a.time, a.account, run=a.submit):
            print(cmd)
        if not a.submit:
            print("(not submitted; add --submit)")
    elif a.cmd == "decide-res":
        rows = load(a.results)
        best, per_res, text = decide_resolution(rows, a.target, a.floor, a.dt_ref, a.sim_time)
        folder = a.results if a.results.is_dir() else a.results.parent
        (folder / "decision_resolution.txt").write_text(text)
        plot_resolution(per_res, rows, best, folder / "decision_resolution.png")
        print(text + f"\nwrote {folder / 'decision_resolution.txt'} and .png")
    elif a.cmd == "sweep":
        lines, text = sweep(load(a.results), a.res, a.mem_factor, a.usable, a.headroom, a.span,
                            a.max_points, a.max_nodes)
        write_conf(lines, a.out, f"stage 2 node sweep for RES_Z = {a.res} (see {a.out.with_suffix('.txt').name})")
        a.out.with_suffix(".txt").write_text(text)
        print(text + f"\nwrote {a.out}")
    elif a.cmd == "decide-nodes":
        best, pts, text = decide_nodes(load(a.results), a.res, a.k, a.dt, a.sim_time)
        folder = a.results if a.results.is_dir() else a.results.parent
        (folder / "decision_nodes.txt").write_text(text)
        plot_nodes(pts, best, folder / "decision_nodes.png", a.k)
        print(text + f"\nwrote {folder / 'decision_nodes.txt'} and .png")


if __name__ == "__main__":
    main()
