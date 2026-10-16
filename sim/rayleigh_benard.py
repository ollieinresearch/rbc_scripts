"""3D horizontally periodic Rayleigh–Bénard convection in Dedalus v3.

No-slip, fixed-temperature plates at z = 0 (T = 1) and z = 1 (T = 0); periodic in
x and y with Lx = Ly = 2 by default. Free-fall units: lengths in H, temperature in
the plate difference, velocity in sqrt(g alpha Delta T H), so that

    kappa = (Ra Pr)^(-1/2),   nu = (Ra / Pr)^(-1/2).

The temperature is split as T = b0 + b with the conductive profile b0 = 1 - z, and
b is the prognostic variable. Launched by ``shell/workflow.sh`` (``srun_sim`` or
``trill_run_sim``); every option below is set from ``params.conf``.

Outputs, all under ``--basepath``:

``state/``
    checkpoints, one write per set, every ``--state_freq`` wall seconds and at the end.
``analysis/``
    volume averages used for Nu and Re (``rbc.analysis.nusselt``).
``horizontal_analysis/``
    horizontally averaged profiles (``rbc.analysis.profiles``).
``spectra/``
    coefficient spectra computed in the run (``rbc.spectra.core.InSimSpectra``).
``snapshots/``
    full fields, only with ``--snapshots``.
``outputs/jobs.csv``
    one row per job: resources used and the stretch of simulation time covered.
"""
from __future__ import annotations

import argparse
import csv
import faulthandler
import logging
import os
import socket
import time
from datetime import datetime
from pathlib import Path

logger = logging.getLogger(__name__)

JOB_FIELDS = ("finished", "job_id", "cluster", "nodes", "ranks", "cpus_per_rank", "mesh",
              "nx", "ny", "nz", "t_start", "t_end", "iter_start", "iter_end",
              "wall_s", "setup_s", "loop_s", "core_hours", "stop_reason")


def get_parser() -> argparse.ArgumentParser:
    """Command-line options of the simulation script.

    Returns
    -------
    argparse.ArgumentParser
        Parser for every option; ``shell/workflow.sh`` (``rb_args``) builds the
        command line from ``params.conf``.
    """
    p = argparse.ArgumentParser(prog="rayleigh_benard.py", allow_abbrev=False,
                                description="3D Rayleigh-Benard convection (Dedalus v3, free-fall units).")
    g = p.add_argument_group("physics and domain")
    g.add_argument("--Ra", type=float, required=True, help="log10 of the Rayleigh number")
    g.add_argument("--Pr", type=float, required=True, help="log10 of the Prandtl number")
    g.add_argument("--Lx", type=float, default=2.0, help="box length in x (H = 1)")
    g.add_argument("--Ly", type=float, default=2.0, help="box length in y")
    g = p.add_argument_group("resolution and parallel layout")
    g.add_argument("--nz", type=int, required=True, help="Chebyshev modes in z")
    g.add_argument("--nx", type=int, help="Fourier modes in x and y (default: gamma * nz)")
    g.add_argument("--gamma", type=float, default=2.0, help="nx / nz when --nx is not given")
    g.add_argument("--meshx", type=int, required=True, help="process mesh size in x")
    g.add_argument("--meshy", type=int, required=True, help="process mesh size in y")
    g = p.add_argument_group("time stepping")
    g.add_argument("--sim_time", type=float, required=True, help="stop at this simulation time")
    g.add_argument("--stepper", default="RK222", help="Dedalus timestepper class name")
    g.add_argument("--dt", type=float, default=1e-5, help="fixed timestep, or the first one with --cfl")
    g.add_argument("--cfl", action="store_true", help="adapt the timestep with a CFL condition")
    g.add_argument("--cfl_safety", type=float, default=0.5, help="CFL safety factor")
    g.add_argument("--cfl_threshold", type=float, default=0.05,
                   help="relative change below which the timestep is left alone")
    g.add_argument("--cfl_cadence", type=int, default=2, help="iterations between CFL updates")
    g.add_argument("--index", type=int, default=-1, help="write index within restart/restart.h5")
    g = p.add_argument_group("output")
    g.add_argument("--basepath", type=Path, required=True, help="run directory")
    g.add_argument("--state_freq", type=float, default=3600, help="wall seconds between checkpoints")
    g.add_argument("--state_layout", choices=("c", "g"), default="c",
                   help="checkpoint in coefficient space (exact, 3.4x smaller) or grid space")
    g.add_argument("--a_freq", type=float, default=0.005, help="simulation time between analysis writes")
    g.add_argument("--snaps_freq", type=float, default=1 / 60, help="simulation time between snapshots")
    g.add_argument("--spec_freq", type=float, default=0.05, help="simulation time between spectra")
    g.add_argument("--snapshots", action="store_true", help="write full-field snapshots")
    g.add_argument("--no_spectra", action="store_true", help="do not compute spectra during the run")
    g.add_argument("--par", default="virtual", choices=("virtual", "mpio", "gather"),
                   help="Dedalus file-handler mode for state and snapshots")
    g.add_argument("--par_small", default="gather", choices=("virtual", "mpio", "gather"),
                   help="file-handler mode for analysis and horizontal_analysis")
    g = p.add_argument_group("filesystem workarounds and checks")
    g.add_argument("--tmp", default="False",
                   help="Dedalus FILEHANDLER_TOUCH_TMPFILE (irrelevant with the rbc.vast patch)")
    g.add_argument("--no_vast_patch", action="store_true",
                   help="use Dedalus' own file-handler code (e.g. with fixed_get_file.py)")
    g.add_argument("--verify_load", action="store_true",
                   help="compare the fast and Dedalus restart loaders, then exit")
    return p


def log_job(path: Path, row: dict) -> None:
    """Append one job's record to ``outputs/jobs.csv`` (call on one rank only).

    Parameters
    ----------
    path : pathlib.Path
        The CSV file; created with a header if absent.
    row : dict
        Values for :data:`JOB_FIELDS`; missing keys are left blank.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    new = not path.exists()
    with open(path, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=JOB_FIELDS, extrasaction="ignore")
        if new:
            w.writeheader()
        w.writerow(row)


def main(argv=None) -> None:
    """Build the problem, load the restart if present, and run the time loop.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments (default: ``sys.argv[1:]``).
    """
    t_script = time.time()
    args = get_parser().parse_args(argv)

    import numpy as np
    from mpi4py import MPI

    import dedalus.public as de
    from dedalus.core import evaluator

    import rbc.vast
    from rbc.shutdown import GracefulStop
    from rbc.spectra.core import InSimSpectra

    comm = MPI.COMM_WORLD
    if not args.no_vast_patch:
        rbc.vast.patch()
    # Read by Dedalus at import, so setting config[...] here would have no effect.
    evaluator.FILEHANDLER_TOUCH_TMPFILE = args.tmp.lower() in ("true", "1", "yes")

    basepath = args.basepath
    Ra, Pr = 10 ** args.Ra, 10 ** args.Pr
    Lx, Ly, Lz = args.Lx, args.Ly, 1.0
    nz = args.nz
    nx = ny = args.nx if args.nx else int(args.gamma * nz)
    mesh = [args.meshx, args.meshy]
    stepper = getattr(de, args.stepper)
    use_cfl = args.cfl
    # Not "dt": the equations are evaluated in locals(), where dt(...) is the time derivative.
    step_dt = args.dt
    message_iters = 500

    # Bases and fields
    dealias = 3 / 2
    coords = de.CartesianCoordinates("x", "y", "z")
    dist = de.Distributor(coords, mesh=mesh, dtype=np.float64)
    xbasis = de.RealFourier(coords["x"], size=nx, bounds=(0, Lx), dealias=dealias)
    ybasis = de.RealFourier(coords["y"], size=ny, bounds=(0, Ly), dealias=dealias)
    zbasis = de.ChebyshevT(coords["z"], size=nz, bounds=(0, Lz), dealias=dealias)
    x, y, z = dist.local_grids(xbasis, ybasis, zbasis)
    ba, ba_p = (xbasis, ybasis, zbasis), (xbasis, ybasis)

    p = dist.Field(name="p", bases=ba)
    b = dist.Field(name="b", bases=ba)
    u = dist.VectorField(coords, name="u", bases=ba)
    tau_p = dist.Field(name="tau_p")
    tau_b1 = dist.Field(name="tau_b1", bases=ba_p)
    tau_b2 = dist.Field(name="tau_b2", bases=ba_p)
    tau_u1 = dist.VectorField(coords, name="tau_u1", bases=ba_p)
    tau_u2 = dist.VectorField(coords, name="tau_u2", bases=ba_p)

    kappa = (Ra * Pr) ** (-1 / 2)
    nu = (Ra / Pr) ** (-1 / 2)
    ex, ey, ez = coords.unit_vector_fields(dist)
    omega = de.Curl(u)
    lift_basis2 = zbasis.derivative_basis(2)
    lift = lambda A: de.Lift(A, zbasis.derivative_basis(1), -1)
    lift2 = lambda A: de.Lift(A, lift_basis2, -1)
    lift2_2 = lambda A: de.Lift(A, lift_basis2, -2)
    b0 = dist.Field(name="b0", bases=zbasis)
    b0["g"] = Lz - z

    problem = de.IVP([p, u, b, tau_p, tau_u1, tau_u2, tau_b1, tau_b2], namespace=locals())
    problem.add_equation("div(u) + lift(tau_p) = 0")
    problem.add_equation("dt(u) - nu*lap(u) + grad(p) + lift2_2(tau_u1) + lift2(tau_u2) - b*ez = cross(u, omega)")
    problem.add_equation("dt(b) + u@grad(b0) - kappa*lap(b) + lift2_2(tau_b1) + lift2(tau_b2) = - (u@grad(b))")
    problem.add_equation("b(z=0) = 0")
    problem.add_equation("u(z=0) = 0", condition="nx != 0 or ny != 0")
    problem.add_equation("p(z=0) = 0", condition="nx == 0 and ny == 0")
    problem.add_equation("ex@(u)(z=0) = 0", condition="nx == 0 and ny == 0")
    problem.add_equation("ey@(u)(z=0) = 0", condition="nx == 0 and ny == 0")
    problem.add_equation("ez@tau_u1 = 0", condition="nx == 0 and ny == 0")
    problem.add_equation("b(z=Lz) = 0")
    problem.add_equation("u(z=Lz) = 0")

    solver = problem.build_solver(stepper)
    solver.stop_sim_time = args.sim_time

    # Watchdog: if no progress for RBC_WATCHDOG seconds, dump the Python stack of every
    # node's first rank (all ranks with RBC_WATCHDOG_ALL=1) to watchdog/ and exit, which
    # makes mpirun end the job. Uses a C thread, so it works while blocked in MPI or HDF5.
    wd_time = float(os.environ.get("RBC_WATCHDOG", "0") or 0)
    wd_file, wd_exit = None, False
    if wd_time > 0:
        node = comm.Split_type(MPI.COMM_TYPE_SHARED)
        wd_exit = node.rank == 0          # leaders dump and exit 30 s after the others have dumped
        if node.rank == 0 or os.environ.get("RBC_WATCHDOG_ALL") == "1":
            (basepath / "watchdog").mkdir(parents=True, exist_ok=True)
            wd_file = open(basepath / "watchdog" / f"rank{comm.rank:05d}_{socket.gethostname()}.txt", "w")
            wd_file.write(f"rank {comm.rank} host {socket.gethostname()} pid {os.getpid()}: "
                          f"stack dumped if no progress for {wd_time:.0f} s\n")
            wd_file.flush()
        node.Free()

    def watchdog(stage=""):
        """(Re)arm the stack dump; ``stage`` is logged to the watchdog file."""
        if wd_file is not None:
            if stage:
                wd_file.write(f"{time.strftime('%H:%M:%S')} armed: {stage}\n")
                wd_file.flush()
            faulthandler.dump_traceback_later(wd_time + 30 * wd_exit, exit=wd_exit, file=wd_file)

    watchdog("restart load / setup")

    restart_path = basepath / "restart" / "restart.h5"
    if restart_path.exists():
        if args.verify_load:
            rbc.vast.verify_load(solver, restart_path, args.index, logger)
            raise SystemExit(0)
        write, last_dt = rbc.vast.load_state(solver, restart_path, args.index, logger)
        if use_cfl:
            step_dt = max(1e-8, float(last_dt))
    else:
        b.fill_random("g", seed=42, distribution="normal", scale=1e-5)
        b["g"] *= z * (Lz - z)
        b.low_pass_filter(scales=0.25)

    # Output
    T = b0 + b
    grad_T = de.grad(T)
    vol, area = 1 / (Lx * Ly * Lz), 1 / (Lx * Ly)
    havg = lambda f: area * de.integ(de.integ(f, "x"), "y")
    for d in ("state", "analysis", "horizontal_analysis") + (("snapshots",) if args.snapshots else ()):
        (basepath / d).mkdir(parents=True, exist_ok=True)

    state = solver.evaluator.add_file_handler(basepath / "state", wall_dt=args.state_freq,
                                              max_writes=1, mode="append", parallel=args.par)
    state.add_tasks(solver.state, layout=args.state_layout)

    if args.snapshots:
        snap = solver.evaluator.add_file_handler(basepath / "snapshots", sim_dt=args.snaps_freq,
                                                 max_writes=30, mode="append", parallel=args.par)
        snap.add_task(T, name="temperature")
        snap.add_task(u, name="u")
        snap.add_task(omega @ omega, name="vorticity")   # |omega|^2

    Ra_f, Pr_f = dist.Field(name="Ra"), dist.Field(name="Pr")
    Ra_f["g"], Pr_f["g"] = Ra, Pr

    an = solver.evaluator.add_file_handler(basepath / "analysis", sim_dt=args.a_freq, mode="append",
                                           parallel=args.par_small)
    an.add_task(Ra_f, name="Ra")
    an.add_task(Pr_f, name="Pr")
    an.add_task(vol * de.integ((u @ ez) * T), name="avg_wT")
    an.add_task(vol * de.integ(grad_T @ grad_T), name="avg_grad_T_sq")
    an.add_task(vol * de.integ(omega @ omega), name="avg_vorticity_sq")
    an.add_task(vol * de.integ(u @ u), name="avg_K")

    an_h = solver.evaluator.add_file_handler(basepath / "horizontal_analysis", sim_dt=args.a_freq,
                                             mode="append", parallel=args.par_small)
    an_h.add_task(havg(T), name="avg_T")
    an_h.add_task(havg((u @ ex) ** 2), name="avg_u_sq")
    an_h.add_task(havg((u @ ey) ** 2), name="avg_v_sq")
    an_h.add_task(havg((u @ ez) ** 2), name="avg_w_sq")
    an_h.add_task(havg(T ** 2), name="avg_T_sq")
    an_h.add_task(havg((u @ ez) * T), name="havg_wT")

    spectra = None
    if not args.no_spectra:
        spectra = InSimSpectra(dist, u, b, args.spec_freq, basepath / "spectra",
                               attrs=dict(Ra=Ra, Pr=Pr, Lx=Lx, Ly=Ly, Lz=Lz), comm=comm)

    flow = de.GlobalFlowProperty(solver, cadence=message_iters)
    flow.add_property(np.sqrt(u @ u) / nu, name="Re")
    CFL = None
    if use_cfl:
        CFL = de.CFL(solver, initial_dt=step_dt, cadence=args.cfl_cadence,
                     safety=args.cfl_safety, threshold=args.cfl_threshold, max_dt=0.1)
        CFL.add_velocity(u)
        step_dt = CFL.compute_timestep()

    stop = GracefulStop(comm, logger=logger)
    good = True
    start_iter, start_time, warmup = solver.iteration, solver.sim_time, 100
    t_start = t_main = time.time()
    dts = []
    logger.info(f"Starting loop: iteration={solver.iteration}, t={solver.sim_time:.5f}, dt={step_dt:.4e}")
    watchdog("first step (all file handlers fire)")
    try:
        while solver.proceed:
            if solver.iteration == start_iter + warmup:
                t_main = time.time()
            solver.step(step_dt)
            if solver.iteration == start_iter + 1:
                logger.info(f"First step done ({time.time() - t_start:.1f} s, includes first handler writes)")
            if solver.iteration % 10 == 0:
                watchdog()
            if spectra is not None:
                spectra.maybe_write(solver)
            if stop.check(solver.iteration):
                break
            if (solver.iteration - 1) % message_iters == 0:
                Re = flow.grid_average("Re")
                if not np.isfinite(Re):
                    good = False
                    logger.error("Re is not finite; stopping.")
                    break
                if dts:
                    logger.info(f"dt changed {len(dts)}x: mean {np.mean(dts):.3e}, "
                                f"min {np.min(dts):.3e}, max {np.max(dts):.3e}")
                logger.info(f"Iteration={solver.iteration}, t={solver.sim_time:.5f}, dt={step_dt:.4e}, Re={Re:.2f}")
                dts = []
            if CFL is not None:
                new_dt = CFL.compute_timestep()
                if new_dt != step_dt:
                    dts.append(new_dt)
                    step_dt = new_dt
    except Exception:
        logger.exception("Exception in main loop")
        comm.Abort(1)    # other ranks would otherwise wait in a collective until the time limit
    finally:
        watchdog("final checkpoint")
        if good and stop.safe_to_write:
            solver.evaluate_handlers(dt=step_dt)      # final checkpoint: next job resumes exactly here
        elif not good:
            logger.error("Solution not finite: no final checkpoint written.")
        stop.record(basepath)
        niter = max(solver.iteration - start_iter - warmup, 1)
        loop = time.time() - t_main
        if comm.rank == 0:
            dof = nx * ny * nz
            print(f"startup {t_main - t_start:.1f} s | loop {loop:.1f} s | {niter} iters | "
                  f"{loop / niter:.4f} s/iter | {dof * niter / (comm.size * loop):.4g} DOF-cycles/cpu-s")
            wall = time.time() - t_script
            cpr = int(os.environ.get("SLURM_CPUS_PER_TASK") or os.environ.get("OMP_NUM_THREADS") or 1)
            log_job(basepath / "outputs" / "jobs.csv", dict(
                finished=datetime.now().isoformat(timespec="seconds"),
                job_id=os.environ.get("SLURM_JOB_ID", ""), cluster=os.environ.get("CC_CLUSTER", ""),
                nodes=os.environ.get("SLURM_JOB_NUM_NODES", ""), ranks=comm.size, cpus_per_rank=cpr,
                mesh=f"{mesh[0]}x{mesh[1]}", nx=nx, ny=ny, nz=nz,
                t_start=f"{start_time:.6g}", t_end=f"{solver.sim_time:.6g}",
                iter_start=start_iter, iter_end=solver.iteration, wall_s=f"{wall:.1f}",
                setup_s=f"{t_start - t_script:.1f}", loop_s=f"{time.time() - t_start:.1f}",
                core_hours=f"{wall * comm.size * cpr / 3600:.2f}", stop_reason=stop.reason or "finished"))
        solver.log_stats()
        if wd_file is not None:
            faulthandler.cancel_dump_traceback_later()
            wd_file.close()
            Path(wd_file.name).unlink(missing_ok=True)   # kept only if the job died stuck


if __name__ == "__main__":
    main()
