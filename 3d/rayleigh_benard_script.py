"""
Dedalus script simulating 3D horizontally-periodic Rayleigh-Benard convection. No slip boundary conditions.

Usage:
    ded_example.py [options]

Options:
    --Ra=<Ra>                         log_{10} of the Rayleigh number
    --Pr=<Pr>                         log_{10} of the Prandtl number
    --nz=<nz>                         Resolution in chebyshev direction
    --dt=<dt>                         Timestep (initial/max for CFL) [default: 0.00001]
    --cfl_safety=<cfl_safety>         Safety factor of CFL (see dedalus docs) [default: 0.5]
    --cfl_threshold=<cfl_threshold>   Threshold to compute new dt (see ded docs) [default: 0.05]
    --cfl_cadence=<cfl_cadence>       Timesteps before checking dt (see ded docs) [default: 2]
    --sim_time=<time>                 Simulation time
    --basepath=<path>                 Base path for output files
    --stepper=<stepper>               Timestepper [default: RK443]
    --Lx=<Lx>                         Length in x [default: 2]
    --Ly=<Ly>                         Length in y [default: 2]
    --meshx=<size>                    2D process mesh size in x
    --meshy=<size>                    2D process mesh size in y
    --index=<index>                   Restart write index [default: -1]
    --snapshots                       Enable snapshots for visualization
    --cfl                             Enable CFL timestepping
"""

from mpi4py import MPI
import numpy as np
import time
from pathlib import Path
import dedalus.public as de
import logging
logger = logging.getLogger(__name__)

from docopt import docopt
args = docopt(__doc__)

comm = MPI.COMM_WORLD
rank = comm.rank
ncpu = comm.size

# Parameters
basepath = Path(str(args['--basepath']))

dtype = np.float64

Ra = 10**float(args['--Ra'])
Pr = 10**float(args['--Pr'])

Lx, Ly, Lz = int(args['--Lx']), int(args['--Ly']), 1
nz = int(args['--nz'])
nx, ny = int(Lx * nz), int(Ly * nz)

meshx, meshy = int(args['--meshx']), int(args['--meshy'])
dealias = 3/2

stop_sim_time = np.float64(args['--sim_time'])
stepper_name = str(args['--stepper'])
timestepper = {
    'RK222': de.RK222,
    'CNAB2': de.CNAB2,
    'MCNAB2': de.MCNAB2,
    'SBDF4': de.SBDF4,
    'SBDF2': de.SBDF2,
    'RK443': de.RK443,
}.get(stepper_name, de.RK443)

restart_index = int(args['--index'])
snapshots_flag = args['--snapshots']
fh_mode = 'append'

# Iteration parameters
state_time = 60*5
snap_time = 1/10
analysis_time = 1/100
message_num_iters = 500

use_cfl = args['--cfl']
arg_dt = np.float64(args['--dt'])
cfl_safety = np.float64(args['--cfl_safety'])
cfl_threshold = np.float64(args['--cfl_threshold'])
cfl_cadence = int(args['--cfl_cadence'])
max_timestep = 1e-1

# Bases
coords = de.CartesianCoordinates('x', 'y', 'z')
dist = de.Distributor(coords, mesh=[meshx,meshy], dtype=dtype)
xbasis = de.RealFourier(coords['x'], size=nx, bounds=(0, Lx), dealias=dealias)
ybasis = de.RealFourier(coords['y'], size=ny, bounds=(0, Ly), dealias=dealias)
zbasis = de.ChebyshevT(coords['z'],  size=nz, bounds=(0, Lz), dealias=dealias)
x = dist.local_grid(xbasis)
y = dist.local_grid(ybasis)
z = dist.local_grid(zbasis)
ba = (xbasis,ybasis,zbasis)
ba_p = (xbasis,ybasis)

# For analysis only
z_an = dist.Field(name='z_an', bases=(zbasis))
z_an['g'] = z


# Fields
p = dist.Field(name='p', bases=ba)
b = dist.Field(name='b', bases=ba)
u = dist.VectorField(coords, name='u', bases=ba)
tau_p = dist.Field(name='tau_p')
tau_b1 = dist.Field(name='tau_b1', bases=ba_p)
tau_b2 = dist.Field(name='tau_b2', bases=ba_p)
tau_u1 = dist.VectorField(coords, name='tau_u1', bases=ba_p)
tau_u2 = dist.VectorField(coords, name='tau_u2', bases=ba_p)

curl = lambda A: de.Curl(A)
omega = curl(u)

# Substitutions
kappa = (Ra * Pr)**(-1/2)
nu = (Ra / Pr)**(-1/2)

ey, ex, ez = coords.unit_vector_fields(dist)

lift_basis = zbasis.derivative_basis(1)
lift = lambda A: de.Lift(A, lift_basis, -1)

lift_basis2 = zbasis.derivative_basis(2)
lift2 = lambda A: de.Lift(A, lift_basis2, -1)
lift2_2 = lambda A: de.Lift(A, lift_basis2, -2)

b0 = dist.Field(name='b0', bases=zbasis)
b0['g'] = Lz - z


# Problem
problem = de.IVP([p, u, b, tau_p, tau_u1, tau_u2, tau_b1, tau_b2], namespace=locals())
problem.add_equation("div(u) + lift(tau_p) = 0")
# TODO: go to cross(u, curl(u)) form of momentum nonlinearity
problem.add_equation("dt(u) - nu*lap(u) + grad(p) + lift2_2(tau_u1) + lift2(tau_u2) - b*ez = cross(u, omega)")
problem.add_equation("dt(b) + u@grad(b0) - kappa*lap(b) + lift2_2(tau_b1) + lift2(tau_b2) = - (u@grad(b))")
problem.add_equation("b(z=0) = 0")
problem.add_equation("u(z=0) = 0", condition="nx != 0 or ny != 0")
problem.add_equation("p(z=0) = 0", condition="nx == 0 and ny == 0") # Pressure gauge
problem.add_equation("ex@(u)(z=0) = 0", condition="nx == 0 and ny == 0")
problem.add_equation("ey@(u)(z=0) = 0", condition="nx == 0 and ny == 0")
problem.add_equation("ez@tau_u1 = 0", condition="nx == 0 and ny == 0")
problem.add_equation("b(z=Lz) = 0")
problem.add_equation("u(z=Lz) = 0")


# Solver
solver = problem.build_solver(timestepper)
solver.stop_sim_time = stop_sim_time

# Initial conditions OR restart
restart_path = basepath / 'restart' / 'restart.h5'
if not restart_path.exists():

    # Start from conductive profile with noise (parallel-friendly)
    b.fill_random('g', seed=42, distribution='normal', scale=1e-5) # Random noise
    b['g'] *= z * (Lz - z) # Damp noise at walls
    b.low_pass_filter(scales=0.25)
    dt = arg_dt
else:
    write, last_dt = solver.load_state(restart_path, index=restart_index)
    dt = last_dt if use_cfl else arg_dt
    logger.info(f"Loaded restart write={write}, dt={last_dt}")




# Helpers for integrals
kappa_xyz = 1/(Lx*Ly*Lz)
kappa_xy = 1/(Lx*Ly)

T = b0 + b

grad_T = de.grad(T)

basepath.mkdir(parents=True, exist_ok=True)




# State (checkpoint) for restart
(basepath / 'state').mkdir(exist_ok=True)
state = solver.evaluator.add_file_handler(str(basepath / 'state'), wall_dt=state_time, max_writes=10, mode=fh_mode)
state.add_tasks(solver.state)
logger.info(f"State tasks added.")



# Snapshots for visualization
if snapshots_flag:
    (basepath / 'snapshots').mkdir(exist_ok=True)
    snap = solver.evaluator.add_file_handler(basepath / 'snapshots', sim_dt=snap_time, max_writes=60, mode=fh_mode)
    snap.add_task(T, name='temperature')
    snap.add_task(omega@omega, name='vorticity')
    # Velocity components & vorticity
    snap.add_task(u@ez, name='w')

    logger.info(f"Snapshots tasks added.")



# Frequent analysis (domain averages for Nu, Re, energies)
(basepath / 'analysis').mkdir(exist_ok=True)

an = solver.evaluator.add_file_handler(basepath / 'analysis', sim_dt=analysis_time, max_writes=None, mode=fh_mode)


Pr_f = dist.Field(name='Pr')
Pr_f['g'] = nu / kappa

Ra_f = dist.Field(name='Ra')
Ra_f['g'] = (kappa * nu) ** -1

an.add_task(Pr_f, name='Pr')
an.add_task(Ra_f, name='Ra')

an.add_task(z_an, name='z_an')

# Kinetic energy and temps
an.add_task(kappa_xyz*de.integ(u@u), name='avg_K')
an.add_task(kappa_xy*de.integ(de.integ(T, 'x'), 'y'), name='avg_T')
an.add_task(kappa_xy*de.integ(de.integ((u@ex)**2, 'x'), 'y'), name='avg_u_sq')
an.add_task(kappa_xy*de.integ(de.integ((u@ey)**2, 'x'), 'y'), name='avg_v_sq')
an.add_task(kappa_xy*de.integ(de.integ((u@ez)**2, 'x'), 'y'), name='avg_w_sq')
an.add_task(kappa_xy*de.integ(de.integ(T**2, 'x'), 'y'), name='avg_T_sq')

# Nusselt proxies
an.add_task(kappa_xyz*de.integ((u@ez) * T), name='avg_wT')
an.add_task(kappa_xyz*de.integ(grad_T @ grad_T), name='avg_grad_T_sq')
an.add_task(kappa_xyz*de.integ(omega @ omega), name='avg_vorticity_sq')

logger.info(f"Analysis tasks added.")


# CFL
CFL = de.CFL(solver, initial_dt=dt, cadence=1, safety=0.5, threshold=0.1,
             max_change=1.5, min_change=0.5, max_dt=max_timestep)
CFL.add_velocity(u)

cadence = 100
# Flow properties
flow = de.GlobalFlowProperty(solver, cadence=cadence)
flow.add_property(np.sqrt(u@u)/nu, name='Re')

startup_iter = 100
# Main loop
try:
    good_solution = True
    logger.info('Starting loop')
    start_time = time.time()
    main_start = time.time()

    start_iter = solver.iteration

    if use_cfl:
        dts = []
        while solver.proceed and good_solution:
            if solver.iteration == start_iter + startup_iter:
                main_start = time.time()
            dt = CFL.compute_timestep()

            solver.step(dt)
            if (solver.iteration-1) % message_num_iters == 0:
                l_dts = len(dts)
                avg_Re = flow.grid_average('Re')
                good_solution = np.isfinite(avg_Re)

                logger.info('-' * 80)
                if l_dts > 0:
                    logger.info(f"The timestep changed {l_dts} time(s) over the last {message_num_iters} iters.")
                    logger.info(f"Average dt={np.mean(dts):.4e}, min={np.min(dts):.4e}, max={np.max(dts):.4e}")
                else:
                    logger.info(f"The timestep has not changed over the last {message_num_iters} iters.")
                
                logger.info('Iteration={:d}, Time={:.5f}, dt={:.4e}, Re={:.2f}'. format(solver.iteration, solver.sim_time, dt, avg_Re))
                logger.info('-' * 80)

                dts = []
                
                
            new_dt = CFL.compute_timestep()
            if new_dt != dt:
                dts.append(new_dt)
                dt = new_dt
    else:
        while solver.proceed:
            solver.step(dt)
            if (solver.iteration - 1) % message_num_iters == 0:
                logger.info('Iteration: %i, Time: %e, dt: %e' % (
                    solver.iteration, solver.sim_time, dt))
                logger.info('-' * 80)
except:
    logger.error('Exception raised, triggering end of main loop.')
    raise
finally:
    solver.evaluate_handlers_now(dt)
    end_time = time.time()

    startup_time = main_start - start_time
    main_loop_time = end_time - main_start
    DOF = nx*ny*nz
    niter = solver.iteration - startup_iter
    if rank==0:
        print('performance metrics:')
        print('    startup time   : {:}'.format(startup_time))
        print('    main loop time : {:}'.format(main_loop_time))
        print('    main loop iter : {:d}'.format(niter))
        print('    wall time/iter : {:f}'.format(main_loop_time/niter))
        print('          iter/sec : {:f}'.format(niter/main_loop_time))
        print('DOF-cycles/cpu-sec : {:}'.format(DOF*niter/(ncpu*main_loop_time)))
        print('scaling:',
              ' {:d} {:d} {:d}'.format(ncpu, nx, nz),
              ' {:12.7g} {:12.7g} {:12.7g} {:12.7g}'.format(startup_time,
                                                            main_loop_time,
                                                            main_loop_time/niter,
                                                            DOF*niter/(ncpu*main_loop_time)))

    solver.log_stats()
    logger.info("mode-stages/DOF = {}".format(solver.total_modes/DOF))