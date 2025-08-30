"""
3D Rayleigh–Bénard convection using Dedalus 3 (d3).
- Periodic in x & y (Real Fourier), no‑slip in z (Chebyshev T).
- Non‑dimensionalized with box height and free‑fall time:
    kappa = (Ra * Pr)**(-1/2)
    nu    = (Ra / Pr)**(-1/2)
- Uses first‑order (lifted) tau formulation recommended in d3 examples.
- Supports fixed or CFL‑limited timestepping, snapshots, analysis, and restarts.

Usage:
    rayleigh_benard_script3d_dedalus3.py [options]

Options:
    --Ra=<Ra>                Rayleigh number
    --Pr_exp=<Pr>            Exponent on the Prandtl number (base 10)
    --res=<res>              Resolution in modes/unit space
    --dt=<dt>                Timestep (initial/max for CFL)
    --sim_time=<time>        Simulation time
    --basepath=<path>        Base path for output files
    --stepper=<stepper>      Timestepper [default: RK443]
    --Lx=<Lx>                Length in x [default: 2]
    --Ly=<Ly>                Length in y [default: 2]
    --meshx=<size>           2D process mesh size in x [default: 16]
    --meshy=<size>           2D process mesh size in y [default: 12]
    --index=<index>          Restart write index [default: -1]
    --snapshots              Enable snapshots for visualization
    --cfl                    Enable CFL timestepping
    --flow_bc=<bc>           Flow BC in z: no-slip | free-slip [default: no-slip]
"""

from docopt import docopt
import numpy as np
from pathlib import Path
import time
import logging

import dedalus.public as d3

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')

args = docopt(__doc__)

# ---------------------
# Parameters
# ---------------------
Lx, Ly, Lz = (int(args['--Lx']), int(args['--Ly']), 1)
Pr = 10 ** np.float64(args['--Pr_exp'])
Ra = np.float64(args['--Ra'])
res = int(args['--res'])
Nx, Ny, Nz = (Lx * res, Ly * res, Lz * res)
stepper_name = str(args['--stepper'])
flow_bc = str(args['--flow_bc'])
meshx, meshy = int(args['--meshx']), int(args['--meshy'])
basepath = Path(str(args['--basepath']))
arg_dt = np.float64(args['--dt'])
stop_sim_time = np.float64(args['--sim_time'])
snapshots_flag = args['--snapshots']
use_cfl = args['--cfl']
restart_index = int(args['--index'])

dtype = np.float64


# Diffusivities in d3 normalization
kappa = (Ra * Pr) ** (-0.5)
nu = (Ra / Pr) ** (-0.5)

# Coordinates, distributor, bases
coords = d3.CartesianCoordinates('x', 'y', 'z')
dist = d3.Distributor(coords, dtype=dtype, mesh=(meshx, meshy))

dealias = 3/2
xbasis = d3.RealFourier(coords['x'], size=Nx, bounds=(0, Lx), dealias=dealias)
ybasis = d3.RealFourier(coords['y'], size=Ny, bounds=(0, Ly), dealias=dealias)
zbasis = d3.ChebyshevT(coords['z'], size=Nz, bounds=(0, Lz), dealias=dealias)

x, y, z = dist.local_grids(xbasis, ybasis, zbasis)
ex, ey, ez = coords.unit_vector_fields(dist)

# Fields and tau fields (first-order reduction)
p = dist.Field(name='p', bases=(xbasis, ybasis, zbasis))
b = dist.Field(name='b', bases=(xbasis, ybasis, zbasis))  # buoyancy/temperature
u_vec = dist.VectorField(coords, name='u', bases=(xbasis, ybasis, zbasis))

# Tau fields: one set on gradients, one set in PDEs, lifted to derivative basis
lift_basis = zbasis.derivative_basis(1)
lift = lambda A: d3.Lift(A, lift_basis, -1)

tau_p = dist.Field(name='tau_p')
# scalar tau for b
tau_b1 = dist.Field(name='tau_b1', bases=(xbasis, ybasis))
tau_b2 = dist.Field(name='tau_b2', bases=(xbasis, ybasis))
# vector taus for velocity
tau_u1 = dist.VectorField(coords, name='tau_u1', bases=(xbasis, ybasis))
tau_u2 = dist.VectorField(coords, name='tau_u2', bases=(xbasis, ybasis))

# First-order gradients
grad_u = d3.grad(u_vec) + ez * lift(tau_u1)
grad_b = d3.grad(b)     + ez * lift(tau_b1)

# Problem (IVP) in first-order, lifted-tau form
problem = d3.IVP([p, b, u_vec, tau_p, tau_b1, tau_b2, tau_u1, tau_u2], namespace=locals())

# Incompressibility with tau (gauge for vertical boundaries)
problem.add_equation("trace(grad_u) + tau_p = 0")

# Buoyancy/temperature
problem.add_equation("dt(b) - kappa*div(grad_b) + lift(tau_b2) = - u@grad(b)")

# Momentum (Boussinesq)
problem.add_equation("dt(u) - nu*div(grad_u) + grad(p) - b*ez + lift(tau_u2) = - u@grad(u)")

# Boundary conditions: Temperature fixed; flow BC selectable
problem.add_equation("b(z=0) = Lz")
problem.add_equation("b(z=Lz) = 0")

if flow_bc == 'no-slip':
    problem.add_equation("u(z=0) = 0")
    problem.add_equation("u(z=Lz) = 0")
elif flow_bc == 'free-slip':
    # Free-slip: normal velocity zero and tangential stress-free -> du/dz tangential = 0
    # Implement using vector condition helpers
    problem.add_equation("u@ez(z=0) = 0")
    problem.add_equation("u@ez(z=Lz) = 0")
    problem.add_equation("(grad(u)@ex)@ez(z=0) = 0")
    problem.add_equation("(grad(u)@ey)@ez(z=0) = 0")
    problem.add_equation("(grad(u)@ex)@ez(z=Lz) = 0")
    problem.add_equation("(grad(u)@ey)@ez(z=Lz) = 0")
else:
    raise ValueError("--flow_bc must be 'no-slip' or 'free-slip'")

# Pressure gauge (removes null space)
problem.add_equation("integ(p) = 0")

# ---------------------
# Solver & timestepping
# ---------------------
stepper = {
    'RK222': d3.RK222,
    'CNAB2': d3.CNAB2,
    'MCNAB2': d3.MCNAB2,
    'SBDF4': d3.SBDF4,
    'SBDF2': d3.SBDF2,
    'RK443': d3.RK443,
}.get(stepper_name, d3.RK443)
print(stepper)
solver = problem.build_solver(stepper)
logger.info(f"Solver {stepper_name} built.")
solver.stop_sim_time = stop_sim_time

# ---------------------
# Initial conditions OR restart
# ---------------------
restart_path = basepath / 'restart' / 'restart.h5'
if not restart_path.exists():
    # Start from conductive profile with noise (parallel-friendly)
    b.fill_random('g', seed=42, distribution='normal', scale=1e-3)
    b['g'] *= z * (Lz - z)  # damp noise at top/bottom
    b['g'] += Lz - z        # add linear background (conductive)
    fh_mode = 'overwrite'
else:
    write, last_dt = solver.load_state(restart_path, index=restart_index)
    logger.info(f"Loaded restart write={write}, dt={last_dt}")
    fh_mode = 'append'

# ---------------------
# Analysis handlers (snapshots, state, analysis)
# ---------------------
basepath.mkdir(parents=True, exist_ok=True)
"""
# Snapshots for visualization
if snapshots_flag:
    (basepath / 'snapshots').mkdir(exist_ok=True)
    snap = solver.evaluator.add_file_handler(str(basepath / 'snapshots'), iter=200, max_writes=50, mode=fh_mode)
    # Temperature anomaly T = b - (Lz - z)
    snap.add_task(b - (Lz - z), name='temp_anomaly')
    # Velocity components & vorticity
    curl_u = d3.curl(u_vec)
    snap.add_task(u_vec@ex, name='u')
    snap.add_task(u_vec@ey, name='v')
    snap.add_task(u_vec@ez, name='w')
    snap.add_task(curl_u@ex, name='vort_x')
    snap.add_task(curl_u@ey, name='vort_y')
    snap.add_task(curl_u@ez, name='vort_z')

# State (checkpoint) for restart
(basepath / 'state').mkdir(exist_ok=True)
state = solver.evaluator.add_file_handler(str(basepath / 'state'), iter=1000, max_writes=25, mode=fh_mode)
print(type(solver.state))

state.add_tasks(solver.state)

# Frequent analysis (domain averages for Nu, Re, energies)
(basepath / 'analysis').mkdir(exist_ok=True)
an = solver.evaluator.add_file_handler(str(basepath / 'analysis'), iter=50, max_writes=1000, mode=fh_mode)

# Helpers for integrals
kappa_xyz = 1/(Lx*Ly*Lz)
kappa_xy = 1/(Lx*Ly)
T_anom = b - (1/2 - z)

an.add_task(nu/kappa, name='Pr')
an.add_task(1/(kappa*nu), name='Ra')

# Kinetic energy and temps
an.add_task(kappa_xyz*d3.integ(u_vec@u_vec), name='avg_K')
an.add_task(kappa_xy*d3.integ(d3.integ(T_anom, 'x'), 'y'), name='avg_T')
an.add_task(kappa_xy*d3.integ(d3.integ((u_vec@ex)**2, 'x'), 'y'), name='avg_u_sq')
an.add_task(kappa_xy*d3.integ(d3.integ((u_vec@ey)**2, 'x'), 'y'), name='avg_v_sq')
an.add_task(kappa_xy*d3.integ(d3.integ((u_vec@ez)**2, 'x'), 'y'), name='avg_w_sq')
an.add_task(kappa_xy*d3.integ(d3.integ(T_anom**2, 'x'), 'y'), name='avg_T_sq')

# Nusselt proxies
an.add_task(kappa_xyz*d3.integ(u_vec@ez * (T_anom)), name='avg_wT')
# Gradients squared
grad_u_sq = 0
for a in (ex, ey, ez):
    comp = d3.grad(u_vec)@a
    grad_u_sq = grad_u_sq + (comp@comp)

an.add_task(kappa_xyz*d3.integ(grad_u_sq), name='avg_grad_u_sq')
an.add_task(kappa_xyz*d3.integ(d3.grad(T_anom)@d3.grad(T_anom)), name='avg_grad_T_sq')
"""
# ---------------------
# CFL & flow properties
# ---------------------
CFL = None
if use_cfl:
    CFL = d3.CFL(solver, initial_dt=arg_dt, cadence=10, safety=0.5, threshold=0.1,
                 max_change=1.1, min_change=0.5, max_dt=arg_dt)
    CFL.add_velocity(u_vec)

flow = d3.GlobalFlowProperty(solver, cadence=10)
flow.add_property((u_vec@u_vec)/(nu**2), name='Re_sq')

# ---------------------
# Main loop
# ---------------------
try:
    logger.info('Starting loop.')
    logger.info('-' * 80)
    start_time = time.time()
    message_num_iters = 500
    start_iter = solver.iteration

    if use_cfl:
        dt = CFL.compute_timestep()
        dts = []
        while solver.proceed:
            solver.step(dt)
            if (solver.iteration - 1) % message_num_iters == 0:
                l_dts = len(dts)
                if l_dts > 0:
                    logger.info(f"The timestep changed {l_dts} time(s) over the last {message_num_iters} iters.")
                    logger.info(f"Average dt={np.mean(dts):.4e}, min={np.min(dts):.4e}, max={np.max(dts):.4e}")
                logger.info('Iteration: %i, Time: %e, dt: %e, max(Re_sq)=%f' % (
                    solver.iteration, solver.sim_time, dt, flow.max('Re_sq')))
                logger.info('-' * 80)
                dts = []
            new_dt = CFL.compute_timestep()
            if new_dt != dt:
                dts.append(new_dt)
                dt = new_dt
    else:
        dt = arg_dt
        while solver.proceed:
            solver.step(dt)
            if (solver.iteration - 1) % message_num_iters == 0:
                logger.info('Iteration: %i, Time: %e, dt: %e, max(Re_sq)=%f' % (
                    solver.iteration, solver.sim_time, dt, flow.max('Re_sq')))
                logger.info('-' * 80)

except Exception:
    logger.exception('Exception raised, triggering end of main loop.')
    raise
finally:
    solver.log_stats()  # prints nice timing summary
    end_time = time.time()
    total_iter = solver.iteration - start_iter
    total_time = end_time - start_time
    dof = Nx * Ny * Nz * len(problem.variables)  # rough dof count (p, b, u_x,u_y,u_z)
    size = dist.comm_cart.size

    logger.info(f'Total iterations: {total_iter}')
    logger.info('Sim end time: %f' % solver.sim_time)
    logger.info('Run time: %.2f sec' % (total_time))
    logger.info('Run time: %f cpu-hr' % ((total_time) / 3600 * size))
    logger.info('Speed (dof*iter / core-sec): %f' % ((dof * total_iter) / (size * total_time)))
