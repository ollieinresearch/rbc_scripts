"""
2D Rayleigh Benard convection using Dedalus. The equations of motion are 
non-dimensionalized using the buoyancy time. Boundary conditions are no slip in 
the z-direction and periodic in the x direction. Using the '--snapshots' 
option, temperature and vorticity can be saved at regular intervals for 
visualization. Various (many) other quantities are saved for analysis. This 
script supports restarts, and can integrate using a variable timestep (--cfl).

Usage:
    rayleigh_benard_script2.py [options]

Options:
    --Ra=<Ra>                Rayleigh number
    --Pr_exp=<Pr>            Exponent on the Prandtl number (base 10)
    --res=<res>              Resolution in modes/unit space
    --dt=<dt>                Timestep [default: 0.02]
    --sim_time=<time>        Simulation time
    --basepath=<path>        basepath for output files
    --stepper=<stepper>      The solver to be used for timestepping [default: RK222]
    --Gamma=<Gamma>          The aspect ratio of the cell [default: 2]
    --index=<index>          The index of the restart file to begin simulations from [default: -1]
    --snapshots              Flag to activate the snapshots filehandler that saves the temperature for visualisations.
    --cfl                    Flag to activate CFL variable time-stepping.
"""

from docopt import docopt
import numpy as np
from mpi4py import MPI # pyright: ignore
import time
from pathlib import Path

from dedalus import public as de
from dedalus.extras import flow_tools

import logging
logger = logging.getLogger(__name__)

from dedalus.tools.config import config
config['logging']['stdout_level'] = 'debug'

args = docopt(__doc__)

# Fix non-dimensional height of the layer at 1
Lx, Lz = (np.float64(args['--Gamma']), 1) 
# Rayleigh and Prandtl number
Pr = 10**np.float64(args['--Pr_exp'])
Ra = np.float64(args['--Ra'])
# Resolution for the simulation
nx, nz = (int(Lx*int(args['--res'])), int(Lz*int(args['--res'])))
# Problem and solver params
stepper = str(args['--stepper'])
flow_bc = "no-slip"
arg_dt = np.float64(args['--dt'])
stop_sim_time = np.float64(args['--sim_time'])
cfl = args['--cfl']
# File parameters
basepath = Path(str(args['--basepath']))
snapshots = args['--snapshots']
restart_index = int(args['--index'])

# Iterations between saves
state_iters = 5000
analysis_iters = 25
field_analysis_iters = 200
snapshots_iters = 200
message_num_iters = 500

if Ra/Pr >= 1e10:
    # The timestep goes very low when Ra big and Pr low; saves should be 
    # less frequent.
    message_num_iters *= 10
    state_iters *= 10
    snapshots_iters *= 2

# Create bases and domain
x_basis = de.Fourier('x', nx, interval=(0, Lx), dealias=3/2)
z_basis = de.Chebyshev('z', nz, interval=(-Lz/2, Lz/2), dealias=3/2)
domain = de.Domain([x_basis, z_basis], grid_dtype=np.float64)

# Problem setup: 2D Boussinesq equations; non-dimensionalized by the buoyancy 
# time and written in terms of the vorticity.
problem = de.IVP(domain, variables=['p','T','u','w','Tz','oy'])

# Boundary conditions
problem.meta['p','T','u','w']['z']['dirichlet'] = True

# Parameters
problem.parameters['P'] = (Ra * Pr)**(-1/2)
problem.parameters['R'] = (Ra / Pr)**(-1/2)
problem.parameters['F'] = F = 1
problem.parameters['kappa_xz'] = 1/(Lx*Lz)
problem.parameters['kappa_x'] = 1/Lx


# Equations describing the PDE, and first order reductions
problem.add_equation("dx(u) + dz(w) = 0")
problem.add_equation("dt(T) - P*(d(T,x=2) + dz(Tz)) - F*w  = -(u*dx(T) + w*Tz)")
problem.add_equation("dt(u) - R*dz(oy) + dx(p)             = -oy*w")
problem.add_equation("dt(w) + R*dx(oy) + dz(p) - T         = oy*u")
problem.add_equation("Tz - dz(T) = 0")
problem.add_equation("oy + dx(w) - dz(u) = 0")


problem.add_bc("left(T) = 0")
problem.add_bc("right(T) = 0")

if flow_bc == "no-slip":
    problem.add_bc("left(u) = 0")
    problem.add_bc("right(u) = 0")

if flow_bc == "stress-free":
    problem.add_bc("left(oy) = 0")
    problem.add_bc("right(oy) = 0")

problem.add_bc("left(w) = 0")
problem.add_bc("right(w) = 0", condition="(nx != 0)")
problem.add_bc("right(p) = 0", condition="(nx == 0)")


# Build solver; timestepper is passed as a CLA
solver=None
if stepper == 'RK222':
    solver = problem.build_solver(de.timesteppers.RK222)
elif stepper == 'CNAB2':
    solver = problem.build_solver(de.timesteppers.CNAB2)
elif stepper == 'MCNAB2':
    solver = problem.build_solver(de.timesteppers.MCNAB2)
elif stepper == 'SBDF4':
    solver = problem.build_solver(de.timesteppers.SBDF4)
else:
    solver = problem.build_solver(de.timesteppers.RK443)
logger.info(f'Solver {stepper} built')

solver.stop_sim_time = stop_sim_time


restart_path = basepath / "restart/restart.h5"
if not restart_path.exists():
    # Start simulation from a static state; initial conditions are parallel
    # friendly
    logger.info("Restart path not found.")
    z = domain.grid(1)
    T = solver.state['T']
    Tz = solver.state['Tz']

    # Random perturbations, initialized globally for same results in parallel
    gshape = domain.dist.grid_layout.global_shape(scales=1)
    slices = domain.dist.grid_layout.slices(scales=1)
    rand = np.random.RandomState(seed=42)
    noise = rand.standard_normal(gshape)[slices] # pyright: ignore

    # Linear background + perturbations damped at walls
    zb, zt = z_basis.interval
    pert =  1e-3 * noise * (zt - z) * (z - zb)
    T['g'] = F*pert
    T.differentiate('z', out=Tz)

    # Timestepping and output
    dt = arg_dt
    fh_mode = 'overwrite'

else:
# restart simulation from a given initial condition, specified as an index in a .h5 file
    # Restart
    write, last_dt = solver.load_state(restart_path, restart_index)

    # Timestepping and output - if the last sim was using a smaller timestep,
    # it's a good idea to use that instead of the user specified.
    if cfl:
        dt = min(last_dt, arg_dt)
    else:
        dt = arg_dt
        
    fh_mode = 'append'


# File handling

# For movie making magic; use command line argument --snapshots to activate
if snapshots:
    snapshots_file = basepath / "snapshots"
    snapshots_file.mkdir(exist_ok=True)

    snapshots = solver.evaluator.add_file_handler(snapshots_file, iter=snapshots_iters, max_writes=1000, mode=fh_mode)
    snapshots.add_task("T-(z-1/2)", name = 'temp')
    snapshots.add_task("oy", name='vorticity')

# States saved as checkpoints for restarting. Can adjust iter as necessary.
# Infrequent saves are better when the simulation runs quickly (speed reasons).
state_file = basepath / 'state'
state_file.mkdir(exist_ok=True)
state = solver.evaluator.add_file_handler(state_file, iter=state_iters, max_writes=1000, mode=fh_mode)
state.add_system(solver.state)

# For Nu and Re calculations - data saved at 1 point in space FREQUENTLY
analysis = basepath / 'analysis'
analysis.mkdir(exist_ok=True)
analysis = solver.evaluator.add_file_handler(analysis, iter=analysis_iters, max_writes=1000, mode=fh_mode)

analysis.add_task('R/P', name='Pr')
analysis.add_task('1/(P*R)', name='Ra')

# For calculating Nu
analysis.add_task("kappa_xz * integ_z(integ_x( (w*(T+1/2-z))))", name='avg_wT')
analysis.add_task("kappa_xz * integ_z( integ_x( (dx(u)**2) + dz(u)**2 + dx(w)**2 + dz(w)**2 ))", name='avg_grad_u_sq')
analysis.add_task("kappa_xz * integ_z( integ_x( (dx(T+1/2-z))**2 + (dz(T+1/2-z)**2) ))", name='avg_grad_T_sq')
analysis.add_task("kappa_xz * integ_z( integ_x( (oy**2) ))", name='avg_oy_sq')

# For calculating Re and plotting profiles
analysis.add_task("kappa_xz * integ_z(integ_x( u**2 + w**2 ))", name='avg_K')
analysis.add_task("kappa_x*integ_x( T )", name='avg_T')
analysis.add_task("kappa_x*integ_x( u**2 )", name='avg_u_sq')
analysis.add_task("kappa_x*integ_x( w**2 )", name='avg_w_sq')
analysis.add_task("kappa_x*integ_x( T**2 )", name='avg_T_sq')


# For analysis - data saved at 1 point in space INFREQUENTLY
field_analysis_file = basepath / 'field_analysis'
field_analysis_file.mkdir(exist_ok=True)
field_analysis = solver.evaluator.add_file_handler(field_analysis_file, iter=field_analysis_iters, max_writes=1000, mode=fh_mode)

field_analysis.add_task('R/P', name='Pr')
field_analysis.add_task('1/(P*R)', name='Ra')

field_analysis.add_task("interp(u, x=0, z=0)", name='u')
field_analysis.add_task("interp(w, x=0, z=0)", name='w')
field_analysis.add_task("interp(T, x=0, z=0)", name='T')
field_analysis.add_task("interp(oy*w+1/2*dx(u**2+w**2), x=0, z=0)", name='u.grad_u')
field_analysis.add_task("interp(-oy*u+1/2*dz(u**2+w**2), x=0, z=0)", name='u.grad_w')
field_analysis.add_task("interp(u*dx(T)+w*Tz, x=0, z=0)", name='u.grad_T')
field_analysis.add_task("interp(dx(p)-1/2*dx(u**2+w**2), x=0, z=0)", name='p_x')
field_analysis.add_task("interp(dz(p)-1/2*dz(u**2+w**2), x=0, z=0)", name='p_z')
field_analysis.add_task("interp(d(u, x=2), x=0, z=0)", name='u_xx')
field_analysis.add_task("interp(d(u, z=2), x=0, z=0)", name='u_zz')
field_analysis.add_task("interp(d(w, x=2), x=0, z=0)", name='w_xx')
field_analysis.add_task("interp(d(w, z=2), x=0, z=0)", name='w_zz')
field_analysis.add_task("interp(d(T, x=2), x=0, z=0)", name='T_xx')
field_analysis.add_task("interp(dz(Tz), x=0, z=0)", name='T_zz')


# Main loop
try:
    logger.info('Starting loop.')
    logger.info("-" * 80)
    start_time = time.time()
    message_num_iters = 500
    start_iter = solver.iteration

    # Variable timestepping:
    if cfl:
        CFL = flow_tools.CFL(
            solver,
            initial_dt=dt,
            cadence=2,
            safety=0.15,
            max_dt=0.1,
            max_change=1.01,
            threshold=0.005
        )
        CFL.add_velocities(('u', 'w'))
        dt = CFL.compute_dt()
        dts = []
        while solver.proceed:
            solver.step(dt)
            if (solver.iteration-1) % message_num_iters == 0:
                l_dts = len(dts)
                if l_dts > 0:
                    logger.info(f"The timestep has changed {l_dts} time(s) over the last {message_num_iters} iterations.")
                    logger.info(f"The average timestep was {np.mean(dts):.4e}.") 
                    logger.info(f"The minimum was {np.min(dts):.4e}.")
                    logger.info(f"The maximum was {np.max(dts):.4e}.")
                elif (solver.iteration-1) == 0:
                    pass
                else:
                    logger.info(f"The timestep did not change over the last {message_num_iters} iterations.")
                
                logger.info('Iteration: %i, Time: %e, dt: %e' %(solver.iteration, solver.sim_time, dt))
                logger.info("-" * 80)
                
                
                dts = []
            new_dt = CFL.compute_dt()
            
            if new_dt != dt:
                dts.append(new_dt)
                dt = new_dt

    # Constant timesteps:
    else:
        while solver.proceed:
            solver.step(dt)
            if (solver.iteration-1) % message_num_iters == 0:
                logger.info('Iteration: %i, Time: %e, dt: %e' %(solver.iteration, solver.sim_time, dt))
except:
    logger.error('Exception raised, triggering end of main loop.')
    raise
finally:
    # Total time and iterations
    end_time = time.time()
    end_iter = solver.iteration
    total_iter = end_iter-start_iter
    total_time = end_time-start_time

    # Degrees of freedom for speed calculation
    dof = nx * nz * len(problem.variables)

    # Write info to job_sim.out file
    logger.info('Iterations: %i' %solver.iteration)
    logger.info('Sim end time: %f' %solver.sim_time)
    logger.info('Run time: %.2f sec' %(total_time))
    logger.info('Run time: %f cpu-hr' %((total_time)/60/60*domain.dist.comm_cart.size))
    logger.info('Speed: %f' %( (dof * total_iter) / (domain.dist.comm_cart.size * total_time)))
