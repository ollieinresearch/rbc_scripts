"""
3D Rayleigh Benard convection using Dedalus. The equations of motion are 
non-dimensionalized using the buoyancy time. Boundary conditions are no slip in 
the z-direction and periodic in the x and y directions. Using the '--snapshots' 
option, temperature and vorticity can be saved at regular intervals for 
visualization. Various (many) other quantities are saved for analysis. This 
script supports restarts, and integrates using a fixed timestep.

Usage:
    rayleigh_benard_script3d.py [options]

Options:
    --Ra=<Ra>                Rayleigh number
    --Pr_exp=<Pr>            Exponent on the Prandtl number (base 10)
    --res=<res>              Resolution in modes/unit space
    --dt=<dt>                Timestep
    --sim_time=<time>        Simulation time
    --basepath=<path>        basepath for output files
    --stepper=<stepper>      The solver to be used for timestepping [default: RK443]
    --Lx=<Lx>                The length of the cell in the x-axis [default: 2]
    --Ly=<Ly>                The length of the cell in the y-axis [default: 2]
    --meshx=<size>           The size of the x axis of the 2D processor mesh [default: 16]
    --meshy=<size>           The size of the y axis of the 2D processor mesh [default: 16]
    --index=<index>          The index of the restart file to begin simulations from [default: -1]
    --snapshots              Flag to activate the snapshots filehandler that saves the temperature for visualisations.
    --cfl                    Flag to activate CFL timestepping.

"""

#TODO: Make iter params 
from docopt import docopt
import numpy as np
from mpi4py import MPI # pyright: ignore
import time
from pathlib import Path

from dedalus import public as de
from dedalus.extras import flow_tools

import logging
logger = logging.getLogger(__name__)


args = docopt(__doc__)

# Parameters
# fix non-dimensional height of the layer at 1
Lx, Ly, Lz = (int(args['--Lx']), int(args['--Ly']), 1)
Pr = 10**np.float64(args['--Pr_exp'])
Ra = np.float64(args['--Ra'])
nx, ny, nz = (Lx*int(args['--res']), Ly*int(args['--res']), Lz*int(args['--res']))
stepper = str(args['--stepper'])
flow_bc = "no-slip"
meshx, meshy = int(args['--meshx']), int(args['--meshy'])
basepath = Path(str(args['--basepath']))
arg_dt = np.float64(args['--dt'])
stop_sim_time = np.float64(args['--sim_time'])
snapshots = args['--snapshots']
cfl = args['--cfl']
restart_index = int(args['--index'])

# Create bases and domain
x_basis = de.Fourier('x', nx, interval=(0, Lx), dealias=3/2)
y_basis = de.Fourier('y', ny, interval=(0, Ly), dealias=3/2)
z_basis = de.Chebyshev('z', nz, interval=(-Lz/2, Lz/2), dealias=3/2)

bases = [x_basis, y_basis, z_basis]
meshes = [meshx, meshy]
domain = de.Domain(bases, grid_dtype=np.float64, mesh=meshes)

# 3D Boussinesq dynamics
variables = ['p','T','u','v','w','Tz','uz','vz','wz']
problem = de.IVP(
    domain,
    variables=variables
)

problem.meta['p','T','u','v','w']['z']['dirichlet'] = True


problem.parameters['P'] = (Ra * Pr)**(-1/2)
problem.parameters['R'] = (Ra / Pr)**(-1/2)
problem.parameters['F'] = F = 1
problem.parameters['kappa_xyz'] = 1/(Lx*Ly*Lz)
problem.parameters['kappa_xy'] = 1/(Lx*Ly)

# Divergence free
problem.add_equation("dx(u) + dy(v) + wz = 0")

# Temperature PDE
problem.add_equation("dt(T) - P*(d(T,x=2) + d(T,y=2) + dz(Tz)) - F*w       = -(u*dx(T) + v*dy(T) + w*Tz)")

# xyz directions NSE
problem.add_equation("dt(u) - R*(d(u,x=2) + d(u,y=2) + dz(uz)) + dx(p)     = -(u*dx(u) + v*dy(u) + w*uz)")
problem.add_equation("dt(v) - R*(d(v,x=2) + d(v,y=2) + dz(vz)) + dy(p)     = -(u*dx(v) + v*dy(v) + w*vz)")
problem.add_equation("dt(w) - R*(d(w,x=2) + d(w,y=2) + dz(wz)) + dz(p) - T = -(u*dx(w) + v*dy(w) + w*wz)")

# Order reduction
problem.add_equation("Tz - dz(T) = 0")
problem.add_equation("uz - dz(u) = 0")
problem.add_equation("vz - dz(v) = 0")
problem.add_equation("wz - dz(w) = 0")


# Boundary conditions
problem.add_bc("left(T) = 0")
problem.add_bc("right(T) = 0")

if flow_bc == "no-slip":
    problem.add_bc("left(u) = 0")
    problem.add_bc("right(u) = 0")
    problem.add_bc("left(v) = 0")
    problem.add_bc("right(v) = 0")

if flow_bc == "free-slip":
    problem.add_bc("left(uz) = 0")
    problem.add_bc("right(uz) = 0")
    problem.add_bc("left(vz) = 0")
    problem.add_bc("right(vz) = 0")

problem.add_bc("left(w) = 0")
problem.add_bc("right(w) = 0", condition="(nx != 0) or (ny != 0)")
problem.add_bc("right(p) = 0", condition="(nx == 0) and (ny == 0)")

# Build solver; time stepper passed as parameter
solver=None
if stepper == 'RK222':
    solver = problem.build_solver(de.timesteppers.RK222)
elif stepper == 'CNAB2':
    solver = problem.build_solver(de.timesteppers.CNAB2)
elif stepper == 'MCNAB2':
    solver = problem.build_solver(de.timesteppers.MCNAB2)
elif stepper == 'SBDF4':
    solver = problem.build_solver(de.timesteppers.SBDF4)
elif stepper == 'SBDF2':
    solver = problem.build_solver(de.timesteppers.SBDF2)
else:
    solver = problem.build_solver(de.timesteppers.RK443)
logger.info(f'Solver {stepper} built.')

restart_path = basepath / "restart/restart.h5"
if not restart_path.exists():
    # start from conductive state with random perturbations; parallel friendly
    # Initial conditions
    z = domain.grid(1)
    T = solver.state['T']
    Tz = solver.state['Tz']

    # Random perturbations, initialized globally for same results in parallel
    gshape = domain.dist.grid_layout.global_shape(scales=1)
    slices = domain.dist.grid_layout.slices(scales=1)
    rand = np.random.RandomState(seed=42)
    noise = rand.standard_normal(gshape)[slices]

    # Linear background + perturbations damped at walls
    zb, zt = z_basis.interval
    pert =  1e-3 * noise * (zt - z) * (z - zb)
    T['g'] = F*pert
    T.differentiate('z', out=Tz)

    fh_mode = 'overwrite'

else:
# Start simulation from initial conditions specified by index in a .h5 file
    # Restart
    write, last_dt = solver.load_state(restart_path, restart_index)

    fh_mode = 'append'

# Integration parameters
solver.stop_sim_time = stop_sim_time

# For movie making magic; use command line argument --snapshots to activate
# TODO: figure out other scalar field for 3d visualization
if snapshots:
    snapshots_file = basepath / 'snapshots'
    snapshots_file.mkdir(exist_ok=True)

    snapshots = solver.evaluator.add_file_handler(snapshots_file, iter=200, max_writes=50, mode=fh_mode)
    snapshots.add_task("T-(z-1/2)", scales=2, name = 'temp')
    snapshots.add_task("w", scales=4, name = 'w')
    

# states saved as checkpoints for restarting. Can adjust iter as necessary.
state_file = basepath / 'state'
state_file.mkdir(exist_ok=True)
state = solver.evaluator.add_file_handler(state_file, iter=1000, max_writes=25, mode=fh_mode)
state.add_system(solver.state)

# For field_analysis - data saved at 1 point in space INFREQUENTLY
field_analysis_file = basepath / 'field_analysis'
field_analysis_file.mkdir(exist_ok=True)
field_analysis = solver.evaluator.add_file_handler(field_analysis_file, iter=100, max_writes=1000, mode=fh_mode)

field_analysis.add_task('R/P', name='Pr')
field_analysis.add_task('1/(P*R)', name='Ra')
field_analysis.add_task("interp(u, x=0, y=0, z=0)", name='u')
field_analysis.add_task("interp(v, x=0, y=0, z=0)", name='v')
field_analysis.add_task("interp(w, x=0, y=0, z=0)", name='w')
field_analysis.add_task("interp(T, x=0, y=0, z=0)", name='T')
field_analysis.add_task("interp(u*dx(u)+v*dy(u)+w*uz, x=0, y=0, z=0)", name='u.grad_u')
field_analysis.add_task("interp(u*dx(v)+v*dy(v)+w*vz, x=0, y=0, z=0)", name='u.grad_v')
field_analysis.add_task("interp(u*dx(w)+v*dy(w)+w*wz, x=0, y=0, z=0)", name='u.grad_w')
field_analysis.add_task("interp(u*dx(T)+v*dy(T)+w*Tz, x=0, y=0, z=0)", name='u.grad_T')
field_analysis.add_task("interp(dx(p), x=0, y=0, z=0)", name='p_x')
field_analysis.add_task("interp(dz(p), x=0, y=0, z=0)", name='p_z')
field_analysis.add_task("interp(d(u, x=2), x=0, y=0, z=0)", name='u_xx')
field_analysis.add_task("interp(d(u, y=2), x=0, y=0, z=0)", name='u_yy')
field_analysis.add_task("interp(dz(uz), x=0, y=0, z=0)", name='u_zz')
field_analysis.add_task("interp(d(v, x=2), x=0, y=0, z=0)", name='v_xx')
field_analysis.add_task("interp(d(v, y=2), x=0, y=0, z=0)", name='v_yy')
field_analysis.add_task("interp(dz(vz), x=0, y=0, z=0)", name='v_zz')
field_analysis.add_task("interp(d(w, x=2), x=0, y=0, z=0)", name='w_xx')
field_analysis.add_task("interp(d(w, y=2), x=0, y=0, z=0)", name='w_yy')
field_analysis.add_task("interp(dz(wz), x=0, y=0, z=0)", name='w_zz')
field_analysis.add_task("interp(d(T, x=2), x=0, y=0, z=0)", name='T_xx')
field_analysis.add_task("interp(d(T, y=2), x=0, y=0, z=0)", name='T_yy')
field_analysis.add_task("interp(dz(Tz), x=0, y=0, z=0)", name='T_zz')


# For calculating Nu and Re - data saved at 1 point in space FREQUENTLY
analysis_file = basepath / 'analysis'
analysis_file.mkdir(exist_ok=True)
analysis = solver.evaluator.add_file_handler(analysis_file, iter=50, max_writes=1000, mode=fh_mode)
analysis.add_task('R/P', name='Pr')
analysis.add_task('1/(P*R)', name='Ra')

# For Re
analysis.add_task("kappa_xyz*integ_z(integ_y(integ_x(u**2+v**2+w**2)))", name='avg_K')
analysis.add_task("kappa_xy*integ_y(integ_x(T))", name='avg_T')
analysis.add_task("kappa_xy*integ_y(integ_x(u**2))", name='avg_u_sq')
analysis.add_task("kappa_xy*integ_y(integ_x(v**2))", name='avg_v_sq')
analysis.add_task("kappa_xy*integ_y(integ_x(w**2))", name='avg_w_sq')
analysis.add_task("kappa_xy*integ_y(integ_x(T**2))", name='avg_T_sq')

# For Nu - 3 methods to compare
analysis.add_task("kappa_xyz*integ_z(integ_y(integ_x(w*(T+1/2-z))))", name='avg_wT')
analysis.add_task("kappa_xyz*integ_z(integ_y(integ_x( (dx(u)**2) + dy(u)**2 + dz(u)**2 + dx(v)**2 + dy(v)**2 + dz(v)**2 + dx(w)**2 + dy(w)**2 + dz(w)**2 )))", name='avg_grad_u_sq')
analysis.add_task("kappa_xyz*integ_z(integ_y(integ_x( (dx(T+1/2-z))**2 + (dy(T+1/2-z)**2) + (dz(T+1/2-z)**2) )))", name='avg_grad_T_sq')


# Main loop
try:
    logger.info('Starting loop.')
    logger.info("-" * 80)
    start_time = time.time()

    # How often to display the message about time(steps)
    message_num_iters = 500
    start_iter = solver.iteration
    if cfl:
        CFL = flow_tools.CFL(
            solver,
            initial_dt=arg_dt,
            cadence=25,
            safety=0.75,
            max_dt=0.1,
            max_change=1.5,
            threshold=0.2
        )
        CFL.add_velocities(('u', 'v', 'w'))
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

    else:
        # Not using CFL
        dt = arg_dt
        while solver.proceed:
            solver.step(dt)
            if (solver.iteration-1) % message_num_iters == 0:
                logger.info('Iteration: %i, Time: %e, dt: %e' %(solver.iteration, solver.sim_time, dt))
                logger.info("-" * 80)

except:
    logger.error('Exception raised, triggering end of main loop.')
    raise
finally:
    end_time = time.time()
    end_iter = solver.iteration

    total_iter = end_iter-start_iter
    total_time = end_time-start_time
    dof = nx * ny * nz * len(problem.variables)

    logger.info(f'Total iterations: {total_iter}')
    logger.info('Sim end time: %f' %solver.sim_time)
    logger.info('Run time: %.2f sec' %(total_time))
    logger.info('Run time: %f cpu-hr' %((total_time)/60/60*domain.dist.comm_cart.size))
    logger.info('Speed: %f' %( (dof * total_iter) / (domain.dist.comm_cart.size * total_time)))
