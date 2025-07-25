"""
2D Rayleigh Benard convection using Dedalus. The equations of motion are non-dimensionalized
using the buoyancy time. Boundary conditions are no slip in the z-direction and periodic in the x direction.
Using the '--snapshots' option, temperature and vorticity can be saved at regular intervals for visualization
Various (many) other quantities are saved for analysis. This script supports restarts, and integrates using a
variable timestep.

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


args = docopt(__doc__)

# Parameters
Lx, Lz = (int(args['--Gamma']), 1) # fix non-dimensional height of the layer at 1
Pr = 10**np.float64(args['--Pr_exp'])
Ra = np.float64(args['--Ra'])
nx, nz = (Lx*int(args['--res']), Lz*int(args['--res']))
stepper = str(args['--stepper'])
flow_bc = "no-slip"
basepath = Path(str(args['--basepath']))
arg_dt = np.float64(args['--dt']) #use constant dt
stop_sim_time = np.float64(args['--sim_time'])
snapshots = args['--snapshots']
cfl = args['--cfl']
restart_index = int(args['--index'])


# Create bases and domain
x_basis = de.Fourier('x', nx, interval=(0, Lx), dealias=3/2)
z_basis = de.Chebyshev('z', nz, interval=(-Lz/2, Lz/2), dealias=3/2)
domain = de.Domain([x_basis, z_basis], grid_dtype=np.float64)

n=10
times = np.ones((2,n))
for i in range(2*n):
    if i%2 == 0: # Fast version (hopefully)
        start_time = time.time_ns()
        # Problem setup: 2D Boussinesq equations; non-dimensionalized by the buoyancy time
        # and written in terms of the vorticity.
        problem = de.IVP(domain, variables=['p','T','u','w','Tz','oy'])

        problem.meta['p','T','u','w']['z']['dirichlet'] = True


        # Parameters
        problem.parameters['P'] = (Ra * Pr)**(-1/2)
        problem.parameters['R'] = (Ra / Pr)**(-1/2)
        problem.parameters['F'] = F = 1
        problem.parameters['kappa_xz'] = 1/(Lx*Lz)
        problem.parameters['kappa_x'] = 1/Lx

        # Equations describing the PDE
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
        logger.info('Solver built')

        restart_path = basepath / "restart/restart.h5"
        if not restart_path.exists():
            # start simulation from a static state; initial conditions are parallel friendly
            # Initial conditions
            logger.info("Restart path not found.")
            x = domain.grid(0)
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

            # Timestepping and output
            if cfl:
                dt = min(last_dt, arg_dt)
            else:
                dt = arg_dt
                
            fh_mode = 'append'
        # Integration parameters
        solver.stop_sim_time = stop_sim_time

        # For movie making magic; use command line argument --snapshots to activate
        if snapshots:
            snapshots_file = basepath / "snapshots"
            snapshots_file.mkdir(exist_ok=True)

            snapshots = solver.evaluator.add_file_handler(snapshots_file, iter=200, max_writes=1000, mode=fh_mode)
            snapshots.add_task("T-(z-1/2)", name = 'temp')
            snapshots.add_task("oy", name='vorticity')

        # states saved as checkpoints for restarting. Can adjust iter as necessary.
        state_file = basepath / 'state'
        state_file.mkdir(exist_ok=True)
        state = solver.evaluator.add_file_handler(state_file, iter=5000, max_writes=100, mode=fh_mode)
        state.add_system(solver.state)

        # For Nu and Re calculations - data saved at 1 point in space FREQUENTLY
        nu_analysis = basepath / 'nu_analysis'
        nu_analysis.mkdir(exist_ok=True)
        nu_analysis = solver.evaluator.add_file_handler(nu_analysis, iter=10, max_writes=1000, mode=fh_mode)
        # easiest way to save Ra and Pr; takes up minimal space
        nu_analysis.add_task('R/P', name='Pr')
        nu_analysis.add_task('1/(P*R)', name='Ra')

        # For calculating Nu
        nu_analysis.add_task("kappa_xz * integ_z(integ_x( (w*(T+1/2-z))))", name='avg_wT')
        nu_analysis.add_task("kappa_xz * integ_z( integ_x( (dx(u)**2) + dz(u)**2 + dx(w)**2 + dz(w)**2 ))", name='avg_grad_u_sq')
        nu_analysis.add_task("kappa_xz * integ_z( integ_x( (dx(T+1/2-z))**2 + (dz(T+1/2-z)**2) ))", name='avg_grad_T_sq')
        nu_analysis.add_task("kappa_xz * integ_z( integ_x( (oy**2) ))", name='avg_oy_sq')
        nu_analysis.add_task("kappa_xz * integ_z( integ_x( (dx(T))**2 + (Tz**2) ))", name='avg_grad_theta_sq')
        nu_analysis.add_task("kappa_xz * integ_z(integ_x( w*T ))", name='avg_wtheta')

        # For calculating Re
        nu_analysis.add_task("kappa_xz * integ_z(integ_x( u**2 + w**2 ))", name='avg_K')
        nu_analysis.add_task("kappa_x*integ_x( T )", name='avg_T')
        nu_analysis.add_task("kappa_x*integ_x( u**2 )", name='avg_u_sq')
        nu_analysis.add_task("kappa_x*integ_x( w**2 )", name='avg_w_sq')
        nu_analysis.add_task("kappa_x*integ_x( T**2 )", name='avg_T_sq')


        # For analysis - data saved at 1 point in space INFREQUENTLY
        analysis_file = basepath / 'analysis'
        analysis_file.mkdir(exist_ok=True)
        analysis = solver.evaluator.add_file_handler(analysis_file, iter=200, max_writes=1000, mode=fh_mode)
        # easiest way to save Ra and Pr; takes up minimal space
        analysis.add_task('R/P', name='Pr')
        analysis.add_task('1/(P*R)', name='Ra')

        # State and other field quantities
        analysis.add_task("interp(u, x=0, z=0)", name='u')
        analysis.add_task("interp(w, x=0, z=0)", name='w')
        analysis.add_task("interp(T, x=0, z=0)", name='T')
        analysis.add_task("interp(oy*w+1/2*dx(u**2+w**2), x=0, z=0)", name='u.grad_u')
        analysis.add_task("interp(-oy*u+1/2*dz(u**2+w**2), x=0, z=0)", name='u.grad_w')
        analysis.add_task("interp(u*dx(T)+w*Tz, x=0, z=0)", name='u.grad_T')
        analysis.add_task("interp(dx(p)-1/2*dx(u**2+w**2), x=0, z=0)", name='p_x')
        analysis.add_task("interp(dz(p)-1/2*dz(u**2+w**2), x=0, z=0)", name='p_z')
        analysis.add_task("interp(d(u, x=2), x=0, z=0)", name='u_xx')
        analysis.add_task("interp(d(u, z=2), x=0, z=0)", name='u_zz')
        analysis.add_task("interp(d(w, x=2), x=0, z=0)", name='w_xx')
        analysis.add_task("interp(d(w, z=2), x=0, z=0)", name='w_zz')
        analysis.add_task("interp(d(T, x=2), x=0, z=0)", name='T_xx')
        analysis.add_task("interp(dz(Tz), x=0, z=0)", name='T_zz')


        # Main loop
        try:
            logger.info('Starting loop')
            # variable timestepping
            if cfl:
                CFL = flow_tools.CFL(
                    solver,
                    initial_dt=dt,
                    cadence=10,
                    safety=0.75,
                    max_dt=0.01,
                    max_change=1.01,
                    threshold=0.05
                )
                CFL.add_velocities(('u', 'w'))
                dt = CFL.compute_dt()
                while solver.proceed:
                    solver.step(dt)
                    if (solver.iteration-1) % 5000 == 0:
                        logger.info('Iteration: %i, Time: %e, dt: %e' %(solver.iteration, solver.sim_time, dt))
                    new_dt = CFL.compute_dt()
                    if new_dt != dt:
                        dt = new_dt
                        logger.info('Iteration: %i, Time: %e, NEW dt: %e' %(solver.iteration, solver.sim_time, dt))

            else:
                while solver.proceed:
                    solver.step(dt)
                    if (solver.iteration-1) % 5000 == 0:
                        logger.info('Iteration: %i, Time: %e, dt: %e' %(solver.iteration, solver.sim_time, dt))
        except:
            logger.error('Exception raised, triggering end of main loop.')
            raise
        finally:
            end_time = time.time_ns()
            run_time = (end_time-start_time) * 1e-9
            logger.info('Iterations: %i' %solver.iteration)
            logger.info('Sim end time: %f' %solver.sim_time)
            logger.info(f'Run time: {run_time:.5f} sec')
            logger.info('Run time: %f cpu-hr' %((run_time)/60/60*domain.dist.comm_cart.size))
        
    else:
        start_time = time.time_ns()
        problem = de.IVP(domain, variables=['p','T','u','w','Tz','oy','avg_wT',
       'avg_K','avg_T','avg_T_sq','avg_u_sq','avg_w_sq','avg_oy_sq',
       'avg_grad_T_sq', 'avg_grad_u_sq', 'avg_wtheta', 'avg_grad_theta_sq'])

        problem.meta['p','T','u','w']['z']['dirichlet'] = True

        problem.meta['avg_T','avg_u_sq','avg_w_sq','avg_T_sq','avg_wT','avg_K',
                    'avg_grad_T_sq','avg_oy_sq', 'avg_grad_u_sq','avg_wtheta',
                    'avg_grad_theta_sq']['x']['constant'] = True

        problem.parameters['P'] = (Ra * Pr)**(-1/2)
        problem.parameters['R'] = (Ra / Pr)**(-1/2)
        problem.parameters['F'] = F = 1
        problem.parameters['kappa_xz'] = 1/(Lx*Lz)
        problem.parameters['kappa_x'] = 1/Lx

        problem.add_equation("dx(u) + dz(w) = 0")
        problem.add_equation("dt(T) - P*(d(T,x=2) + dz(Tz)) - F*w  = -(u*dx(T) + w*Tz)")
        problem.add_equation("dt(u) - R*dz(oy) + dx(p)             = -oy*w")
        problem.add_equation("dt(w) + R*dx(oy) + dz(p) - T         = oy*u")
        problem.add_equation("Tz - dz(T) = 0")
        problem.add_equation("oy + dx(w) - dz(u) = 0")

        # time integrals of various quantities
        # for calculating Nu
        problem.add_equation("dt(avg_wT)            = kappa_xz*integ_z(integ_x( w*(T+1/2-z) ))")
        problem.add_equation("dt(avg_wtheta)        = kappa_xz*integ_z(integ_x( w*T ))")
        problem.add_equation("dt(avg_oy_sq)         = kappa_xz*integ_z(integ_x( (oy**2) ))")
        problem.add_equation("dt(avg_grad_u_sq)     = kappa_xz*integ_z(integ_x( (dx(u)**2) + dz(u)**2 + dx(w)**2 + dz(w)**2 ))")
        problem.add_equation("dt(avg_grad_T_sq)     = kappa_xz*integ_z(integ_x( (dx(T+1/2-z))**2 + (dz(T+1/2-z)**2) ))")
        problem.add_equation("dt(avg_grad_theta_sq) = kappa_xz*integ_z(integ_x( (dx(T))**2 + (Tz**2) ))")

        # for calculating Re
        problem.add_equation("dt(avg_K)    = kappa_xz*integ_z(integ_x( u**2 + w**2 ))")
        problem.add_equation("dt(avg_T)    = kappa_x*integ_x( T )")
        problem.add_equation("dt(avg_u_sq) = kappa_x*integ_x( u**2 )")
        problem.add_equation("dt(avg_w_sq) = kappa_x*integ_x( w**2 )")
        problem.add_equation("dt(avg_T_sq) = kappa_x*integ_x( T**2 )")


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
        logger.info('Solver built')

        restart_path = basepath / "restart/restart.h5"
        if not restart_path.exists():
            # start simulation from a static state; initial conditions are parallel friendly
            # Initial conditions
            logger.info("Restart path not found.")
            x = domain.grid(0)
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

            # Timestepping and output
            if cfl:
                dt = min(last_dt, arg_dt)
            else:
                dt = arg_dt
                
            fh_mode = 'append'
        # Integration parameters
        solver.stop_sim_time = stop_sim_time

        # For movie making magic; use command line argument --snapshots to activate
        if snapshots:
            snapshots_file = basepath / "snapshots"
            snapshots_file.mkdir(exist_ok=True)

            snapshots = solver.evaluator.add_file_handler(snapshots_file, iter=200, max_writes=1000, mode=fh_mode)
            snapshots.add_task("T-(z-1/2)", name = 'temp')
            snapshots.add_task("oy", name='vorticity')

        # states saved as checkpoints for restarting. Can adjust iter as necessary.
        state_file = basepath / 'state'
        state_file.mkdir(exist_ok=True)
        state = solver.evaluator.add_file_handler(state_file, iter=5000, max_writes=100, mode=fh_mode)
        state.add_system(solver.state)

        # For analysis - data saved at 1 point in space
        analysis_file = basepath / 'analysis'
        analysis_file.mkdir(exist_ok=True)
        analysis = solver.evaluator.add_file_handler(analysis_file, iter=100, max_writes=1000, mode=fh_mode)
        # easiest way to save Ra and Pr; takes up minimal space
        analysis.add_task('R/P', name='Pr')
        analysis.add_task('1/(P*R)', name='Ra')

        analysis.add_task("interp(u, x=0, z=0)", name='u')
        analysis.add_task("interp(w, x=0, z=0)", name='w')
        analysis.add_task("interp(T, x=0, z=0)", name='T')
        analysis.add_task("interp(oy*w+1/2*dx(u**2+w**2), x=0, z=0)", name='u.grad_u')
        analysis.add_task("interp(-oy*u+1/2*dz(u**2+w**2), x=0, z=0)", name='u.grad_w')
        analysis.add_task("interp(u*dx(T)+w*Tz, x=0, z=0)", name='u.grad_T')
        analysis.add_task("interp(dx(p)-1/2*dx(u**2+w**2), x=0, z=0)", name='p_x')
        analysis.add_task("interp(dz(p)-1/2*dz(u**2+w**2), x=0, z=0)", name='p_z')
        analysis.add_task("interp(d(u, x=2), x=0, z=0)", name='u_xx')
        analysis.add_task("interp(d(u, z=2), x=0, z=0)", name='u_zz')
        analysis.add_task("interp(d(w, x=2), x=0, z=0)", name='w_xx')
        analysis.add_task("interp(d(w, z=2), x=0, z=0)", name='w_zz')
        analysis.add_task("interp(d(T, x=2), x=0, z=0)", name='T_xx')
        analysis.add_task("interp(dz(Tz), x=0, z=0)", name='T_zz')
        analysis.add_task("interp(avg_wT, x=0, z=0)", name='avg_wT')
        analysis.add_task("interp(avg_K, x=0, z=0)", name='avg_K')
        analysis.add_task("interp(avg_T, x=0)", name='avg_T')
        analysis.add_task("interp(avg_u_sq, x=0)", name='avg_u_sq')
        analysis.add_task("interp(avg_w_sq, x=0)", name='avg_w_sq')
        analysis.add_task("interp(avg_u_sq + avg_w_sq, x=0)", name='h_avg_K')
        analysis.add_task("interp(avg_T_sq, x=0)", name='avg_T_sq')
        analysis.add_task("interp(avg_grad_u_sq, x=0, z=0)", name='avg_grad_u_sq')
        analysis.add_task("interp(avg_grad_T_sq, x=0, z=0)", name='avg_grad_T_sq')
        analysis.add_task("interp(avg_oy_sq, x=0, z=0)", name='avg_oy_sq')
        analysis.add_task("interp(avg_grad_theta_sq, x=0, z=0)", name='avg_grad_theta_sq')
        analysis.add_task("interp(avg_wtheta, x=0, z=0)", name='avg_wtheta')



        # Main loop
        try:
            logger.info('Starting loop')
            # variable timestepping
            if cfl:
                CFL = flow_tools.CFL(
                    solver,
                    initial_dt=dt,
                    cadence=10,
                    safety=0.75,
                    max_dt=0.01,
                    max_change=1.01,
                    threshold=0.05
                )
                CFL.add_velocities(('u', 'w'))
                dt = CFL.compute_dt()
                while solver.proceed:
                    solver.step(dt)
                    if (solver.iteration-1) % 5000 == 0:
                        logger.info('Iteration: %i, Time: %e, dt: %e' %(solver.iteration, solver.sim_time, dt))
                    new_dt = CFL.compute_dt()
                    if new_dt != dt:
                        dt = new_dt
                        logger.info('Iteration: %i, Time: %e, NEW dt: %e' %(solver.iteration, solver.sim_time, dt))

            else:
                while solver.proceed:
                    solver.step(dt)
                    if (solver.iteration-1) % 5000 == 0:
                        logger.info('Iteration: %i, Time: %e, dt: %e' %(solver.iteration, solver.sim_time, dt))
        except:
            logger.error('Exception raised, triggering end of main loop.')
            raise
        finally:
            end_time = time.time_ns()
            run_time = (end_time-start_time) * 1e-9
            logger.info('Iterations: %i' %solver.iteration)
            logger.info('Sim end time: %f' %solver.sim_time)
            logger.info(f'Run time: {run_time:.5f} sec')
            logger.info('Run time: %f cpu-hr' %((run_time)/60/60*domain.dist.comm_cart.size))
        
    times[i%2, int(i/2)] = run_time
import matplotlib.pyplot as plt

plt.plot(times[0,:], 'ro', label="saving instantaneous")
plt.plot(times[1,:], 'bo', label='saving cumulative')
plt.legend()
plt.xlabel("Run Number")
plt.ylabel("Seconds to run")
plt.hlines(y=np.mean(times, axis=-1), xmin=0, xmax=n+1, colors=['r', 'b'], linestyles='dashed')
plt.savefig('times.png')