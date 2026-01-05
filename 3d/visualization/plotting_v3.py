"""
Script to perform analysis tasks on data obtained from running a 3D
Rayleigh-Benard convection simulation. Accepts an analysis file as a command
line argument and takes time averages from the specified time. All output is
written to the folders 'outputs' and 'preliminary_outputs'. Outputs plots of 
the Nusselt number, various profiles and information about the simulation to a 
text document.

Usage:
    isos_test.py <files>... [--basepath=<dir>] [--max_vort=<max_vort>] [--max_vert=<max_vert>] [--nu=<nu>]
    isos_test.py <files>...

Options:
    --basepath=<dir>  Path to parent folder for output [default: ./analysis]
    --max_vort=<max_vort>  Maximum vorticity [default: 100]
    --max_vert=<max_vert>  Maximum vertical velocity [default: 100]
    --nu=<nu> Nusselt from this sim for boundary layer [default: 5]
"""

#TODO: choose consistent scales and then reduce number of grids.

import h5py
import numpy as np
import pyvista as pv
import os
from pathlib import Path
from matplotlib.scale import AsinhTransform

# For vorticity scaling
asinh = AsinhTransform(linear_width=1.5)

def main(h5_file, start, count, output_dir, mvort, mvert, nu):
    #mvort = asinh.transform(mvort)
    #mvert = asinh.transform(mvert)
    bl = (4.0*nu)
    # Load temp data and grid
    pv.start_xvfb()
    pv.global_theme.allow_empty_mesh = True
    with h5py.File(h5_file, 'r') as f:
        nt, nTx, nTy, nTz = f['tasks']['temperature'].shape
        _, nωx, nωy, nωz = f['tasks']['vorticity'].shape
        _, nwx, nwy, nwz = f['tasks']['w'].shape
        domain_sizes = (2.0, 2.0, 1.0)
        # Center of domain
        xmid = domain_sizes[0] / 2.0
        ymid = domain_sizes[1] / 2.0
        zmid = domain_sizes[2] / 2.0

        T = np.array(f['tasks']['temperature'][start:start+count], dtype=np.float64)
        T_1 = np.array(f['tasks']['temperature'][start:start+count, :, :, int(nTz/2)], dtype=np.float64)
        T_2 = np.array(f['tasks']['temperature'][start:start+count, :, :, int(nTz/(bl))], dtype=np.float64)

        ω = np.array(f['tasks']['vorticity'][start:start+count], dtype=np.float64)
        w = np.array(f['tasks']['w'][start:start+count], dtype=np.float64)

        w_1 = np.array(f['tasks']['w'][start:start+count, :, :, int(nwz/2)], dtype=np.float64)
        w_2 = np.array(f['tasks']['w'][start:start+count, :, :, int(nwz/(bl))], dtype=np.float64)

        Tdx = domain_sizes[0] / (nTx - 1)
        Tdy = domain_sizes[1] / (nTy - 1)
        Tdz = domain_sizes[2] / (nTz - 1)

        ωdx = domain_sizes[0] / (nωx - 1)
        ωdy = domain_sizes[1] / (nωy - 1)
        ωdz = domain_sizes[2] / (nωz - 1)

        wdx = domain_sizes[0] / (nwx - 1)
        wdy = domain_sizes[1] / (nwy - 1)
        wdz = domain_sizes[2] / (nwz - 1)
        
        
        origin = (0.0, 0.0, 0.0)
        
        # Rotating view radius
        rot_radius = 4.5

        # Ensure output directory exists
        os.makedirs(output_dir, exist_ok=True)

        # Loop through each time step
        for t in range(0, count):
            temp = T[t]
            #temp_1 = T_1[t]
            #temp_2 = T_2[t]

            vort = ω[t]
            vert = w[t]

            #vert_1 = w_1[t]
            #vert_2 = w_2[t]
            
            flat_temp = temp.flatten(order='F')
            flat_temp_1 = temp[:,:,int(nTz/2)].flatten(order='F')
            flat_temp_2 = temp[:,:,int(nTz/bl)].flatten(order='F')

            flat_vort = vort.flatten(order='F')
            flat_vert = vert.flatten(order='F')
            flat_vert_1 = vert[:,:,int(nwz/2)].flatten(order='F')
            flat_vert_2 = vert[:,:,int(nwz/bl)].flatten(order='F')

            # Create uniform grid
            t_grid = pv.ImageData(
                dimensions=(nTx, nTy, nTz),
                spacing=(Tdx, Tdy, Tdz),
                origin=origin
            )

            t_grid_1 = pv.Plane(
                center=origin,
                direction=(0, 0, 1),
                i_size=2.0,
                j_size=2.0,
                i_resolution=nTx-1,
                j_resolution=nTy-1
            )
            
            t_grid_2 = pv.Plane(
                center=origin,
                direction=(0, 0, 1),
                i_size=2.0,
                j_size=2.0,
                i_resolution=nTx-1,
                j_resolution=nTy-1
            )

            vort_grid = pv.ImageData(
                dimensions=(nωx, nωy, nωz),
                spacing=(ωdx, ωdy, ωdz),
                origin=origin
            )

            vert_grid = pv.ImageData(
                dimensions=(nwx, nwy, nwz),
                spacing=(wdx, wdy, wdz),
                origin=origin
            )

            vert_grid_1 = pv.Plane(
                center=origin,
                direction=(0, 0, 1),
                i_size=2.0,
                j_size=2.0,
                i_resolution=nwx-1,
                j_resolution=nwy-1
            )

            vert_grid_2 = pv.Plane(
                center=origin,
                direction=(0, 0, 1),
                i_size=2.0,
                j_size=2.0,
                i_resolution=nwx-1,
                j_resolution=nwy-1
            )

            t_grid.point_data['temp'] = flat_temp
            t_grid_1.point_data['temp'] = flat_temp_1
            t_grid_2.point_data['temp'] = flat_temp_2
            
            vort_grid.point_data['vort'] = asinh.transform(flat_vort)
            vert_grid.point_data['vert'] = asinh.transform(flat_vert)

            vert_grid_1.point_data['vert'] = asinh.transform(flat_vert_1)
            vert_grid_2.point_data['vert'] = asinh.transform(flat_vert_2)

            
            t_opacity = np.linspace(0, 1, 255)

            #Plug into desmos to see
            t_opacity_tf = 1.0-np.exp(-100000000.0*(t_opacity-0.5)**12)

            w_opacity = np.linspace(-1/2, 1/2, 255)
            w_opacity_tf = 1.0-np.exp(-3500.0*(w_opacity)**8)
            ω_opacity_tf = 1.0-np.exp(-3500.0*(w_opacity)**8)

            
            #inv_opacity = np.linspace(0, 1, 255)
            #inv_opacity_tf = np.exp(-100000.0*(opacity-0.5)**4)    
                    

            #grid.point_data['opacity'] = opacity
            #grid.point_data['inv_opacity'] = inv_opacity


            time_text = f"t = {f['scales/sim_time'][start+t]:.3f}"
            write_num = f['scales/write_number'][start+t]

            # Set up a 2x2 subplot
            plotter = pv.Plotter(shape=(3, 2), off_screen=True, border=False)

            # Top-left: diagonal view of T
            plotter.subplot(0, 0)
            actor1 = plotter.add_volume(
                t_grid, scalars='temp', opacity=t_opacity_tf, cmap='jet', clim=[0, 1]#, shade=False
            )
            actor1.mapper.interpolate_before_map = False
            plotter.camera_position = [(4, 4, 0.5), (xmid, ymid, zmid), (0, 0, 1)]
            plotter.add_text(time_text, position='upper_right', font_size=20, color='black')

            """            
            # Top-right: diagonal view of ω
            plotter.subplot(0, 1)
            actor1 = plotter.add_volume(
                vort_grid, scalars='vort', opacity=w_opacity_tf, cmap='jet', clim=[-0.8*mvort, mvort]#, shade=False
            )
            actor1.mapper.interpolate_before_map = False
            plotter.camera_position = [(4, 4, 0.5), (xmid, ymid, zmid), (0, 0, 1)]
            """

            # Top-right: diagonal view of w
            plotter.subplot(0, 1)
            actor1 = plotter.add_volume(
                vert_grid, scalars='vert', opacity=w_opacity_tf, cmap='RdBu_r', clim=[-mvert, mvert]#, shade=False
            )
            actor1.mapper.interpolate_before_map = False
            plotter.camera_position = [(4, 4, 0.5), (xmid, ymid, zmid), (0, 0, 1)]

            """
            # Top-right: rotating diagonal view
            angle = 2 * np.pi * (write_num % 150) / 150.0
            cam_x = xmid + rot_radius * np.cos(angle)
            cam_y = ymid + rot_radius * np.sin(angle)
            plotter.subplot(0, 1)
            actor2 = plotter.add_volume(
                wgrid, scalars='w', opacity=w_opacity_tf, cmap='jet', shade=False
            )
            actor2.mapper.interpolate_before_map = False
            plotter.camera_position = [(cam_x, cam_y, zmid), (xmid, ymid, zmid), (0, 0, 1)]
            """

            # middle-left: temp midplane
            plotter.subplot(1, 0)
            actor3 = plotter.add_mesh(
                t_grid_1, scalars='temp', cmap='RdBu_r'
            )
            actor3.mapper.interpolate_before_map = False
            plotter.enable_2d_style()
            plotter.enable_parallel_projection()
            plotter.view_xy()

            """
            # middle-right: temp BL - no clim!
            plotter.subplot(1, 1)
            actor5 = plotter.add_mesh(
                t_grid_2, scalars='temp', cmap='RdBu_r'
            )
            actor5.mapper.interpolate_before_map = False
            plotter.enable_2d_style()
            plotter.enable_parallel_projection()
            plotter.view_xy()
            """

            # middle-right: diagonal view of ω
            plotter.subplot(1, 1)
            actor1 = plotter.add_volume(
                vort_grid, scalars='vort', opacity=ω_opacity_tf, cmap='jet', clim=[-0.8*mvort, mvort]#, shade=False
            )
            actor1.mapper.interpolate_before_map = False
            plotter.camera_position = [(4, 4, 0.5), (xmid, ymid, zmid), (0, 0, 1)]
            

            # Bottom-left: w midplane
            plotter.subplot(2, 0)
            actor4 = plotter.add_mesh(
                vert_grid_1, scalars='vert', cmap='RdBu_r'#, clim=[-mvert, mvert]
            )
            actor4.mapper.interpolate_before_map = False
            plotter.enable_2d_style()
            plotter.enable_parallel_projection()
            plotter.view_xy()

            # Bottom-right: w bl - no clim!
            plotter.subplot(2, 1)
            actor6 = plotter.add_mesh(
                vert_grid_2, scalars='vert', cmap='RdBu_r'#, clim=[-mvert, mvert]
            )
            actor6.mapper.interpolate_before_map = False
            plotter.enable_2d_style()
            plotter.enable_parallel_projection()
            plotter.view_xy()
            plotter.reset_camera()

            
            # Save snapshot
            plotter.set_background('white')
            fname = os.path.join(output_dir, f"write_{write_num:06}.png")
            plotter.image_scale = 2
            plotter.screenshot(fname)
            plotter.close()

# Example usage
if __name__ == '__main__':
    
    from docopt import docopt
    from dedalus.tools import logging
    from dedalus.tools import post

    args = docopt(__doc__)

    basepath = Path(args["--basepath"])
    output_dir = basepath / "visualization"
    output_dir.mkdir(exist_ok=True)

    mvort = float(args["--max_vort"])
    mvert = float(args["--max_vert"])
    nu = float(args["--nu"])

    post.visit_writes(
        args["<files>"],
        main,
        output_dir=output_dir,
        mvort=mvort,
        mvert=mvert,
        nu=nu
    )

