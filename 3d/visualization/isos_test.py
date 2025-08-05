"""
Script to perform analysis tasks on data obtained from running a 3D
Rayleigh-Benard convection simulation. Accepts an analysis file as a command
line argument and takes time averages from the specified time. All output is
written to the folders 'outputs' and 'preliminary_outputs'. Outputs plots of 
the Nusselt number, various profiles and information about the simulation to a 
text document.

Usage:
    isos_test.py <files>... [--basepath=<dir>]
    isos_test.py <files>...

Options:
    --basepath=<dir>  Path to parent folder for output [default: ./analysis]
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

def main(h5_file, start, count, output_dir, iso_temp=None, slice_axis='z'):
    # Load temp data and grid
    pv.start_xvfb()
    pv.global_theme.allow_empty_mesh = True
    with h5py.File(h5_file, 'r') as f:
        T = np.array(f['tasks']['temp'][start:start+count], dtype=np.float32)
        w = np.array(f['tasks']['w'][start:start+count], dtype=np.float32)
        x_vort = np.array(f['tasks']['x_vort'][start:start+count], dtype=np.float32)
        y_vort = np.array(f['tasks']['y_vort'][start:start+count], dtype=np.float32)

        nt, nx, ny, nz = T.shape
        _, nwx, nwy, nwz = w.shape
        _, nxvx, nxvy, nxvz = x_vort.shape
        _, nyvx, nyvy, nyvz = y_vort.shape

        domain_sizes = (2.0, 2.0, 1.0)
        dx = domain_sizes[0] / (nx - 1)
        dy = domain_sizes[1] / (ny - 1)
        dz = domain_sizes[2] / (nz - 1)

        wdx = domain_sizes[0] / (nwx - 1)
        wdy = domain_sizes[1] / (nwy - 1)
        wdz = domain_sizes[2] / (nwz - 1)

        xvdx = domain_sizes[0] / (nxvx - 1)
        xvdy = domain_sizes[1] / (nxvy - 1)
        xvdz = domain_sizes[2] / (nxvz - 1)

        yvdx = domain_sizes[0] / (nyvx - 1)
        yvdy = domain_sizes[1] / (nyvy - 1)
        yvdz = domain_sizes[2] / (nyvz - 1)

        origin = (0.0, 0.0, 0.0)
        # Center of domain
        xmid = domain_sizes[0] / 2.0
        ymid = domain_sizes[1] / 2.0
        zmid = domain_sizes[2] / 2.0

        # Rotating view radius
        rot_radius = 4.5

        # Ensure output directory exists
        os.makedirs(output_dir, exist_ok=True)

        # Loop through each time step
        for t in range(0, count):
            temp = T[t]  # (nx, ny, nz)
            vert_velocity = w[t]
            xv = x_vort[t]
            yv = y_vort[t]

            flat_temp = temp.flatten(order='F')
            flat_velocity = vert_velocity.flatten(order='F')
            flat_xv = xv.flatten(order='F')
            flat_yv = yv.flatten(order='F')

            # Create uniform grid
            tgrid = pv.ImageData(
                dimensions=(nx, ny, nz),
                spacing=(dx, dy, dz),
                origin=origin
            )

            wgrid = pv.ImageData(
                dimensions=(nwx, nwy, nwz),
                spacing=(wdx, wdy, wdz),
                origin=origin
            )

            xv_grid = pv.ImageData(
                dimensions=(nxvx, nxvy, nxvz),
                spacing=(xvdx, xvdy, xvdz),
                origin=origin
            )

            yv_grid = pv.ImageData(
                dimensions=(nyvx, nyvy, nyvz),
                spacing=(yvdx, yvdy, yvdz),
                origin=origin
            )

            tgrid.point_data['temp'] = flat_temp
            wgrid.point_data['w'] = flat_velocity

            xv_grid.point_data['xv'] = asinh.transform(flat_xv)
            yv_grid.point_data['yv'] = asinh.transform(flat_yv)
            
            t_opacity = np.linspace(0, 1, 255)
            t_opacity_tf = 1.0-np.exp(-100000000.0*(t_opacity-0.5)**12)

            w_opacity = np.linspace(-1/2, 1/2, 255)
            w_opacity_tf = 1.0-np.exp(-10000.0*(w_opacity)**6)

            
            #inv_opacity = np.linspace(0, 1, 255)
            #inv_opacity_tf = np.exp(-100000.0*(opacity-0.5)**4)    
                    

            #grid.point_data['opacity'] = opacity
            #grid.point_data['inv_opacity'] = inv_opacity


            time_text = f"t = {f['scales/sim_time'][start+t]:.3f}"
            write_num = f['scales/write_number'][start+t]

            # Set up a 2x2 subplot
            plotter = pv.Plotter(shape=(2, 2), off_screen=True, border=False)

            # Top-left: diagonal view - CHANGED REMOVE INTERPOLATE TO BE IN ADD_MESH
            plotter.subplot(0, 0)
            actor1 = plotter.add_volume(
                tgrid, scalars='temp', opacity=t_opacity_tf, cmap='jet', clim=[0, 1], shade=False
            )
            actor1.mapper.interpolate_before_map = False
            plotter.camera_position = [(4, 4, 0.5), (xmid, ymid, zmid), (0, 0, 1)]
            plotter.add_text(time_text, position='upper_left', font_size=14, color='black')

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




            # Bottom-left: xvort
            plotter.subplot(1, 0)
            actor3 = plotter.add_volume(
                xv_grid, scalars='xv', cmap='RdBu_r', opacity=w_opacity_tf, shade=False
            )
            actor3.mapper.interpolate_before_map = False
            plotter.camera_position = [(4, 4, 0.5), (xmid, ymid, zmid), (0, 0, 1)]


            # Bottom-right: yvort
            plotter.subplot(1, 1)
            actor4 = plotter.add_volume(
                yv_grid, scalars='yv', cmap='RdBu_r', opacity=w_opacity_tf,shade=False
            )
            actor4.mapper.interpolate_before_map = False
            plotter.camera_position = [(cam_x, cam_y, zmid), (xmid, ymid, zmid), (0, 0, 1)]

            

            # Save snapshot
            plotter.set_background('white')
            fname = os.path.join(output_dir, f"write_{write_num:06}.png")
            plotter.show(screenshot=fname)
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

    post.visit_writes(
        args["<files>"],
        main,
        output_dir=output_dir
    )

