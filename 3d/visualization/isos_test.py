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


import h5py
import numpy as np
import pyvista as pv
import os
from pathlib import Path

def main(h5_file, start, count, output_dir, iso_temp=None, slice_axis='z'):
    # Load temp data and grid
    pv.start_xvfb()
    pv.global_theme.allow_empty_mesh = True
    with h5py.File(h5_file, 'r') as f:
        T = f['tasks']['temp']
        w = f['tasks']['w']
        nt, nx, ny, nz = T.shape
        
        domain_sizes = (2.0, 2.0, 1.0)
        dx = domain_sizes[0] / (nx - 1)
        dy = domain_sizes[1] / (ny - 1)
        dz = domain_sizes[2] / (nz - 1)
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
        for t in range(start, start+count):
            temp = T[t]  # (nx, ny, nz)
            vert_velocity = w[t]

            flat_temp = temp.flatten(order='F')
            flat_velocity = vert_velocity.flatten(order='F')

            # Create uniform grid
            grid = pv.ImageData(
                dimensions=(nx, ny, nz),
                spacing=(dx, dy, dz),
                origin=origin
            )

            grid.point_data['temp'] = flat_temp
            grid.point_data['w'] = flat_velocity
            
            opacity = np.linspace(0, 1, 255)
            opacity_tf = 1.0-np.exp(-100000000.0*(opacity-0.5)**12)

            velocity_opacity = np.linspace(0, 1, 255)
            velocity_opacity_tf = 1.0-np.exp(-10000000.0*(opacity-0.5)**6)

            #inv_opacity = np.linspace(0, 1, 255)
            #inv_opacity_tf = np.exp(-100000.0*(opacity-0.5)**4)    
                    

            #grid.point_data['opacity'] = opacity
            #grid.point_data['inv_opacity'] = inv_opacity


            time_text = f"t = {f['scales/sim_time'][t]:.3f}"
            write_num = f['scales/write_number'][t]

            # Set up a 2x2 subplot
            plotter = pv.Plotter(shape=(2, 2), off_screen=True, border=False)

            # Top-left: diagonal view
            plotter.subplot(0, 0)
            actor1 = plotter.add_volume(
                grid, scalars='temp', opacity=opacity_tf, cmap='jet', clim=[0, 1], shade=False
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
                grid, scalars='temp', opacity=opacity_tf, cmap='jet', clim=[0, 1], shade=False
            )
            actor2.mapper.interpolate_before_map = False
            plotter.camera_position = [(cam_x, cam_y, zmid), (xmid, ymid, zmid), (0, 0, 1)]




            # Bottom-left: diagonal view
            plotter.subplot(1, 0)
            actor3 = plotter.add_volume(
                grid, scalars='w', cmap='jet', opacity=opacity_tf, clim=[-0.6,0.6],shade=False
            )
            actor3.mapper.interpolate_before_map = False
            plotter.camera_position = [(4, 4, 0.5), (xmid, ymid, zmid), (0, 0, 1)]


            # Bottom-right: rotating diagonal view
            plotter.subplot(1, 1)
            actor4 = plotter.add_volume(
                grid, scalars='w', cmap='jet', opacity=opacity_tf, clim=[-0.6,0.6],shade=True
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

