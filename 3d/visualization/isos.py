"""
Script to perform analysis tasks on data obtained from running a 3D
Rayleigh-Benard convection simulation. Accepts an analysis file as a command
line argument and takes time averages from the specified time. All output is
written to the folders 'outputs' and 'preliminary_outputs'. Outputs plots of 
the Nusselt number, various profiles and information about the simulation to a 
text document.

Usage:
    isos.py <files>... [--basepath=<dir>]
    isos.py <files>...

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
    with h5py.File(h5_file, 'r') as f:
        T = f['tasks']['temp']
        nt, nx, ny, nz = T.shape
        coords = {
            'x': np.linspace(0, 2, nx),
            'y': np.linspace(0, 2, ny),
            'z': np.linspace(-1/2, 1/2, nz)
        }
        x, y, z = coords['x'], coords['y'], coords['z']

        # Ensure output directory exists
        os.makedirs(output_dir, exist_ok=True)

        # Loop through each time step
        for t in range(start, start+count):
            print(t)
            temp = T[t]  # (nx, ny, nz)
            grid = pv.StructuredGrid(*np.meshgrid(x, y, z, indexing='ij'))
            grid['temp'] = temp.flatten(order='F')

            p = pv.Plotter(off_screen=True)
            p.add_axes()

            time_text = f"t = {f['scales/sim_time'][t]:.3f}"
            p.add_text(time_text, position='lower_left', font_size=14, color='black')
            
            #  Add planar slices
            axes = ['x', 'y', 'z', 'z']
            origins = [(1e-16,0,0), (0,1e-16,0), (0,0,(-1/2)+1e-16), (0,0,1/2-1e-16)]
            for ax, og in zip(axes,origins,ks):
                slice_mesh = grid.slice(normal=ax, origin=og)
                p.add_mesh(slice_mesh, scalars='temp', cmap='RdBu_r',
                        show_scalar_bar=False, opacity=1.0)

            
            iso = grid.contour(isosurfaces=[0.4, 0.6], scalars='temp')
            p.add_mesh(iso, cmap='RdBu_r', clim=[0.2, 0.8] opacity=0.75)

            # 4 Finalize view and save image
            p.set_background('white')
            p.camera_position = [(4, 4, 0), (0, 0, 0), (0, 0, 1)]
            fname = os.path.join(output_dir, f"write_{t:06d}.png")
            p.show(screenshot=fname)
            p.close()

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

