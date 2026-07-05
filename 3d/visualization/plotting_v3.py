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
    --max_vort=<max_vort>  Maximum vorticity [default: 1000]
    --max_vert=<max_vert>  Maximum vertical velocity [default: 5000]
    --nu=<nu> Nusselt from this sim for boundary layer [default: 5]
"""

#TODO: choose consistent scales and then reduce number of grids.
# TODO: change the x,y,z to be pulled from the scales
# TODO: change the asinh transform to be a function on the scalars rather than change the data?

import h5py
import numpy as np
import pyvista as pv
import os
from pathlib import Path
from matplotlib.scale import AsinhTransform

# For vorticity scaling
asinh = AsinhTransform(linear_width=1.5)

def main(h5_file, start, count, output_dir, mvort, mvert, nu):
    mvort = asinh.transform(mvort)
    mvert = asinh.transform(mvert)
    bl = (4.0*nu)
    # Load temp data and grid
    pv.start_xvfb()
    pv.global_theme.allow_empty_mesh = True
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    with h5py.File(h5_file, 'r') as f:
        
        _, nTx, nTy, nTz = f['tasks']['temperature'].shape
        _, nωx, nωy, nωz = f['tasks']['vorticity'].shape
        try:
            _, nwx, nwy, nwz = f['tasks']['w'].shape
        except:
            _, nwx, nwy, nwz = f['tasks']['u'][0].shape
        domain_sizes = (2.0, 2.0, 1.0)
        # Center of domain
        xmid = domain_sizes[0] / 2.0
        ymid = domain_sizes[1] / 2.0
        zmid = domain_sizes[2] / 2.0

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
        
        # Loop through each time step
        for ti in range(0, count):
            time = f['scales']['sim_time'][ti]

            temp = np.array(f['tasks']['temperature'][ti], dtype=np.float64)
            vort = np.array(f['tasks']['vorticity'][ti], dtype=np.float64)
            try:
                vert = np.array(f['tasks']['w'][ti], dtype=np.float64)
            except:
                vert = np.array(f['tasks']['u'][ti, 2], dtype=np.float64)
            #T_1 = np.array(f['tasks']['temperature'][ti, :, :, int(nTz/2)], dtype=np.float64)
            #T_2 = np.array(f['tasks']['temperature'][ti, :, :, int(nTz/(bl))], dtype=np.float64)

            

            #w_1 = np.array(f['tasks']['w'][ti, :, :, int(nwz/2)], dtype=np.float64)
            #w_2 = np.array(f['tasks']['w'][ti, :, :, int(nwz/(bl))], dtype=np.float64)


            #temp_1 = T_1[t]
            #temp_2 = T_2[t]

            #vert_1 = w_1[t]
            #vert_2 = w_2[t]
            
            flat_temp = temp.flatten(order='F')
            flat_vort = vort.flatten(order='F')
            flat_vert = vert.flatten(order='F')
            
            #flat_temp_1 = temp[:,:,int(nTz/2)].flatten(order='F')
            #flat_temp_2 = temp[:,:,int(nTz/bl)].flatten(order='F')

            #flat_vert_1 = vert[:,:,int(nwz/2)].flatten(order='F')
            #flat_vert_2 = vert[:,:,int(nwz/bl)].flatten(order='F')

            # Create uniform grid
            t_grid = pv.ImageData(
                dimensions=(nTx, nTy, nTz),
                spacing=(Tdx, Tdy, Tdz),
                origin=origin
            )
            """
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
            """

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
            """
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
            """
            t_grid.point_data['temp'] = flat_temp
            #t_grid_1.point_data['temp'] = flat_temp_1
            #t_grid_2.point_data['temp'] = flat_temp_2
            
            vort_grid.point_data['vort'] = asinh.transform(flat_vort)
            vert_grid.point_data['vert'] = asinh.transform(flat_vert)
            #vert_grid_1.point_data['vert'] = asinh.transform(flat_vert_1)
            #vert_grid_2.point_data['vert'] = asinh.transform(flat_vert_2)
                

            
            t_opacity = np.linspace(0, 1, 255)

            #Plug into desmos to see
            t_opacity_tf = 1.0-np.exp(-100000000.0*(t_opacity-0.5)**12)

            w_opacity = np.linspace(-1/2, 1/2, 255)
            w_opacity_tf =  np.exp(-3500.0*(w_opacity)**8)
            min_vort = np.min(asinh.transform(flat_vort))
            max_vort = np.max(asinh.transform(flat_vort))
            lin_color = np.linspace(min_vort, max_vort, 256)
            zero_start = np.searchsorted(lin_color, 0)

            ω_opacity_tf = 1.0-np.exp(-300000000000.0*(w_opacity)**22)
            # ω_opacity_tf = 1.0-np.exp(-300000000000.0*(w_opacity)**26)
            if zero_start < 125:
                if zero_start > 0:
                    ω_opacity_tf[zero_start-1:125] = np.zeros(125-zero_start+1)
                else:
                    ω_opacity_tf[zero_start:125] = np.zeros(125-zero_start)
            
            #inv_opacity = np.linspace(0, 1, 255)
            #inv_opacity_tf = np.exp(-100000.0*(opacity-0.5)**4)    
                    

            #grid.point_data['opacity'] = opacity
            #grid.point_data['inv_opacity'] = inv_opacity


            time_text = f"t = {f['scales/sim_time'][start+ti]:.3f}"
            write_num = f['scales/write_number'][start+ti]

            # Set up a 2x2 subplot
            plotter = pv.Plotter(shape=(1, 2), off_screen=True, border=False)
            plotter.window_size = [840,300]
            plotter.add_text(time_text, position='upper_edge', font_size=36, color='k')
            # Top-left: diagonal view of T
            plotter.subplot(0, 0)
            actor1 = plotter.add_volume(
                t_grid, scalars='temp', opacity=t_opacity_tf, cmap='jet', clim=[0, 1], show_scalar_bar=False#, shade=False
            )
            plotter.add_scalar_bar(
                title='Temperature',
                width=0.75,
                height=0.1,
                position_x=0.125,
                position_y=0.04,

            )
            actor1.mapper.interpolate_before_map = False
            plotter.camera_position = [(4, 4, 1/2), (xmid, ymid, zmid-0.15), (0, 0, 1)]
            #plotter.add_text(time_text, position='upper_edge', font_size=36, color='k')


            
            # Top-middle: diagonal view of w
            plotter.subplot(0, 1)
            actor1 = plotter.add_volume(
                vort_grid, scalars='vort', opacity=ω_opacity_tf, cmap='jet', show_scalar_bar=False#, clim=[-mvert, mvert]#, shade=False
            )
            plotter.add_scalar_bar(
                title='Magnitude of Vorticity',
                width=0.75,
                height=0.1,
                position_x=0.125,
                position_y=0.04,

            )
            actor1.mapper.interpolate_before_map = False
            plotter.camera_position = [(4, 4, 1/2), (xmid, ymid, zmid-0.15), (0, 0, 1)]
            #plotter.add_text(time_text, position='upper_left', font_size=36, color='red')
            
            
            # Save snapshot
            plotter.set_background('white')
            fname = os.path.join(output_dir, f"write_{write_num:06}.png")
            plotter.image_scale = 4
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

