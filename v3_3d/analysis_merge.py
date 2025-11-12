from pathlib import Path
import h5py as h5
import numpy as np

fp = Path("/home/ollie/links/scratch/rbc_scripts/v3_3d/analysis")
fs = fp.glob('*.h5')
f0 = next(fs)

with h5.File(f0, 'r') as f:
    Ra = f['tasks']['Ra'][-1]
    Pr = f['tasks']['Pr'][-1]
    z_an = f['tasks']['z_an'][:]

    full_time = f['scales']['sim_time'][:]
    avg_K = f['tasks']['avg_K'][:]
    avg_wT = f['tasks']['avg_wT'][:]
    avg_vorticity_sq = f['tasks']['avg_vorticity_sq'][:]
    avg_grad_T_sq = f['tasks']['avg_grad_T_sq'][:]
    avg_T = f['tasks']['avg_T'][:]
    avg_u_sq = f['tasks']['avg_u_sq'][:]
    avg_v_sq = f['tasks']['avg_v_sq'][:]
    avg_w_sq = f['tasks']['avg_w_sq'][:]
    print(full_time.shape)

for fi in fs:
    with h5.File(fi, 'r') as f:
        full_time = np.append(full_time, f['scales']['sim_time'][:])
        avg_K = np.append(avg_K, f['tasks']['avg_K'][:])
        avg_wT = np.append(avg_wT, f['tasks']['avg_wT'][:])
        avg_vorticity_sq = np.append(avg_vorticity_sq, f['tasks']['avg_vorticity_sq'][:])
        avg_grad_T_sq = np.append(avg_grad_T_sq, f['tasks']['avg_grad_T_sq'][:])
        avg_T = np.append(avg_T, f['tasks']['avg_T'][:])
        avg_u_sq = np.append(avg_u_sq, f['tasks']['avg_u_sq'][:])
        avg_v_sq = np.append(avg_v_sq, f['tasks']['avg_v_sq'][:])
        avg_w_sq = np.append(avg_w_sq, f['tasks']['avg_w_sq'][:])
        print(full_time.shape)

