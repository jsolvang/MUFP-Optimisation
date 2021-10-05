from scipy import interpolate
import numpy as np


# Interpolation function
def interpolate_panel_pressure(sim, column_pressure, f_rad):
    panel_pressure_interp = np.zeros(shape=(2, len(f_rad), len(column_pressure[0, 0, :, 0]), 21))

    for ii, _ in enumerate(column_pressure[:, 0, 0, 0]):
        for jj, _ in enumerate(column_pressure[0, 0, :, 0]):
            for kk, _ in enumerate(column_pressure[0, 0, 0, :]):
                func_f = interpolate.interp1d(sim.wave_disc[:, 4], column_pressure[ii, :, jj, kk])
                panel_pressure_interp[ii, :, jj, kk] = func_f(f_rad)
    return panel_pressure_interp
