import numpy as np
import os
from scipy import interpolate
import pickle
import sys
dirname = os.path.dirname(__file__)
sys.path.append(os.path.join(dirname, 'Functions/Data Scraping'))
sys.path.append(os.path.join(dirname, 'Functions/Hydrostatic Analysis'))
sys.path.append(os.path.join(dirname, 'Functions/Hydrodynamic Analysis'))
sys.path.append(os.path.join(dirname, 'Functions/Spectrum Analysis'))
from FloaterParameters import FloaterParameters
from Environment import Environment
from Buoy import Buoy
from Mass import Mass
from Density import Density
from Area import Area
from GlobalCoordinateSystem import GlobalCoordinateSystem
from SystemMatrices import MatrixCalculation
from CheckInterpolation import InterpolateParameters
from WindForceSpectrum import WindForceSpectrum
from WaveSpectrumAnalysis import wavespectrumanalysis
from WindSpectrumAnalysis import windspectrumanalysis
from HydroDInterpolation import LinearInterpolate, Interpolate2d, Interpolate4d


def ObjectiveFunction(x_space, y_space, dia_column, draft):
    # Objective function calculated the Significant Amplitude (SA) response, Most Probable Maxima (MPM) response,
    # eigenfrequencies and mass of a Multi-Unit Floating Platform

    # Defining system GCS and hydrostatic system matrices
    plot = 0
    env = Environment()
    rho = Density()
    mufp = FloaterParameters(x_space, y_space, dia_column, draft)
    csa = Area(mufp)
    buoy = Buoy(mufp, csa, rho)
    mass = Mass(mufp, csa, buoy, rho)
    coord = GlobalCoordinateSystem(mufp, csa, mass, rho, buoy, env)
    matrix = MatrixCalculation(coord, mass, mufp, rho, env, csa, buoy)

    # Adding artificial stiffness to degrees of freedom without hydrostatic stiffness
    matrix.stiffness[0, 0] = 1e6
    matrix.stiffness[1, 1] = 1e6
    matrix.stiffness[5, 5] = 1e9

    # Defining desired statistical return period
    Tr = 3 * 3600  # 3 hour return period (s)

    # Defining iteration design
    iter_design = [x_space, y_space, dia_column, draft]

    # Finding the nearest HydroD simulations for surrogate model
    interp_inputs = np.array(
        [[np.floor(x_space / 10) * 10, np.floor(y_space / 10) * 10, np.floor(dia_column), np.floor(draft / 2) * 2],
         [np.ceil(x_space / 10) * 10, np.floor(y_space / 10) * 10, np.floor(dia_column), np.floor(draft / 2) * 2],
         [np.floor(x_space / 10) * 10, np.ceil(y_space / 10) * 10, np.floor(dia_column), np.floor(draft / 2) * 2],
         [np.ceil(x_space / 10) * 10, np.ceil(y_space / 10) * 10, np.floor(dia_column), np.floor(draft / 2) * 2],
         [np.floor(x_space / 10) * 10, np.floor(y_space / 10) * 10, np.ceil(dia_column), np.floor(draft / 2) * 2],
         [np.ceil(x_space / 10) * 10, np.floor(y_space / 10) * 10, np.ceil(dia_column), np.floor(draft / 2) * 2],
         [np.floor(x_space / 10) * 10, np.ceil(y_space / 10) * 10, np.ceil(dia_column), np.floor(draft / 2) * 2],
         [np.ceil(x_space / 10) * 10, np.ceil(y_space / 10) * 10, np.ceil(dia_column), np.floor(draft / 2) * 2],
         [np.floor(x_space / 10) * 10, np.floor(y_space / 10) * 10, np.floor(dia_column), np.ceil(draft / 2) * 2],
         [np.ceil(x_space / 10) * 10, np.floor(y_space / 10) * 10, np.floor(dia_column), np.ceil(draft / 2) * 2],
         [np.floor(x_space / 10) * 10, np.ceil(y_space / 10) * 10, np.floor(dia_column), np.ceil(draft / 2) * 2],
         [np.ceil(x_space / 10) * 10, np.ceil(y_space / 10) * 10, np.floor(dia_column), np.ceil(draft / 2) * 2],
         [np.floor(x_space / 10) * 10, np.floor(y_space / 10) * 10, np.ceil(dia_column), np.ceil(draft / 2) * 2],
         [np.ceil(x_space / 10) * 10, np.floor(y_space / 10) * 10, np.ceil(dia_column), np.ceil(draft / 2) * 2],
         [np.floor(x_space / 10) * 10, np.ceil(y_space / 10) * 10, np.ceil(dia_column), np.ceil(draft / 2) * 2],
         [np.ceil(x_space / 10) * 10, np.ceil(y_space / 10) * 10, np.ceil(dia_column), np.ceil(draft / 2) * 2]])

    # Loading simulations from binary database
    results = []
    for i in np.linspace(0, len(interp_inputs[:, 1]) - 1, len(interp_inputs)).astype(int):
        file_loc = r'E:\pickle_files'
        file_name = "\sim_x_%d_y_%d_D%d_dr%d" % (
            interp_inputs[i, 0], interp_inputs[i, 1], interp_inputs[i, 2], interp_inputs[i, 3])
        file_path = file_loc + file_name
        infile = open(file_path, 'rb')
        results.append(pickle.load(infile))
        infile.close()
    sim = Interpolate4d(results, interp_inputs, iter_design)

    # Adding 5% critical damping to account for viscous damping
    if matrix.stiffness[4, 4] > 0:
        sim.DAMPING += (1 / 20) * 2 * np.sqrt(matrix.stiffness * (matrix.mass + sim.ADDEDMASS))

    # Adding linearised aerodynamic damping
    B_aero = np.zeros(shape=(6, 6))
    B_aero[0, 0] = 2 * 7.880e+04
    B_aero[0, 4] = 2 * mufp.hub_height * 7.880e+04
    B_aero[4, 0] = 2 * mufp.hub_height * 7.880e+04
    B_aero[4, 4] = 2 * (mufp.hub_height ** 2) * 7.880e+04
    B_aero[5, 5] = 2 * (mufp.y_space / 2) ** 2 * 7.880e+04
    sim.DAMPING += B_aero

    # Defining frequencies evaluations (radians)
    df1 = (2 * np.pi) / 36000
    df2 = (2 * np.pi) / 3600
    df3 = (2 * np.pi) / 360

    f1 = np.arange(1e-7, 2, df1)
    f2 = np.arange(2, 2.5, df2)
    f3 = np.arange(2.5, 5, df3)

    df = np.concatenate([np.ones(len(f1)) * df1, np.ones(len(f2)) * df2, np.ones(len(f3)) * df3])
    f = np.concatenate([f1, f2, f3])

    # Spectral Moment Analysis - Wave Loading
    sa_wave, mpm_wave = wavespectrumanalysis(sim, matrix, plot, mufp, coord, df, f, Tr)
    sa_wave = np.nan_to_num(sa_wave)
    mpm_wave = np.nan_to_num(mpm_wave)

    # Spectral Moment Analysis - Wave Loading
    sa_wind, mpm_wind = windspectrumanalysis(sim, matrix, plot, mufp, coord, df, f, Tr)
    sa_wind = np.nan_to_num(sa_wind)
    mpm_wind = np.nan_to_num(mpm_wind)

    # Summing wind and wave spectral responses
    mpm = np.zeros(mpm_wave.shape)
    sa = np.zeros(mpm_wave.shape)

    mpm[0, :, [0, 2, 4, 5]] = mpm_wind[:, [0, 2, 4, 5]].transpose()
    mpm[0, :, :] += mpm_wave[0, :, :]

    sa[0, :, [0, 2, 4, 5]] = sa_wind[:, [0, 2, 4, 5]].transpose()
    sa[0, :, :] += sa_wave[0, :, :]

    mpm[1, :, [1, 3]] = mpm_wind[:, [1, 3]].transpose()
    mpm[1, :, :] += mpm_wave[1, :, :]

    sa[1, :, [1, 3]] = sa_wind[:, [1, 3]].transpose()
    sa[1, :, :] += mpm_wave[1, :, :]

    # Converting from radians to degrees
    mpm_wave[:, 3:6] *= (180 / np.pi)
    sa_wave[:, 3:6] *= (180 / np.pi)
    mpm_wind[:, 3:6] *= (180 / np.pi)
    sa_wind[:, 3:6] *= (180 / np.pi)
    mpm[:, :, 3:6] *= (180 / np.pi)
    sa[:, :, 3:6] *= (180 / np.pi)

    return sa, mpm, matrix.stiffness, mass.steel_mass
