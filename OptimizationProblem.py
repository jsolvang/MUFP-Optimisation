import numpy as np
from scipy import interpolate
import pickle
import sys
from FloaterParameters import FloaterParameters
from Environment import Environment
from Buoy import Buoy
from Mass import Mass
from Density import Density
from Area import Area
from GeneralisedCoordinateSystem import GeneralisedCoordinateSystem
from SystemMatrices import MatrixCalculation
from CheckInterpolation import InterpolateParameters
from WindForceSpectrum import WindForceSpectrum
from WaveSpectrumAnalysis import wavespectrumanalysis
from WindSpectrumAnalysis import windspectrumanalysis
from HydroDInterpolation import Interpolate2d
sys.path.append(r'C:\Users\Joar\Documents\1_Education\NTNU\OneDrive - NTNU\Thesis\Modelling\FD Model')


def OptimizationProblem(x_space, y_space, dia_column):
    # Defining static system matrices
    # x_space = x_space[0]
    plot = 0
    env = Environment()
    rho = Density()
    mufp = FloaterParameters(x_space, y_space, dia_column)
    csa = Area(mufp)
    buoy = Buoy(mufp, csa, rho)
    mass = Mass(mufp, csa, buoy, rho)
    coord = GeneralisedCoordinateSystem(mufp, csa, mass, rho, buoy, env)
    matrix = MatrixCalculation(coord, mass, mufp, rho, env, csa, buoy)
    print(matrix.mass[4,4])
    print(matrix.stiffness[4,4])
    # Adding artificial stiffness in DOF without any
    matrix.stiffness[0, 0] = 1e6
    matrix.stiffness[1, 1] = 1e6
    matrix.stiffness[5, 5] = 1e9

    # Defining return period for statistical evaluation
    Tr = 3 * 3600  # 3 hour return period

    # Finding the nearest HydroD simulations for surrogate model
    interp_inputs = np.array([[np.ceil(x_space / 10) * 10, np.ceil(y_space / 10) * 10, np.ceil(dia_column)],
                              [np.floor(x_space / 10) * 10, np.floor(y_space / 10) * 10, np.ceil(dia_column)],
                              [np.ceil(x_space / 10) * 10, np.ceil(y_space / 10) * 10, np.floor(dia_column)],
                              [np.floor(x_space / 10) * 10, np.floor(y_space / 10) * 10, np.floor(dia_column)]])

    # Interpolating surrogate models to get estimate
    design = [x_space, y_space, dia_column]
    results = []

    for i in np.linspace(0, len(interp_inputs[:, 1]) - 1, len(interp_inputs)).astype(int):
        file_loc = r'C:\Users\Joar\Documents\1_Education\NTNU\pickle_files'
        file_name = "\sim_x_%d_y_%d_D%d" % (interp_inputs[i, 0], interp_inputs[i, 1], interp_inputs[i, 2])
        file_path = file_loc + file_name
        infile = open(file_path, 'rb')
        results.append(pickle.load(infile))
        infile.close()

    sim = Interpolate2d(results, interp_inputs, design)

    # Defining frequencies of interest
    df = 1 / 3600
    df_rad = (2 * np.pi) * df
    f_rad = np.arange(0.01, 5, df_rad)

    # Spectral Moment Analysis
    tz_wave, sa_wave, n_mpm_wave, mpm_wave = wavespectrumanalysis(sim, matrix, plot, mufp, coord, df_rad, f_rad, Tr)
    sa_wave = np.nan_to_num(sa_wave)
    mpm_wave = np.nan_to_num(mpm_wave)


    # Wind Loading Spectrum Analysis
    # Pulling Wind Force Spectrum from SIMA simulation
    f_xx, P_wind, _, _ = WindForceSpectrum()
    f_xx = f_xx * 2 * np.pi             # converting to rad/s
    P_wind = P_wind / (2 * np.pi)       # converting to rad/s
    f = np.arange(0, f_xx[-1], df_rad)  # Defining new frequency range
    func_f = interpolate.interp1d(f_xx, P_wind)
    P_wind = func_f(f)

    # Applying force to surge and pitch
    s_wind = np.zeros(shape=(len(f), 8))
    s_wind[:, 0] = P_wind
    s_wind[:, 4] = P_wind
    tz_wind, sa_wind, n_mpm_wind, mpm_wind, _, _, _, _ = windspectrumanalysis(sim, matrix, plot, mufp, coord, df_rad, f, Tr, s_wind)
    sa_wind = np.nan_to_num(sa_wind)
    mpm_wind = np.nan_to_num(mpm_wind)

    mpm = mpm_wind + mpm_wave
    sa = sa_wind + sa_wave

    # Significant_Amplitude = Significant_Amplitude + Significant_Amplitude_wind
    # MPM = MPM + MPM_wind

    mpm_wave[:, 3:6] *= (180 / np.pi)
    sa_wave[:, 3:6] *= (180 / np.pi)
    mpm_wind[:, 3:6] *= (180 / np.pi)
    sa_wind[:, 3:6] *= (180 / np.pi)
    mpm[:, 3:6] *= (180 / np.pi)
    sa[:, 3:6] *= (180 / np.pi)

    return tz_wind, n_mpm_wind, sa_wind, mpm_wind, tz_wave, n_mpm_wave, sa_wave, mpm_wave, sa, mpm
