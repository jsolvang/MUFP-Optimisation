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
from GlobalCoordinateSystem import GlobalCoordinateSystem
from SystemMatrices import MatrixCalculation
from CheckInterpolation import InterpolateParameters
from WindForceSpectrum import WindForceSpectrum
from WaveSpectrumAnalysis import wavespectrumanalysis
from WindSpectrumAnalysis import windspectrumanalysis
from HydroDInterpolation import LinearInterpolate, Interpolate2d, Interpolate3d, Interpolate4d
from scipy import linalg
sys.path.append(r'C:\Users\Joar\Documents\1_Education\NTNU\OneDrive - NTNU\Thesis\Modelling\FD Model')


def OptimizationProblem(x_space, y_space, dia_column, draft):
    # Defining static system matrices
    plot = 0
    env = Environment()
    rho = Density()
    mufp = FloaterParameters(x_space, y_space, dia_column, draft)
    csa = Area(mufp)
    buoy = Buoy(mufp, csa, rho)
    mass = Mass(mufp, csa, buoy, rho)
    coord = GlobalCoordinateSystem(mufp, csa, mass, rho, buoy, env)
    matrix = MatrixCalculation(coord, mass, mufp, rho, env, csa, buoy)

    # Adding artificial stiffness in DOF without any hydrostatic stiffness
    matrix.stiffness[0, 0] = 1e6
    matrix.stiffness[1, 1] = 1e6
    matrix.stiffness[5, 5] = 1e9

    # Defining return period for statistical evaluation
    Tr = 3 * 3600  # 3 hour return period

    # Finding the nearest HydroD simulations for surrogate model
    interp_inputs = np.array(
        [[np.floor(x_space / 10) * 10, np.floor(y_space / 10) * 10, np.floor(dia_column), np.floor(draft / 2) * 2],
         [np.ceil(x_space / 10) * 10, np.floor(y_space / 10) * 10, np.floor(dia_column), np.floor(draft / 2) * 2],
         [np.floor(x_space / 10) * 10, np.ceil(y_space / 10) * 10, np.floor(dia_column), np.floor(draft / 2) * 2],
         [np.ceil(x_space / 10) * 10, np.ceil(y_space / 10) * 10, np.floor(dia_column), np.floor(draft / 2) * 2]])
                              #[np.ceil(x_space / 10) * 10, np.floor(y_space / 10) * 10, np.ceil(dia_column), np.ceil(draft / 2) * 2],
                              #[np.floor(x_space / 10) * 10, np.floor(y_space / 10) * 10, np.ceil(dia_column), np.ceil(draft / 2) * 2],
                              #[np.ceil(x_space / 10) * 10, np.ceil(y_space / 10) * 10, np.ceil(dia_column), np.ceil(draft / 2) * 2],
                              #[np.floor(x_space / 10) * 10, np.ceil(y_space / 10) * 10, np.ceil(dia_column), np.ceil(draft / 2) * 2]])

    # Interpolating surrogate models to get estimate
    design = [x_space, y_space, dia_column, draft]
    results = []

    for i in np.linspace(0, len(interp_inputs[:, 1]) - 1, len(interp_inputs)).astype(int):
        file_loc = r'F:\pickle_files'
        file_name = "\sim_x_%d_y_%d_D%d_dr%d" % (
            interp_inputs[i, 0], interp_inputs[i, 1], interp_inputs[i, 2], interp_inputs[i, 3])
        file_path = file_loc + file_name
        infile = open(file_path, 'rb')
        results.append(pickle.load(infile))
        infile.close()

    #sim = LinearInterpolate(results, interp_inputs, design, 0)
    sim = results[0] #Interpolate2d(results, interp_inputs, design, [0,1])

    # Adding 5% critical damping
    if matrix.stiffness[4, 4] > 0 and matrix.stiffness[3, 3] > 0:
        sim.DAMPING += (1/20)*2*np.sqrt(matrix.stiffness * (matrix.mass+sim.ADDEDMASS))
    else:
        print('Negative Stiffness - No Critical Damping')

    B_aero = np.zeros(shape=(6, 6))
    B_aero[0, 0] = 2*7.880219917522345e+04
    B_aero[0, 4] = 2*mufp.hub_height * 7.880219917522345e+04
    B_aero[4, 0] = 2*mufp.hub_height * 7.880219917522345e+04
    B_aero[4, 4] = 2*(mufp.hub_height ** 2) * 7.880219917522345e+04
    B_aero[5, 5] = 2 * (mufp.y_space/2)**2 * 7.880219917522345e+04

    sim.DAMPING += B_aero

    eig_freq, eig_vec = linalg.eig(matrix.stiffness,matrix.mass+sim.ADDEDMASS[5,:,:], right=True)
    eig_freq = np.sqrt(eig_freq)*np.sqrt(1-0.05)

    # Defining frequencies of interest
    # Breaking these up into three regions
    df_rad1 = (2 * np.pi) / 36000
    f_rad1 = np.arange(1e-7, 0.1, df_rad1)
    df_rad2 = (2 * np.pi) / 3600
    f_rad2 = np.arange(0.1, 1, df_rad2)
    df_rad3 = (2 * np.pi) / 3600
    f_rad3 = np.arange(1, 2.5, df_rad3)
    df_rad4 = (2 * np.pi) / 360
    f_rad4 = np.arange(2.5, 5, df_rad4)
    df_rad = np.concatenate(
        [np.ones(len(f_rad1)) * df_rad1, np.ones(len(f_rad2)) * df_rad2, np.ones(len(f_rad3)) * df_rad3, np.ones(len(f_rad4)) * df_rad4])
    f_rad = np.concatenate([f_rad1, f_rad2, f_rad3, f_rad4])

    # Spectral Moment Analysis
    tz_wave, sa_wave, n_mpm_wave, mpm_wave, rao_wave, _, _, _, wave_response = wavespectrumanalysis(sim, matrix, plot, mufp, coord, df_rad,
                                                                              f_rad, Tr)
    sa_wave = np.nan_to_num(sa_wave)
    mpm_wave = np.nan_to_num(mpm_wave)

    # Wind Loading Spectrum Analysis
    # Pulling Wind Force Spectrum from SIMA simulation
    f_xx, sx_wind, sy_wind, smz_wind = WindForceSpectrum()

    df_rad1 = (2 * np.pi) / 36000
    f_rad1 = np.arange(1e-7, 0.1, df_rad1)
    df_rad2 = (2 * np.pi) / 3600
    f_rad2 = np.arange(0.1, 1, df_rad2)
    df_rad3 = (2 * np.pi) / 3600
    f_rad3 = np.arange(1, 2.5, df_rad3)
    df_rad4 = (2 * np.pi) / 360
    f_rad4 = np.arange(2.5, f_xx[-1], df_rad4)
    df_rad = np.concatenate(
        [np.ones(len(f_rad1)) * df_rad1, np.ones(len(f_rad2)) * df_rad2, np.ones(len(f_rad3)) * df_rad3, np.ones(len(f_rad4)) * df_rad4])
    f = np.concatenate([f_rad1, f_rad2, f_rad3, f_rad4])

    func_f = interpolate.interp1d(f_xx, sx_wind)
    sx_wind = func_f(f)

    func_f = interpolate.interp1d(f_xx, sy_wind)
    sy_wind = func_f(f)

    func_f = interpolate.interp1d(f_xx, smz_wind)
    smz_wind = func_f(f)

    # Applying full force to surge and pitch & misalignment force to sway and roll
    s_wind = np.zeros(shape=(len(f), 8))
    s_wind[:, 0] = sx_wind
    s_wind[:, 1] = sy_wind
    s_wind[:, 3] = sy_wind
    s_wind[:, 4] = sx_wind
    s_wind[:, 5] = smz_wind

    tz_wind, sa_wind, n_mpm_wind, mpm_wind, wind_response, rao_wind, Y, H = windspectrumanalysis(sim, matrix,
                                                                                                           plot, mufp,
                                                                                                           coord,
                                                                                                           df_rad, f,
                                                                                                           Tr, s_wind)
    sa_wind = np.nan_to_num(sa_wind)
    mpm_wind = np.nan_to_num(mpm_wind)

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

    mpm_wave[:, 3:6] *= (180 / np.pi)
    sa_wave[:, 3:6] *= (180 / np.pi)
    mpm_wind[:, 3:6] *= (180 / np.pi)
    sa_wind[:, 3:6] *= (180 / np.pi)
    mpm[:, :, 3:6] *= (180 / np.pi)
    sa[:, :, 3:6] *= (180 / np.pi)

    return sa, mpm, rao_wave, rao_wind, matrix.stiffness, eig_freq, eig_vec, mass.steel_mass, wind_response, tz_wind, wave_response


def constraint3(x0):
    dia = x0[2]/10
    draft = x0[3]/10
    mufp = FloaterParameters(x0[0],x0[1],dia,draft)
    env = Environment()
    rho = Density()
    csa = Area(mufp)
    buoy = Buoy(mufp, csa, rho)
    mass = Mass(mufp, csa, buoy, rho)
    coord = GlobalCoordinateSystem(mufp, csa, mass, rho, buoy, env)
    matrix = MatrixCalculation(coord, mass, mufp, rho, env, csa, buoy)
    heave = (mass.total - buoy.total)/matrix.stiffness[2,2]
    return -heave