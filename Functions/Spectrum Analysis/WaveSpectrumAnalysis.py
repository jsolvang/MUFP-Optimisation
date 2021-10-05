from jonswap import jonswap
from RAO_Calculation import calulate_RAOs
from SpectrumStatistics import spectrum_response_statistics, spectrum_response
from InterpolateHydroDResults import InterpolateResults
import numpy as np

def wavespectrumanalysis(sim, matrix, plot, mufp, coord, df, f, Tr):
    # Defining Jonswap variables
    Hs = 4.1
    Tp = 10.2
    gammaJS = 3.3

    # Interpolating HydroD results to smaller intervals
    sim_interp = InterpolateResults(sim, f)

    # Deriving system RAOs using hand calculated static matrices and interpolated hydrodynamic coefficients
    RAO, Y, H, Y_inv = calulate_RAOs(sim_interp.wave_disc, matrix.mass, sim_interp.ADDEDMASS, sim_interp.DAMPING,
                                 matrix.stiffness,
                                 sim_interp.WAVEEX, plot, mufp, coord)

    # Wave Loading Frequency Domain Analysis
    s, _ = jonswap(Hs, Tp, df, gammaJS, f, plot)
    s_mat = np.ones(shape=(len(s), 8))

    for ii in np.arange(0, 8, 1):
        s_mat[:, ii] = s

    # Creating empty variables
    Tz = np.zeros(shape=(2, 3, 8))
    Significant_Amplitude = np.zeros(shape=(2, 3, 8))
    N_mpm = np.zeros(shape=(2, 3, 8))
    MPM = np.zeros(shape=(2, 3, 8))

    # Calculating response spectrum from RAO and excitation spectrum
    response = np.zeros(shape=(2, len(RAO[0, :, 0, 0]), len(s), 8))
    response[0, :, :, :] = spectrum_response(RAO[0, :, :, :], s_mat)
    response[1, :, :, :] = spectrum_response(RAO[1, :, :, :], s_mat)

    # Deriving system response characteristics from response spectrum
    Tz[0, :, :], Significant_Amplitude[0, :, :], N_mpm[0, :, :], MPM[0, :, :], _, _, _ = spectrum_response_statistics(
        Tr, response[0, :, :, :], f, df)
    Tz[1, :, :], Significant_Amplitude[1, :, :], N_mpm[1, :, :], MPM[1, :, :], _, _, _ = spectrum_response_statistics(
        Tr, response[1, :, :, :], f, df)

    return Significant_Amplitude, MPM
