from jonswap import jonswap
from RAO_Calculation import calulate_RAOs
from SpectrumStatistics import spectrum_response_statistics, spectrum_response
from InterpolateHydroDResults import InterpolateResults

import numpy as np


def wavespectrumanalysis(sim, matrix, plot, mufp, coord, df, f, Tr):
    # Defining Jonswap variables
    Hs = 6
    Tp = 10
    gammaJS = 3.3

    # Interpolating HydroD results to smaller intervals
    sim_interp = InterpolateResults(sim.wave_disc, sim.WAVEEX, sim.MOTIONS, sim.ADDEDMASS, sim.DAMPING, f)

    # Deriving system RAOs using hand calculated static matrices and interpolated hydrodynamic coefficients
    RAO, _, _, _ = calulate_RAOs(sim_interp.wave_disc, matrix.mass, sim_interp.ADDEDMASS, sim_interp.DAMPING, matrix.stiffness,
                    sim_interp.WAVEEX, plot, mufp, coord)

    # Wave Loading Frequency Domain Analysis
    s, _ = jonswap(Hs, Tp, df, gammaJS, f, plot)
    s_mat = np.ones(shape=(len(s), 8))

    for ii in np.arange(0, 8, 1):
        s_mat[:, ii] = s

    # Creating empty variables
    Tz = np.zeros(shape=(3, 8))
    Significant_Amplitude = np.zeros(shape=(3, 8))
    N_mpm = np.zeros(shape=(3, 8))
    MPM = np.zeros(shape=(3, 8))

    # Calculating response spectrum from RAO and excitation spectrum
    response = spectrum_response(RAO, s_mat)

    # Deriving system response characteristics from response spectrum
    Tz[:, :], Significant_Amplitude[:, :], N_mpm[:, :], MPM[:, :], _, _, _ = spectrum_response_statistics(Tr,
                                                                                                      response, f, df)
    return Tz, Significant_Amplitude, N_mpm, MPM