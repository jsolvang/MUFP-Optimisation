from RAO_Calculation import calulate_RAOs
from SpectrumStatistics import spectrum_response_statistics, spectrum_response
from InterpolateHydroDResults import InterpolateResults

import numpy as np


def windspectrumanalysis(sim, matrix, plot, mufp, coord, df, f, Tr, s_wind):
    # Declaring variable unit force of unit [N/N] and [Nm/m]
    unit_force = np.ones(shape=(2, len(f), 6, 4))
    unit_force[:, :, 3:6, 0] = mufp.hub_height
    unit_force[:, :, 3:6, 1] = mufp.hub_height

    # Interpolating Results to new frequency range
    sim_interp = InterpolateResults(sim.wave_disc, sim.WAVEEX, sim.MOTIONS, sim.ADDEDMASS, sim.DAMPING, f)

    # Calculating RAO per unit force
    RAO_per_unit_force, Y, H, _ = calulate_RAOs(sim_interp.wave_disc, matrix.mass, sim_interp.ADDEDMASS, sim_interp.DAMPING,
                                       matrix.stiffness,
                                       unit_force, plot, mufp, coord)

    # Calculating wind response spectrum
    wind_response = spectrum_response(RAO_per_unit_force, s_wind)

    # Creating shell variables
    Tz = np.zeros(shape=(3, 8))
    Significant_Amplitude = np.zeros(shape=(3, 8))
    N_mpm = np.zeros(shape=(3, 8))
    MPM = np.zeros(shape=(3, 8))

    # Spectral Moment Analysis
    Tz[:, :], Significant_Amplitude[:, :], N_mpm[:, :], MPM[:, :], _, _, _ = spectrum_response_statistics(Tr, wind_response, f, df)

    return Tz, Significant_Amplitude, N_mpm, MPM, wind_response, RAO_per_unit_force, Y, H
