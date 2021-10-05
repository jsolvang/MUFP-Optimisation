from RAO_Calculation import calulate_RAOs
from SpectrumStatistics import spectrum_response_statistics, spectrum_response
from InterpolateHydroDResults import InterpolateResults
from WindForceSpectrum import WindForceSpectrum
import numpy as np


def windspectrumanalysis(sim, matrix, plot, mufp, coord, df, f, Tr):
    # Pulling Wind Force Spectrum from SIMA simulation
    f_xx, sx_wind, sy_wind, smz_wind = WindForceSpectrum()

    # Interpolating loading spectrum to evaluated discretisation
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

    # Declaring variable unit force of unit [N/N] and [Nm/m]
    unit_force = np.ones(shape=(2, len(f), 6, 4)) + 0j
    unit_force[:, :, 3:5, 0] = mufp.hub_height + 0j
    unit_force[:, :, 3:5, 1] = mufp.hub_height + 0j

    # Interpolating Results to new frequency range
    sim_interp = InterpolateResults(sim, f)

    # Calculating RAO per unit force
    RAO_per_unit_force, Y, H, _ = calulate_RAOs(sim_interp.wave_disc, matrix.mass, sim_interp.ADDEDMASS, sim_interp.DAMPING,
                                       matrix.stiffness,
                                       unit_force, plot, mufp, coord)

    # Calculating wind response spectrum
    wind_response = spectrum_response(RAO_per_unit_force[0,:,:,:], s_wind)

    # Creating shell variables
    Tz = np.zeros(shape=(3, 8))
    Significant_Amplitude = np.zeros(shape=(3, 8))
    N_mpm = np.zeros(shape=(3, 8))
    MPM = np.zeros(shape=(3, 8))

    # Spectral Moment Analysis
    Tz[:, :], Significant_Amplitude[:, :], N_mpm[:, :], MPM[:, :], _, _, _ = spectrum_response_statistics(Tr, wind_response, f, df)

    return Significant_Amplitude, MPM
