import numpy as np
import subprocess
import os
import matplotlib.pyplot as plt
import pickle


class LinearInterpolate:
    def __init__(self, results, interp_inputs, interp_at, variable):
        # Variable = desired interpolation variable, 0 = X, 1 = Y, 2 = Diameter
        # Function to accept two "results"
        interp_ADDEDMASS = np.zeros(shape=(int(results[1].numwavelengths), 6, 6))
        interp_DAMPING = np.zeros(shape=(int(results[1].numwavelengths), 6, 6))
        interp_WAVEEX = np.zeros(shape=(2, int(results[1].numwavelengths), 6, 6))
        interp_motions = np.zeros(shape=(int(results[1].numheadangles), int(results[1].numwavelengths), 6, 4))

        for nn in np.linspace(0, int(results[1].numwavelengths) - 1, int(results[1].numwavelengths)).astype(int):
            for mm in np.linspace(0, 5, 6).astype(int):
                for kk in np.linspace(0, 5, 6).astype(int):
                    interp_ADDEDMASS[nn, mm, kk] = np.interp(interp_at[variable],
                                                             [interp_inputs[0, variable], interp_inputs[1, variable]],
                                                             [results[0].ADDEDMASS[nn, mm, kk],
                                                              results[1].ADDEDMASS[nn, mm, kk]])
                    interp_DAMPING[nn, mm, kk] = np.interp(interp_at[variable],
                                                           [interp_inputs[0, variable], interp_inputs[1, variable]],
                                                           [results[0].DAMPING[nn, mm, kk],
                                                            results[1].DAMPING[nn, mm, kk]])

        for nn in np.linspace(0, int(results[0].numwavelengths) - 1, int(results[0].numwavelengths)).astype(int):
            for mm in np.linspace(0, 5, 6).astype(int):
                for kk in np.linspace(0, 3, 4).astype(int):
                    interp_WAVEEX[0, nn, mm, kk] = np.interp(interp_at[variable],
                                                             [interp_inputs[0, variable], interp_inputs[1, variable]],
                                                             [results[0].WAVEEX[0, nn, mm, kk],
                                                              results[1].WAVEEX[0, nn, mm, kk]])
                    interp_WAVEEX[1, nn, mm, kk] = np.interp(interp_at[variable],
                                                             [interp_inputs[0, variable], interp_inputs[1, variable]],
                                                             [results[0].WAVEEX[1, nn, mm, kk],
                                                              results[1].WAVEEX[1, nn, mm, kk]])

        for pp in np.linspace(0, len(results[0].MOTIONS[:,0,0,0])-1, len(results[0].MOTIONS[:,0,0,0])).astype(int):
            for nn in np.linspace(0, int(results[0].numwavelengths) - 1, int(results[0].numwavelengths)).astype(int):
                for mm in np.linspace(0, 5, 6).astype(int):
                    for kk in np.linspace(0, 3, 4).astype(int):
                        interp_motions[pp, nn, mm, kk] = np.interp(interp_at[variable],
                                                                 [interp_inputs[0, variable], interp_inputs[1, variable]],
                                                                 [results[0].MOTIONS[pp, nn, mm, kk],
                                                                  results[1].MOTIONS[pp, nn, mm, kk]])
        self.numwavelengths = results[1].numwavelengths
        self.numheadangles = results[1].numheadangles
        self.ADDEDMASS=interp_ADDEDMASS
        self.DAMPING=interp_DAMPING
        self.WAVEEX=interp_WAVEEX
        self.MOTIONS=interp_motions
        self.wave_disc = results[1].wave_disc

class Interpolate2d:
    def __init__(self, results, interp_inputs, interp_at):
        sim_interp1 = LinearInterpolate(results[0:2], interp_inputs, interp_at, 1)
        sim_interp2 = LinearInterpolate(results[2:4], interp_inputs, interp_at, 1)
        sim_interp = LinearInterpolate([sim_interp1,sim_interp2], interp_inputs, interp_at, 2)
        self.numwavelengths = sim_interp.numwavelengths
        self.numheadangles = sim_interp.numheadangles
        self.ADDEDMASS=sim_interp.ADDEDMASS
        self.DAMPING=sim_interp.DAMPING
        self.WAVEEX=sim_interp.WAVEEX
        self.MOTIONS=sim_interp.MOTIONS
        self.wave_disc = sim_interp1.wave_disc