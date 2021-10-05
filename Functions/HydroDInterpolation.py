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
        interp_WAVEEX = np.zeros(shape=(2, int(results[1].numwavelengths), 6, 4))
        interp_motions = np.zeros(shape=(int(results[1].numheadangles), int(results[1].numwavelengths), 6, 4))
        self.front_column_force = np.zeros(results[0].front_column_force.shape) + 0j
        self.left_column_force = np.zeros(results[0].front_column_force.shape) + 0j
        self.right_column_force = np.zeros(results[0].front_column_force.shape) + 0j

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

        for ii in np.linspace(0, len(results[0].front_column_force[:, 0, 0]) - 1, len(results[0].front_column_force[:, 0, 0])).astype(int):
            for jj in np.linspace(0, len(results[0].front_column_force[0, :, 0]) - 1, len(results[0].front_column_force[0, :, 0])).astype(int):
                for kk in np.linspace(0, len(results[0].front_column_force[0, 0, :]) - 1, len(results[0].front_column_force[0, 0, :])).astype(int):
                    self.front_column_force[ii, jj, kk] = np.interp(interp_at[variable], [interp_inputs[0, variable], interp_inputs[1, variable]], [results[0].front_column_force[ii,jj,kk], results[1].front_column_force[ii,jj,kk]])
                    self.left_column_force[ii, jj, kk] = np.interp(interp_at[variable], [interp_inputs[0, variable], interp_inputs[1, variable]], [results[0].left_column_force[ii,jj,kk], results[1].left_column_force[ii,jj,kk]])
                    self.right_column_force[ii,jj, kk] = np.interp(interp_at[variable], [interp_inputs[0, variable], interp_inputs[1, variable]], [results[0].right_column_force[ii,jj,kk], results[1].right_column_force[ii,jj,kk]])


        self.numwavelengths = results[1].numwavelengths
        self.numheadangles = results[1].numheadangles
        self.ADDEDMASS=interp_ADDEDMASS
        self.DAMPING=interp_DAMPING
        self.WAVEEX=interp_WAVEEX
        self.MOTIONS=interp_motions
        self.wave_disc = results[1].wave_disc

class Interpolate2d:
    def __init__(self, results, interp_inputs, interp_at, interp_var):
        sim_interp1 = LinearInterpolate(results[0:2], np.array([interp_inputs[0,:],interp_inputs[1,:]]), interp_at, interp_var[0])
        sim_interp2 = LinearInterpolate(results[2:4], np.array([interp_inputs[2,:],interp_inputs[3,:]]), interp_at, interp_var[0])
        sim_interp = LinearInterpolate([sim_interp1, sim_interp2], np.array([interp_inputs[1,:],interp_inputs[2,:]]), interp_at,  interp_var[1])
        self.numwavelengths = sim_interp.numwavelengths
        self.numheadangles = sim_interp.numheadangles
        self.ADDEDMASS=sim_interp.ADDEDMASS
        self.DAMPING=sim_interp.DAMPING
        self.WAVEEX=sim_interp.WAVEEX
        self.MOTIONS=sim_interp.MOTIONS
        self.front_column_force = sim_interp.front_column_force
        self.left_column_force = sim_interp.left_column_force
        self.right_column_force = sim_interp.right_column_force
        self.wave_disc = sim_interp1.wave_disc

class Interpolate3d:
    def __init__(self, results, interp_inputs, interp_at):
        # Methodology interpolated between each of the design variables
        # 2d interpolate is used to interpolate between x & y, used in interp_func1,2,4,5
        # Linear interp used to interpolate between the x&y interp result and diameter
        # Final interp is between the drafts
        # For first draft, interpolating x, y, and then diameter
        interp_func1 = Interpolate2d(results, interp_inputs, interp_at, [0, 1])
        interp_func2 = Interpolate2d(results[4:8], interp_inputs[4:8, :], interp_at, [0, 1])
        interp_func3 = LinearInterpolate([interp_func1, interp_func2], interp_inputs[[0, 4], :], interp_at, 2)

        self.numwavelengths = interp_func3.numwavelengths
        self.numheadangles = interp_func3.numheadangles
        self.ADDEDMASS=interp_func3.ADDEDMASS
        self.DAMPING=interp_func3.DAMPING
        self.WAVEEX=interp_func3.WAVEEX
        self.MOTIONS=interp_func3.MOTIONS
        self.wave_disc = interp_func3.wave_disc
        self.front_column_force = interp_func3.front_column_force
        self.left_column_force = interp_func3.left_column_force
        self.right_column_force = interp_func3.right_column_force


class Interpolate4d:
    def __init__(self, results, interp_inputs, interp_at):
        # Methodology interpolated between each of the design variables
        # 2d interpolate is used to interpolate between x & y, used in interp_func1,2,4,5
        # Linear interp used to interpolate between the x&y interp result and diameter
        # Final interp is between the drafts
        # For first draft, interpolating x, y, and then diameter
        interp_func1 = Interpolate2d(results, interp_inputs, interp_at, [0, 1])
        interp_func2 = Interpolate2d(results[4:8], interp_inputs[4:8, :], interp_at, [0, 1])
        interp_func3 = LinearInterpolate([interp_func1, interp_func2], interp_inputs[[0, 4], :], interp_at, 2)

        interp_func4 = Interpolate2d(results[8:12], interp_inputs[8:12, :], interp_at, [0, 1])
        interp_func5 = Interpolate2d(results[12:16], interp_inputs[12:16, :], interp_at, [0, 1])
        interp_func6 = LinearInterpolate([interp_func4, interp_func5], interp_inputs[[8, 12], :], interp_at, 2)

        interp_func7 = LinearInterpolate([interp_func3, interp_func6], interp_inputs[[0, 8], :], interp_at, 3)
        self.numwavelengths = interp_func7.numwavelengths
        self.numheadangles = interp_func7.numheadangles
        self.ADDEDMASS=interp_func7.ADDEDMASS
        self.DAMPING=interp_func7.DAMPING
        self.WAVEEX=interp_func7.WAVEEX
        self.MOTIONS=interp_func7.MOTIONS
        self.wave_disc = interp_func7.wave_disc
        self.front_column_force = interp_func7.front_column_force
        self.left_column_force = interp_func7.left_column_force
        self.right_column_force = interp_func7.right_column_force
