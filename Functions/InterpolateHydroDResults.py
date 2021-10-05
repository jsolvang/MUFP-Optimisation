import numpy as np
from scipy import interpolate

class InterpolateResults:
    def __init__(self, sim, f):
        self.numheadangles = len(sim.WAVEEX[:,0,0,0])
        self.wave_disc = np.zeros(shape=(len(f), 5))
        self.WAVEEX = np.zeros(shape=(int(self.numheadangles), len(f), 6, 4))
        self.MOTIONS = np.zeros(shape=(int(self.numheadangles), len(f), 6, 4))
        self.DAMPING = np.zeros(shape=(len(f), 6, 6))
        self.ADDEDMASS = np.zeros(shape=(len(f), 6, 6))
        self.front_column_force = np.zeros(shape=(2, len(f), 9)) + 0j
        self.left_column_force = np.zeros(shape=(2, len(f), 9)) + 0j
        self.right_column_force = np.zeros(shape=(2, len(f), 9)) + 0j
        self.DOF = [0 ,1 ,2 ,3 ,4 ,5]
        self.directions = [0 ,1]

        self.wave_disc[:, 4] = f


        for direct in self.directions:
            for dof in self.DOF:
                for jj in np.linspace(0,  len(self.WAVEEX[0, 0, 0, :]) -1,  len(self.WAVEEX[0, 0, 0, :])).astype(int):
                    self.WAVEEX[direct, :, dof, jj] = np.interp(f, sim.wave_disc[:, 4], sim.WAVEEX[direct, :, dof, jj])
                    self.MOTIONS[direct, :, dof, jj] = np.interp(f, sim.wave_disc[:, 4], sim.MOTIONS[direct, :, dof, jj])

        for ii in np.linspace(0, len(self.DAMPING[0, :, 0]) - 1, len(self.DAMPING[0, :, 0])).astype(int):
            for jj in np.linspace(0, len(self.DAMPING[0, 0, :]) - 1, len(self.DAMPING[0, :, 0])).astype(int):
                self.DAMPING[:, ii, jj] = np.interp(f, sim.wave_disc[:, 4], sim.DAMPING[:, ii, jj])
                self.ADDEDMASS[:, ii, jj] = np.interp(f, sim.wave_disc[:, 4], sim.DAMPING[:, ii, jj])

        for ii in np.linspace(0, len(sim.front_column_force[:, 0, 0]) - 1, len(sim.front_column_force[:, 0, 0])).astype(int):
            for jj in np.linspace(0, len(sim.front_column_force[0, 0, :]) - 1, len(sim.front_column_force[0, 0, :])).astype(int):
                self.front_column_force[ii,:, jj] = np.interp(f, sim.wave_disc[:, 4], sim.front_column_force[ii,:, jj])
                self.left_column_force[ii,:, jj] = np.interp(f, sim.wave_disc[:, 4], sim.left_column_force[ii,:, jj])
                self.right_column_force[ii,:, jj] = np.interp(f, sim.wave_disc[:, 4], sim.right_column_force[ii,:, jj])
