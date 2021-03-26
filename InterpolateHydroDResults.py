import numpy as np


class InterpolateResults:
    def __init__(self, wave_disc_og, WAVEEX_og, MOTIONS_og, ADDEDMASS_og, DAMPING_og, f):
        self.numheadangles = len(WAVEEX_og[:,0,0,0])
        self.wave_disc = np.zeros(shape=(len(f), 5))
        self.WAVEEX = np.zeros(shape=(int(self.numheadangles), len(f), 6, 4))
        self.MOTIONS = np.zeros(shape=(int(self.numheadangles), len(f), 6, 4))
        self.DAMPING = np.zeros(shape=(len(f), 6, 6))
        self.ADDEDMASS = np.zeros(shape=(len(f), 6, 6))
        self.DOF = [0 ,1 ,2 ,3 ,4 ,5]
        self.directions = [0 ,1]

        for jj in np.linspace(0, len(self.wave_disc[0, :] ) -1, len(self.wave_disc[0, :])).astype(int):
            self.wave_disc[:, jj] = np.interp(f, wave_disc_og[:, 4], wave_disc_og[:, jj])

        for direct in self.directions:
            for dof in self.DOF:
                for jj in np.linspace(0,  len(self.WAVEEX[0, 0, 0, :]) -1,  len(self.WAVEEX[0, 0, 0, :])).astype(int):
                    self.WAVEEX[direct, :, dof, jj] = np.interp(f, wave_disc_og[:, 4], WAVEEX_og[direct, :, dof, jj])
                    self.MOTIONS[direct, :, dof, jj] = np.interp(f, wave_disc_og[:, 4], MOTIONS_og[direct, :, dof, jj])

        for ii in np.linspace(0, len(self.DAMPING[0, :, 0]) - 1, len(self.DAMPING[0, :, 0])).astype(int):
            for jj in np.linspace(0, len(self.DAMPING[0, 0, :]) - 1, len(self.DAMPING[0, :, 0])).astype(int):
                self.DAMPING[:, ii, jj] = np.interp(f, wave_disc_og[:, 4], DAMPING_og[:, ii, jj])
                self.ADDEDMASS[:, ii, jj] = np.interp(f, wave_disc_og[:, 4], ADDEDMASS_og[:, ii, jj])