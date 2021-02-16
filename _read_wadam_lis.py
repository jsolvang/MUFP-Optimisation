import numpy as np


class ReadWadamLis:

    def __init__(self, HydroD_result):
        self.path = HydroD_result + "\WADAM1.LIS"
        self._extract_results()

    def _extract_results(self):
        hydroD_results = open(self.path, 'r').readlines()
        ii = 1
        found = 0

        while found == 0:
            k = hydroD_results[ii].find('2.8 ENVIRONMENTAL DATA:\n')
            ii += 1
            if k != -1:
                found = 1
                self.waterdepth = float(hydroD_results[ii + 2].split()[3])
                self.numwavelengths = float(hydroD_results[ii + 3].split()[5])
                self.numheadangles = float(hydroD_results[ii + 4].split()[5])
                self.wavedata = np.zeros(shape=(int(self.numwavelengths), int(self.numheadangles), 5))

                for jj in np.linspace(1, int(self.numheadangles), int(self.numheadangles)).astype(int) - 1:
                    for kk in np.linspace(1, int(self.numwavelengths), int(self.numwavelengths)).astype(int) - 1:
                        self.wavedata[kk, jj, :] = np.float_(hydroD_results[ii + kk + 12].split())

        ii = 1
        found = 0

        self.wave_disc = np.zeros(shape=(int(self.numwavelengths), 5))

        while found == 0:
            k = hydroD_results[ii].find('     WAVE DESCRIPTION:\n')
            ii += 1
            if k != -1:
                found = 1
                for jj in np.linspace(0, int(self.numwavelengths) -1, int(self.numwavelengths)).astype(int):
                    self.wave_disc[jj, :] = np.float_(hydroD_results[ii + 4 + jj].split())

        ii = 1
        found = 0
        while found == 0:
            k = hydroD_results[ii].find(
                '    THE OUTPUT IS NON-DIMENSIONALIZED USING -\n')
            ii += 1
            if k != -1:
                found = 1
                self.RO = float(hydroD_results[ii + 7].split()[2])
                self.G = float(hydroD_results[ii + 8].split()[2])
                self.VOL = float(hydroD_results[ii + 9].split()[2])
                self.L = float(hydroD_results[ii + 10].split()[2])
                self.WA = float(hydroD_results[ii + 11].split()[2])

        ii = 0

        self.ADDMASS = np.zeros(shape=(int(self.numwavelengths), 6, 7))
        for nn in np.linspace(0, int(self.numwavelengths) - 1, int(self.numwavelengths)).astype(int):
            found = 0
            while found == 0 and ii < len(hydroD_results):
                k = hydroD_results[ii].find(
                    '    ADDED MASS MATRIX                                                                                                              \n')
                ii += 1
                if k != -1:
                    found = 1
                    self.ADDMASS[nn, 0, :] = np.float_(hydroD_results[ii + 3].split())
                    self.ADDMASS[nn, 1, :] = np.float_(hydroD_results[ii + 4].split())
                    self.ADDMASS[nn, 2, :] = np.float_(hydroD_results[ii + 5].split())
                    self.ADDMASS[nn, 3, :] = np.float_(hydroD_results[ii + 6].split())
                    self.ADDMASS[nn, 4, :] = np.float_(hydroD_results[ii + 7].split())
                    self.ADDMASS[nn, 5, :] = np.float_(hydroD_results[ii + 8].split())
                    #print(self.ADDMASS[nn,0,:])
                    break
        ii = 0
        self.ADDEDMASS = self.ADDMASS[:, :, 1:7]
        self.ADDEDMASS[:, 0:3, 0:3] = self.ADDEDMASS[:, 0:3, 0:3] * self.RO * self.VOL
        self.ADDEDMASS[:, 0:3, 3:6] = self.ADDEDMASS[:, 0:3, 3:6] * self.RO * self.VOL * self.L
        self.ADDEDMASS[:, 3:6, 0:3] = self.ADDEDMASS[:, 3:6, 0:3] * self.RO * self.VOL * self.L
        self.ADDEDMASS[:, 3:6, 3:6] = self.ADDEDMASS[:, 3:6, 3:6] * self.RO * self.VOL * self.L * self.L

        self.DAMPING = np.zeros(shape=(int(self.numwavelengths), 6, 7))
        for nn in np.linspace(0, int(self.numwavelengths) - 1, int(self.numwavelengths)).astype(int):
            found = 0
            while found == 0 and ii < len(hydroD_results):
                k = hydroD_results[ii].find(
                    '    DAMPING MATRIX                                                                                                                 \n')
                ii += 1
                if k != -1:
                    found = 1;
                    self.DAMPING[nn, 0, :] = np.float_(hydroD_results[ii + 3].split())
                    self.DAMPING[nn, 1, :] = np.float_(hydroD_results[ii + 4].split())
                    self.DAMPING[nn, 2, :] = np.float_(hydroD_results[ii + 5].split())
                    self.DAMPING[nn, 3, :] = np.float_(hydroD_results[ii + 6].split())
                    self.DAMPING[nn, 4, :] = np.float_(hydroD_results[ii + 7].split())
                    self.DAMPING[nn, 5, :] = np.float_(hydroD_results[ii + 8].split())

        ii = 1

        self.DAMPING = self.DAMPING[:,:,1:7]

        self.DAMPING[:, 0: 3, 0: 3] = self.DAMPING[:, 0: 3, 0: 3] * self.RO * self.VOL * np.sqrt(np.divide(self.G, self.L))
        self.DAMPING[:, 0: 3, 3: 6] = self.DAMPING[:, 0: 3, 3: 6] * self.RO * self.VOL * np.sqrt(self.G * self.L)
        self.DAMPING[:, 3: 6, 0: 3] = self.DAMPING[:, 3: 6, 0: 3] * self.RO * self.VOL * np.sqrt(self.G * self.L)
        self.DAMPING[:, 3: 6, 3: 6] = self.DAMPING[:, 3: 6, 3: 6] * self.RO * self.VOL * self.L * np.sqrt(self.G * self.L)

        self.WAVEEX = np.zeros(shape=(int(self.numheadangles), int(self.numwavelengths), 6, 4))
        self.MOTIONS = np.zeros(shape=(int(self.numheadangles), int(self.numwavelengths), 6, 4))

        for nn in np.linspace(0, int(self.numwavelengths) - 1, int(self.numwavelengths)).astype(int):
            for mm in np.linspace(0, int(self.numheadangles) - 1, int(self.numheadangles)).astype(int):
                found = 0
                while found == 0 and ii < len(hydroD_results):
                    k = hydroD_results[ii].find(
                        "    EXCITING FORCES AND MOMENTS FROM THE HASKIND'S RELATIONS                                                                       \n")
                    ii += 1
                    if k != -1:
                        found = 1
                        for qq in np.linspace(0, 5, 6).astype(int):
                            self.WAVEEX[mm, nn, qq, :] = np.float_(hydroD_results[ii + 4 + 2 * qq].split()[1:5])
                            self.MOTIONS[mm, nn, qq, :] = np.float_(hydroD_results[ii + 22 + 2 * qq].split()[1:5])
