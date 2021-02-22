import numpy as np


class ReadWadamLis:

    def __init__(self, HydroD_result):
        self.path = HydroD_result + "\WADAM1.LIS"
        self._extract_results()

    def _extract_results(self):
        hydroD_results = open(self.path, 'r').readlines()

        ii = 1
        found = 0
        self.mass_matrix = np.zeros(shape=(6, 6))

        while found == 0:
            k = hydroD_results[ii].find('     MASS PROPERTIES AND STRUCTURAL DATA:\n')
            ii += 1
            if k != -1:
                found = 1
                self.mass = float(hydroD_results[ii+2].split()[6])
                self.COG = [float(hydroD_results[ii+4].split()[5]), float(hydroD_results[ii+5].split()[2]), float(hydroD_results[ii+6].split()[2])]
                self.RoG = [float(hydroD_results[ii+7].split()[6]), float(hydroD_results[ii+8].split()[6]), float(hydroD_results[ii+9].split()[6])]
                self.CFM = [float(hydroD_results[ii+10].split()[4].replace("=","")), float(hydroD_results[ii+11].split()[5]), float(hydroD_results[ii+12].split()[5])]

        ii = 1
        found = 0

        self.mass_matrix[0, 0] = self.mass
        self.mass_matrix[1, 1] = self.mass
        self.mass_matrix[2, 2] = self.mass
        self.mass_matrix[3, 3] = self.mass*np.square(self.RoG[0])
        self.mass_matrix[4, 4] = self.mass*np.square(self.RoG[1])
        self.mass_matrix[5, 5] = self.mass*np.square(self.RoG[2])

        self.mass_matrix[3, 4] = -self.mass * self.CFM[0]
        self.mass_matrix[4, 3] = -self.mass * self.CFM[0]

        self.mass_matrix[3, 5] = -self.mass * self.CFM[1]
        self.mass_matrix[5, 3] = -self.mass * self.CFM[1]

        self.mass_matrix[0, 4] = self.mass*self.COG[2]
        self.mass_matrix[4, 0] = self.mass*self.COG[2]

        self.mass_matrix[1, 3] = -self.mass*self.COG[2]
        self.mass_matrix[3, 1] = -self.mass*self.COG[2]

        self.mass_matrix[1, 5] = self.mass*self.COG[0]
        self.mass_matrix[5, 1] = self.mass*self.COG[0]

        self.mass_matrix[2, 4] = -self.mass*self.COG[0]
        self.mass_matrix[4, 2] = -self.mass*self.COG[0]



        self.stiffness_matrix =  np.zeros(shape=(6, 6))

        while found == 0:
            k = hydroD_results[ii].find('     HYDROSTATIC DATA:\n')
            ii += 1
            if k != -1:
                found = 1
                print(ii)
                self.displaced_volume = float(hydroD_results[ii+3].split()[6])
                self.WPA = float(hydroD_results[ii+5].split()[5])
                self.COB = [float(hydroD_results[ii+7].split()[5]), float(hydroD_results[ii+8].split()[2]),
                            float(hydroD_results[ii+9].split()[1].replace("=", ""))]

                self.meta_centric_height_transverse = float(hydroD_results[ii+10].split()[5])
                self.meta_centric_height_longitudinal = float(hydroD_results[ii+11].split()[5])

                self.stiffness_matrix[2, 2] = float(hydroD_results[ii+12].split()[5])
                self.stiffness_matrix[2, 3] = float(hydroD_results[ii+13].split()[5])
                self.stiffness_matrix[2, 4] = float(hydroD_results[ii+14].split()[4].replace("=", ""))
                self.stiffness_matrix[4, 2] = self.stiffness_matrix[2, 4]
                self.stiffness_matrix[3, 3] = float(hydroD_results[ii+15].split()[5])
                self.stiffness_matrix[4, 4] = float(hydroD_results[ii+16].split()[5])
                self.stiffness_matrix[3, 4] = float(hydroD_results[ii+17].split()[5])
                self.stiffness_matrix[3, 5] = float(hydroD_results[ii+18].split()[5])
                self.stiffness_matrix[4, 5] = float(hydroD_results[ii+19].split()[5])

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
