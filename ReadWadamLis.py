import numpy as np


class ReadWadamLis:

    def __init__(self, HydroD_result, sectional_loads, reflection, fixed):
        self.path = HydroD_result + "\WADAM1.LIS"
        self._extract_results(sectional_loads, reflection, fixed)

    def _extract_results(self, sectional_loads, reflection, fixed):
        hydroD_results = open(self.path, 'r').readlines()
        ii = 1
        found = 0
        self.mass_matrix = np.zeros(shape=(6, 6))
        self.COG = np.zeros(3)
        self.CFM = np.zeros(3)

        while found == 0:
            k = hydroD_results[ii].find('     MASS PROPERTIES AND STRUCTURAL DATA:\n')
            ii += 1
            if k != -1:
                found = 1
                self.mass = float(hydroD_results[ii + 2].split()[6])
                if hydroD_results[ii + 4].split()[4] == "=":
                    self.COG[0] = float(hydroD_results[ii + 4].split()[5])
                else:
                    self.COG[0] = float(hydroD_results[ii + 4].split()[4].replace("=", ""))
                if hydroD_results[ii + 5].split()[1] == "=":
                    self.COG[1] = float(hydroD_results[ii + 5].split()[2])
                else:
                    self.COG[1] = float(hydroD_results[ii + 5].split()[1].replace("=", ""))
                if hydroD_results[ii + 6].split()[1] == "=":
                    self.COG[2] = float(hydroD_results[ii + 6].split()[2])
                else:
                    self.COG[2] = float(hydroD_results[ii + 6].split()[1].replace("=", ""))

                self.RoG = [float(hydroD_results[ii + 7].split()[6]), float(hydroD_results[ii + 8].split()[6]),
                            float(hydroD_results[ii + 9].split()[6])]

                if hydroD_results[ii + 10].split()[4] == "=":
                    self.CFM[0] = float(hydroD_results[ii + 10].split()[5])
                else:
                    self.CFM[0] = float(hydroD_results[ii + 10].split()[4].replace("=", ""))
                if hydroD_results[ii + 11].split()[4] == "=":
                    self.CFM[1] = float(hydroD_results[ii + 11].split()[5])
                else:
                    self.CFM[1] = float(hydroD_results[ii + 11].split()[4].replace("=", ""))
                if hydroD_results[ii + 12].split()[4] == "=":
                    self.CFM[2] = float(hydroD_results[ii + 12].split()[5])
                else:
                    self.CFM[2] = float(hydroD_results[ii + 12].split()[4].replace("=", ""))

        if fixed == 0:
            ii = 1
            found = 0
            while found == 0:
                k = hydroD_results[ii].find(' 4.2 STATIC  RESULTS\n')
                ii += 1
                if k != -1:
                    found = 1
                    M1 = np.float_(hydroD_results[ii + 8].split()[1:7])
                    M2 = np.float_(hydroD_results[ii + 9].split()[1:7])
                    M3 = np.float_(hydroD_results[ii + 10].split()[1:7])
                    M4 = np.float_(hydroD_results[ii + 11].split()[1:7])
                    M5 = np.float_(hydroD_results[ii + 12].split()[1:7])
                    M6 = np.float_(hydroD_results[ii + 13].split()[1:7])
                    self.mass_matrix = np.array([M1, M2, M3, M4, M5, M6])

                    K1 = np.float_(hydroD_results[ii + 20].split()[1:7])
                    K2 = np.float_(hydroD_results[ii + 21].split()[1:7])
                    K3 = np.float_(hydroD_results[ii + 22].split()[1:7])
                    K4 = np.float_(hydroD_results[ii + 23].split()[1:7])
                    K5 = np.float_(hydroD_results[ii + 24].split()[1:7])
                    K6 = np.float_(hydroD_results[ii + 25].split()[1:7])
                    self.stiffness_matrix = np.array([K1, K2, K3, K4, K5, K6])

            ii = 1
            found = 0

            self.COB = np.zeros(3)

            while found == 0:
                k = hydroD_results[ii].find('     HYDROSTATIC DATA:\n')
                ii += 1
                if k != -1:
                    found = 1
                    self.displaced_volume = float(hydroD_results[ii + 3].split()[6])
                    self.WPA = float(hydroD_results[ii + 5].split()[5])
                    if hydroD_results[ii + 7].split()[4] == "=":
                        self.COB[0] = float(hydroD_results[ii + 7].split()[5])
                    else:
                        self.COB[0] = float(hydroD_results[ii + 7].split()[4].replace("=", ""))
                    if hydroD_results[ii + 8].split()[1] == "=":
                        self.COB[1] = float(hydroD_results[ii + 8].split()[2])
                    else:
                        self.COB[1] = float(hydroD_results[ii + 8].split()[1].replace("=", ""))
                    if hydroD_results[ii + 9].split()[1] == "=":
                        self.COB[2] = float(hydroD_results[ii + 9].split()[2])
                    else:
                        self.COB[2] = float(hydroD_results[ii + 9].split()[1].replace("=", ""))

                    if hydroD_results[ii + 10].split()[4] == "=":
                        self.meta_centric_height_transverse = float(hydroD_results[ii + 10].split()[5])
                    else:
                        self.meta_centric_height_transverse = float(hydroD_results[ii + 10].split()[4].replace("=", ""))
                    if hydroD_results[ii + 11].split()[4] == "=":
                        self.meta_centric_height_longitudinal = float(hydroD_results[ii + 11].split()[5])
                    else:
                        self.meta_centric_height_longitudinal = float(hydroD_results[ii + 11].split()[4].replace("=", ""))

        ii = 1
        found = 0

        while found == 0:
            k = hydroD_results[ii].find('2.8 ENVIRONMENTAL DATA:\n')
            ii += 1
            if k != -1:
                found = 1
                self.waterdepth = float(hydroD_results[ii + 2].split()[3])
                self.numwavelengths = int(hydroD_results[ii + 3].split()[5])
                self.numheadangles = int(hydroD_results[ii + 4].split()[5])
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
                for jj in np.linspace(0, int(self.numwavelengths) - 1, int(self.numwavelengths)).astype(int):
                    self.wave_disc[jj, :] = np.float_(hydroD_results[ii + 4 + jj].split())

        ii = 1
        found = 0
        while found == 0 and ii < len(hydroD_results):
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

        ADDMASS = np.zeros(shape=(int(self.numwavelengths), 6, 7))
        for nn in np.linspace(0, int(self.numwavelengths) - 1, int(self.numwavelengths)).astype(int):
            found = 0
            while found == 0 and ii < len(hydroD_results):
                k = hydroD_results[ii].find(
                    '    ADDED MASS MATRIX                                                                                                              \n')
                ii += 1
                if k != -1:
                    found = 1
                    ADDMASS[nn, 0, :] = np.float_(hydroD_results[ii + 3].split())
                    ADDMASS[nn, 1, :] = np.float_(hydroD_results[ii + 4].split())
                    ADDMASS[nn, 2, :] = np.float_(hydroD_results[ii + 5].split())
                    ADDMASS[nn, 3, :] = np.float_(hydroD_results[ii + 6].split())
                    ADDMASS[nn, 4, :] = np.float_(hydroD_results[ii + 7].split())
                    ADDMASS[nn, 5, :] = np.float_(hydroD_results[ii + 8].split())
                    # print(self.ADDMASS[nn,0,:])
                    break
        ii = 0
        self.ADDEDMASS = ADDMASS[:, :, 1:7]


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

        self.DAMPING = self.DAMPING[:, :, 1:7]


        self.WAVEEX = np.zeros(shape=(int(self.numheadangles), int(self.numwavelengths), 6, 4))
        self.MOTIONS = np.zeros(shape=(int(self.numheadangles), int(self.numwavelengths), 6, 4))

        for nn in np.linspace(0, int(self.numwavelengths) - 1, int(self.numwavelengths)).astype(int):
            for mm in np.linspace(0, int(self.numheadangles) - 1, int(self.numheadangles)).astype(int):
                found = 0
                while found == 0 and ii < len(hydroD_results):
                    k = hydroD_results[ii].find(
                        "    EXCITING FORCES AND MOMENTS FROM INTEGRATION OF PRESSURE                                                                       \n")
                    ii += 1
                    if k != -1:
                        found = 1
                        for qq in np.linspace(0, 5, 6).astype(int):
                            self.WAVEEX[mm, nn, qq, :] = np.float_(hydroD_results[ii + 4 + 2 * qq].split()[1:5])
                            self.MOTIONS[mm, nn, qq, :] = np.float_(hydroD_results[ii + 22 + 2 * qq].split()[1:5])

        if fixed == 0:
            self.stiffness_matrix[0: 3, 0: 3] = self.stiffness_matrix[0: 3, 0: 3] * self.RO * self.VOL * self.G / self.L
            self.stiffness_matrix[3: 6, 0: 3] = self.stiffness_matrix[3: 6, 0: 3] * self.RO * self.VOL * self.G
            self.stiffness_matrix[0: 3, 3: 6] = self.stiffness_matrix[0: 3, 3: 6] * self.RO * self.VOL * self.G
            self.stiffness_matrix[3: 6, 3: 6] = self.stiffness_matrix[3: 6, 3: 6] * self.RO * self.VOL * self.G * self.L

        self.mass_matrix[0: 3, 0: 3] = self.mass_matrix[0: 3, 0: 3] * self.RO * self.VOL
        self.mass_matrix[3: 6, 0: 3] = self.mass_matrix[3: 6, 0: 3] * self.RO * self.VOL * self.L
        self.mass_matrix[0: 3, 3: 6] = self.mass_matrix[0: 3, 3: 6] * self.RO * self.VOL * self.L
        self.mass_matrix[3: 6, 3: 6] = self.mass_matrix[3: 6, 3: 6] * self.RO * self.VOL * self.L * self.L

        self.ADDEDMASS[:, 0:3, 0:3] = self.ADDEDMASS[:, 0:3, 0:3] * self.RO * self.VOL
        self.ADDEDMASS[:, 0:3, 3:6] = self.ADDEDMASS[:, 0:3, 3:6] * self.RO * self.VOL * self.L
        self.ADDEDMASS[:, 3:6, 0:3] = self.ADDEDMASS[:, 3:6, 0:3] * self.RO * self.VOL * self.L
        self.ADDEDMASS[:, 3:6, 3:6] = self.ADDEDMASS[:, 3:6, 3:6] * self.RO * self.VOL * self.L * self.L

        self.DAMPING[:, 0: 3, 0: 3] = self.DAMPING[:, 0: 3, 0: 3] * self.RO * self.VOL * np.sqrt(
            np.divide(self.G, self.L))
        self.DAMPING[:, 0: 3, 3: 6] = self.DAMPING[:, 0: 3, 3: 6] * self.RO * self.VOL * np.sqrt(self.G * self.L)
        self.DAMPING[:, 3: 6, 0: 3] = self.DAMPING[:, 3: 6, 0: 3] * self.RO * self.VOL * np.sqrt(self.G * self.L)
        self.DAMPING[:, 3: 6, 3: 6] = self.DAMPING[:, 3: 6, 3: 6] * self.RO * self.VOL * self.L * np.sqrt(
            self.G * self.L)

        self.WAVEEX[:, :, 0: 3, 0: 3] = self.WAVEEX[:, :, 0: 3, 0: 3] * self.RO * self.VOL * self.G * np.divide(self.WA, self.L)
        self.WAVEEX[:, :, 3: 6, 0: 3] = self.WAVEEX[:, :, 3: 6, 0: 3] * self.RO * self.VOL * self.G * self.WA

        self.MOTIONS[:, :, 0: 3, 0: 3] = self.MOTIONS[:, :, 0: 3, 0: 3] * self.WA
        self.MOTIONS[:, :, 3: 6, 0: 3] = self.MOTIONS[:, :, 3: 6, 0: 3] * np.divide(self.WA, self.L)

        ii = 1
        found = 0

        while found == 0:
            k = hydroD_results[ii].find('     SUMMARY FROM THE PANEL MODEL GENERATOR FOR BODY NUMBER :   1\n')
            ii += 1
            if k != -1:
                found = 1
                self.num_wet_panels = int(hydroD_results[ii + 1].split()[-1])
                # self.num_panels = int(hydroD_results[ii + 3].split()[10])
                # self.num_dry_panels = int(hydroD_results[ii + 4].split()[8])

        if reflection == 1:
            self.panel_pressure = np.zeros(shape=(self.numheadangles, self.numwavelengths, 2 * self.num_wet_panels, 15)) + 0j
        else:
            self.panel_pressure = np.zeros(shape=(self.numheadangles, self.numwavelengths, self.num_wet_panels, 15)) + 0j

        self.panel_disc = np.zeros(shape=(5, self.num_wet_panels))

        found = 0
        kk = 0

        while found == 0:
            k = hydroD_results[ii].find(
                '   PI   PANO IVER X         Y         Z          NX      NY      NZ     XC        YC        ZC        AREA      DIAG\n')
            ii += 1
            if k != -1:
                ii += 3
                for jj in np.arange(0, 15, 1):
                    if hydroD_results[ii + 5 * jj] == '1\n':
                        break
                    self.panel_disc[0, kk] = hydroD_results[ii + 5 * jj].split()[0]
                    self.panel_disc[1, kk] = hydroD_results[ii + 5 * jj].split()[-2]
                    self.panel_disc[2, kk] = float(''.join(list(hydroD_results[ii + 5 * jj])[47:70]).split()[0])
                    self.panel_disc[3, kk] = float(''.join(list(hydroD_results[ii + 5 * jj])[47:70]).split()[1])
                    self.panel_disc[4, kk] = float(''.join(list(hydroD_results[ii + 5 * jj])[47:70]).split()[2])
                    kk += 1
                    if int(hydroD_results[ii + 5 * jj].split()[0]) == self.num_wet_panels:
                        found = 1
                        break

        found = 0
        kk = 0
        ii = 0

        ###### Array indexes : [numwavedir ,numwavefreq ,panels ,panel data]
        ######  panel data:  [PI, PANO, SET, IND,  XC, YC, ZC, REAL PART, IMAG.PART ,  ABS.VALUE   PHASE ANGLE(DEG),Panel Area, Real Force, Imag Force, Abs Force]

        for bb in np.arange(0, self.numwavelengths, 1):
            for aa in np.arange(0, self.numheadangles, 1):
                kk = 0
                found = 0
                while found == 0:
                    k = hydroD_results[ii].find('BASIC PART  ')
                    ii += 1
                    if k != -1:
                        found = 1
                        for jj in np.arange(0, 41, 1):
                            if hydroD_results[ii + 5 + jj] == '1\n':
                                break
                            self.panel_pressure[aa, bb, kk, 0:11] =  np.float_(hydroD_results[ii + 5 + jj].split())
                            self.panel_pressure[aa, bb, kk, 11] = self.panel_disc[1, int(self.panel_pressure[aa, bb, kk, 0]) - 1]
                            # [real + 1j*[imag] ] * area * normal direction
                            self.panel_pressure[aa, bb, kk, 12] = (self.panel_pressure[aa, bb, kk, 7] + 1j*self.panel_pressure[aa, bb, kk, 8]) * self.panel_pressure[aa, bb, kk, 11] * self.panel_disc[2, int(self.panel_pressure[aa, bb, kk, 0]) - 1]
                            self.panel_pressure[aa, bb, kk, 13] = (self.panel_pressure[aa, bb, kk, 7] + 1j*self.panel_pressure[aa, bb, kk, 8]) * self.panel_pressure[aa, bb, kk, 11] * self.panel_disc[3, int(self.panel_pressure[aa, bb, kk, 0]) - 1]
                            self.panel_pressure[aa, bb, kk, 14] = (self.panel_pressure[aa, bb, kk, 7] + 1j*self.panel_pressure[aa, bb, kk, 8]) * self.panel_pressure[aa, bb, kk, 11] * self.panel_disc[4, int(self.panel_pressure[aa, bb, kk, 0]) - 1]
                            kk += 1

                found = 0

                while found == 0:
                    k = hydroD_results[ii].find(
                        '    BASIC PART  CONT.                                                                                                              \n')
                    ii += 1
                    if k != -1:
                        for jj in np.arange(0, 41, 1):
                            if hydroD_results[ii + 5 + jj] == '1\n':
                                break
                            self.panel_pressure[aa, bb, kk, 0:11] =  np.float_(hydroD_results[ii + 5 + jj].split())
                            self.panel_pressure[aa, bb, kk, 11] = self.panel_disc[1, int(self.panel_pressure[aa, bb, kk, 0]) - 1]
                            # [real, imag, abs]* area * normal direction
                            self.panel_pressure[aa, bb, kk, 12] = (self.panel_pressure[aa, bb, kk, 7] + 1j*self.panel_pressure[aa, bb, kk, 8]) * self.panel_pressure[aa, bb, kk, 11] * self.panel_disc[2, int(self.panel_pressure[aa, bb, kk, 0]) - 1]
                            self.panel_pressure[aa, bb, kk, 13] = (self.panel_pressure[aa, bb, kk, 7] + 1j*self.panel_pressure[aa, bb, kk, 8]) * self.panel_pressure[aa, bb, kk, 11] * self.panel_disc[3, int(self.panel_pressure[aa, bb, kk, 0]) - 1]
                            self.panel_pressure[aa, bb, kk, 14] = (self.panel_pressure[aa, bb, kk, 7] + 1j*self.panel_pressure[aa, bb, kk, 8]) * self.panel_pressure[aa, bb, kk, 11] * self.panel_disc[4, int(self.panel_pressure[aa, bb, kk, 0]) - 1]
                            kk += 1
                            if int(hydroD_results[ii + 5 + jj].split()[0]) == self.num_wet_panels:
                                found = 1
                                break

                found = 0
                if int(hydroD_results[ii + 5 + jj].split()[0]) == self.num_wet_panels and reflection == 0:
                    found = 1

                while found == 0:
                    k = hydroD_results[ii].find('REFLECTION')
                    ii += 1
                    if ii == len(hydroD_results) - 1:
                        found = 1
                    if k != -1:
                        for jj in np.arange(0, 60, 1):
                            if hydroD_results[ii + 5 + jj] == '1\n':
                                found = 1
                                break
                            self.panel_pressure[aa, bb, kk, 0:11] =  np.float_(hydroD_results[ii + 5 + jj].split())
                            self.panel_pressure[aa, bb, kk, 11] = self.panel_disc[1, int(self.panel_pressure[aa, bb, kk, 0]) - 1]
                            # [real, imag, abs]* area * normal direction
                            self.panel_pressure[aa, bb, kk, 12] = (self.panel_pressure[aa, bb, kk, 7] + 1j*self.panel_pressure[aa, bb, kk, 8]) * self.panel_pressure[aa, bb, kk, 11] * self.panel_disc[2, int(self.panel_pressure[aa, bb, kk, 0]) - 1]
                            self.panel_pressure[aa, bb, kk, 13] = (self.panel_pressure[aa, bb, kk, 7] + 1j*self.panel_pressure[aa, bb, kk, 8]) * self.panel_pressure[aa, bb, kk, 11] * -self.panel_disc[3, int(self.panel_pressure[aa, bb, kk, 0]) - 1]
                            self.panel_pressure[aa, bb, kk, 14] = (self.panel_pressure[aa, bb, kk, 7] + 1j*self.panel_pressure[aa, bb, kk, 8]) * self.panel_pressure[aa, bb, kk, 11] * self.panel_disc[4, int(self.panel_pressure[aa, bb, kk, 0]) - 1]
                            kk += 1

                found = 0
                if int(hydroD_results[ii + 5 + jj].split()[0]) == self.num_wet_panels and reflection == 0:
                    found = 1
                
                while found == 0:
                    k = hydroD_results[ii].find(
                        '    1. REFLECTION  CONT.                                                                                                           \n')
                    ii += 1
                    if k != -1:
                        for jj in np.arange(0, 60, 1):
                            if hydroD_results[ii + 5 + jj] == '1\n':
                                break
                            self.panel_pressure[aa, bb, kk, 0:11] =  np.float_(hydroD_results[ii + 5 + jj].split())
                            self.panel_pressure[aa, bb, kk, 11] = self.panel_disc[1, int(self.panel_pressure[aa, bb, kk, 0]) - 1]
                            # [real, imag, abs]* area * normal direction
                            self.panel_pressure[aa, bb, kk, 12] = (self.panel_pressure[aa, bb, kk, 7] + 1j*self.panel_pressure[aa, bb, kk, 8]) * self.panel_pressure[aa, bb, kk, 11] * self.panel_disc[2, int(self.panel_pressure[aa, bb, kk, 0]) - 1]
                            self.panel_pressure[aa, bb, kk, 13] = (self.panel_pressure[aa, bb, kk, 7] + 1j*self.panel_pressure[aa, bb, kk, 8]) * self.panel_pressure[aa, bb, kk, 11] * - self.panel_disc[3, int(self.panel_pressure[aa, bb, kk, 0]) - 1]
                            self.panel_pressure[aa, bb, kk, 14] = (self.panel_pressure[aa, bb, kk, 7] + 1j*self.panel_pressure[aa, bb, kk, 8]) * self.panel_pressure[aa, bb, kk, 11] * self.panel_disc[4, int(self.panel_pressure[aa, bb, kk, 0]) - 1]
                            kk += 1
                            if int(hydroD_results[ii + 5 + jj].split()[0]) == self.num_wet_panels:
                                found = 1
                                break

        self.panel_pressure[:, :, :, 7:10] *= (self.RO * self.G * self.WA)
        self.panel_pressure[:, :, :, 12:15] *= (self.RO * self.G * self.WA)

        self.front_column = self.panel_pressure[:, :, np.where(self.panel_pressure[0, 0, :, 4] < -5)[0], :].copy()
        self.left_column = self.panel_pressure[:, :, np.where(self.panel_pressure[0, 0, :, 5] < -15)[0], :].copy()
        self.right_column = self.panel_pressure[:, :, np.where(self.panel_pressure[0, 0, :, 5] > 15)[0], :].copy()

        self.front_column_force = np.zeros(shape=(2, 84, 3)) + 0j
        self.left_column_force = np.zeros(shape=(2, 84, 3)) + 0j
        self.right_column_force = np.zeros(shape=(2, 84, 3)) + 0j

        for ii in np.arange(0, self.numheadangles, 1).astype(int):
            for jj in np.arange(0, self.numwavelengths, 1).astype(int):
                self.front_column_force[ii, jj, :] = np.sum(self.front_column[ii, jj, :, 12:15], axis=0)

        for ii in np.arange(0, self.numheadangles, 1).astype(int):
            for jj in np.arange(0, self.numwavelengths, 1).astype(int):
                self.left_column_force[ii, jj, :] = np.sum(self.left_column[ii, jj, :, 12:15], axis=0)

        for ii in np.arange(0, self.numheadangles, 1).astype(int):
            for jj in np.arange(0, self.numwavelengths, 1).astype(int):
                self.right_column_force[ii, jj, :] = np.sum(self.right_column[ii, jj, :, 12:15], axis=0)

        if sectional_loads == 1:
            found = 0
            ii = 0
            while found == 0:
                k = hydroD_results[ii].find('            SECTIONAL LOADS ON THE PANEL MODEL ARE CALCULATED\n')
                ii += 1
                if k != -1:
                    self.numsectloads = np.int(hydroD_results[ii].split()[1])
                    found = 1

            self.stillwater_secloads = np.zeros(shape=(self.numsectloads, 6, 3))


            found = 0
            ii = 0
            jj = 0

            while found == 0:
                k = hydroD_results[ii].find('       STILLWATER SECTIONAL LOADS ON THE PANEL MODEL:\n')
                ii += 1
                if k != -1:
                    while jj < self.numsectloads:
                        k = 0
                        k = hydroD_results[ii].find('     ---------------------\n')
                        if k != -1:
                            for ll in np.linspace(0, 5, 6).astype(int):
                                self.stillwater_secloads[jj, ll, :] = hydroD_results[ii + 3 + ll].split()[1:4]
                            jj += 1
                        ii += 1
                    found = 1

            if fixed == 0:
                self.secloads = np.zeros(shape=(self.numheadangles, self.numwavelengths, self.numsectloads, 6, 1)) + 0j
                for nn in np.linspace(0, int(self.numwavelengths) - 1, int(self.numwavelengths)).astype(int):
                    for mm in np.linspace(0, int(self.numheadangles) - 1, int(self.numheadangles)).astype(int):
                        for jj in np.linspace(0, int(self.numsectloads) - 1, int(self.numsectloads)).astype(int):
                            found = 0
                            while found == 0:
                                ii += 1
                                k = hydroD_results[ii].find(
                                    '              REAL PART       IMAGINARY PART      REAL PART       IMAGINARY PART      ABSOLUTE VALUE  PHASE (DEGREES)\n')
                                if k != -1:
                                    for ll in np.linspace(0, 5, 6).astype(int):
                                        self.secloads[mm, nn, jj, ll, 0] = float(hydroD_results[ii + 2 + ll].split()[1]) + 1j* float(hydroD_results[ii + 2 + ll].split()[2])
                                    found = 1
            else:
                self.secloads = np.zeros(shape=(self.numheadangles, self.numwavelengths, self.numsectloads, 6, 1)) + 0j
                for nn in np.linspace(0, int(self.numwavelengths) - 1, int(self.numwavelengths)).astype(int):
                    for mm in np.linspace(0, int(self.numheadangles) - 1, int(self.numheadangles)).astype(int):
                        for jj in np.linspace(0, int(self.numsectloads) - 1, int(self.numsectloads)).astype(int):
                            found = 0
                            while found == 0:
                                ii += 1
                                k = hydroD_results[ii].find(
                                    '              REAL PART       IMAGINARY PART      ABSOLUTE VALUE  PHASE (DEGREES)    \n')
                                if k != -1:
                                    for ll in np.linspace(0, 5, 6).astype(int):
                                        self.secloads[mm, nn, jj, ll, :] = float(hydroD_results[ii + 2 + ll].split()[1]) + 1j* float(hydroD_results[ii + 2 + ll].split()[2])
                                    found = 1

            self.secloads[:, :, :, 0:3, :] *= self.RO*self.VOL*self.G*(self.WA/self.L)
            self.secloads[:, :, :, 3:6, :] *= self.RO * self.VOL * self.G * self.WA
            self.stillwater_secloads[:, 0:3, :] *= self.RO * self.VOL * self.G * (self.WA / self.L)
            self.stillwater_secloads[:, 3:6, :] *= self.RO * self.VOL * self.G * self.WA

