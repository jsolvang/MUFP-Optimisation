import openmdao.api as om

import numpy as np
import subprocess
import os
from scipy.interpolate import interp1d

class comp_hydrocoeff():

    def initialize(self):
        self.options.declare('GeniE_dict', types=dict)
        self.options.declare('HydroD_dict', types=dict)
        self.options.declare('freqs_dict', types=dict)

    def setup(self):
        # GeniE = self.options['GeniE_dict']
        self.GeniE_path = r"C:\Program Files\DNVGL\GeniE V7.11-07\Program\GenieR.exe"
        self.GeniE_database = r"C:\Users/Joar/Documents/1_Education/NTNU/test.gnx"
        self.GeniE_license = r"C:\Program Files\DNVGL\license.lic"
        self.GeniE_JScommand = r"C:\Users\Joar\Documents\1_Education\NTNU\3_column_semisub.js"

        # HydroD = self.options['HydroD_dict']
        self.HydroD_path = r'C:\Program Files (x86)\DNVGL\HydroD V4.10-01\Program\HydroD.exe'
        self.HydroD_result = r'C:\Users\Joar\Documents\1_Education\NTNU\WadamRun1'
        self.HydroD_database = r'C:\Users/Joar/Documents/1_Education/NTNU/test.hyd'
        self.HydroD_license = r'C:\Program Files\DNVGL\license.lic'
        self.HydroD_JScommand = r'C:\Users\Joar\Documents\1_Education\NTNU\Python_Hyd_script.js'
        self.HydroD_w = 'HydroD_w'


    def compute(self, coord, mass, floater):
        omega = [0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2, 2.1, 2.2, 2.25, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3, 3.2, 3.4, 3.6, 3.8, 4, 4.2, 4.6, 5]
        N_omega = len(omega)
        HydroD_w = self.HydroD_w
        N_HydroD_w = len(HydroD_w)


        # Write the new value of D and draft in the GeniE file
        with open(self.GeniE_JScommand, 'r') as file:
            lines = file.readlines()
        with open(self.GeniE_JScommand, 'w') as file:
            for kk in range(len(lines)):
                if kk == 10:
                    file.write("Diameter2 = %fm; // side column diameter\n" % floater.dia_column)
                elif kk == 19:
                    file.write("x2 = %fm; 	// center-center x distance for side column \n" % floater.x_space)
                elif kk == 20:
                    file.write("y2 = %fm; 	// center-center x distance for side column \n" % floater.y_space)
                else:
                    file.write(lines[kk])

        # Write the total mass and the CoG location in the HydroD file
        with open(self.HydroD_JScommand, 'r') as file:
            lines = file.readlines()
        with open(self.HydroD_JScommand, 'w') as file:
            for kk in range(len(lines)):
                if kk == 97:
                    file.write("MassModel1.setTotalMass(%f Kg); \n" % mass.total)
                elif kk == 98:
                    file.write("MassModel1.setCOG(Point(%f m,%f m,%f m)); \n" % (-coord.COM[0], coord.COM[1], coord.COM[2]))
                    print(coord.COM)
                elif kk == 99:
                    file.write("MassModel1.setRadiusGyration(Vector3d(%f m,%f m,%f m)); \n" % (coord.RoG, coord.RoG, coord.RoG))
                    print(coord.RoG)
                else:
                    file.write(lines[kk])

        # Run GeniE
        subprocess.run(
            self.GeniE_path + " " + self.GeniE_database + " /new " + self.GeniE_license + " /com=" + self.GeniE_JScommand + " /exit")
        os.remove(self.GeniE_database)

        # Run HydroD
        subprocess.run(
            self.HydroD_path + " " + self.HydroD_database + " /new " + self.HydroD_license + " /com=" + self.HydroD_JScommand + " /exit")
        os.remove(self.HydroD_database)

        def _results_retrieval(self):
            # Reads results
            folder1 = self.HydroD_result + "WADAM1.LIS"
            with open(self.HydroD_JScommand, 'r') as hydroD_results:
                lines = hydroD_results.readlines()



            A_tot, B_tot, Fexc_tot, Phexc_tot = read_WAMIT2(folder1, folder2, N_HydroD_w)

            # Interpolate the results in the more refined frequency vector
            A_HydroD = A_tot[np.ix_(np.arange(N_HydroD_w), [0, 4], [0, 4])]
            fA11 = interp1d(HydroD_w, A_HydroD[:, 0, 0])
            fA15 = interp1d(HydroD_w, A_HydroD[:, 0, 1])
            fA51 = interp1d(HydroD_w, A_HydroD[:, 1, 0])
            fA55 = interp1d(HydroD_w, A_HydroD[:, 1, 1])
            outputs['A'][:, 0, 0] = fA11(omega)
            outputs['A'][:, 0, 1] = fA15(omega)
            outputs['A'][:, 1, 0] = fA51(omega)
            outputs['A'][:, 1, 1] = fA55(omega)

            fA33 = interp1d(HydroD_w, A_tot[:, 2, 2])
            outputs['A33'] = fA33(omega)

            B_HydroD = B_tot[np.ix_(np.arange(N_HydroD_w), [0, 4], [0, 4])]
            fB11 = interp1d(HydroD_w, B_HydroD[:, 0, 0])
            fB15 = interp1d(HydroD_w, B_HydroD[:, 0, 1])
            fB51 = interp1d(HydroD_w, B_HydroD[:, 1, 0])
            fB55 = interp1d(HydroD_w, B_HydroD[:, 1, 1])
            outputs['B'][:, 0, 0] = fB11(omega)
            outputs['B'][:, 0, 1] = fB15(omega)
            outputs['B'][:, 1, 0] = fB51(omega)
            outputs['B'][:, 1, 1] = fB55(omega)

            Fexc_HydroD = Fexc_tot[np.ix_(np.arange(N_HydroD_w), [0, 4])]
            fF1 = interp1d(HydroD_w, Fexc_HydroD[:, 0])
            fF5 = interp1d(HydroD_w, Fexc_HydroD[:, 1])
            outputs['Fexc'][:, 0] = fF1(omega)
            outputs['Fexc'][:, 1] = fF5(omega)

            Phexc_HydroD = Phexc_tot[np.ix_(np.arange(N_HydroD_w), [0, 4])]
            fP1 = interp1d(HydroD_w, Phexc_HydroD[:, 0])
            fP5 = interp1d(HydroD_w, Phexc_HydroD[:, 1])
            outputs['Phexc'][:, 0] = fP1(omega)
            outputs['Phexc'][:, 1] = fP5(omega)

