import openmdao.api as om

import numpy as np
import subprocess
import os
from scipy.interpolate import interp1d


class CompHydroCoefficient:

    def __init__(self):
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
        self.HydroD_JScommand = r'C:\Users\Joar\Documents\1_Education\NTNU\python_hyd_script_com.js'
        self.HydroD_w = 'HydroD_w'

    def compute(self, coord, mass, floater, run):
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
                if kk == 68:
                    file.write("PanelModel1.setModelTranslation(Vector3d(%f m,0 m,0 m)); \n" % coord.COM[0])
                elif kk == 80:
                    file.write("MassModel1.setTotalMass(%f Kg); \n" % mass.total)
                elif kk == 81:
                    file.write("MassModel1.setCOG(Point(0 m,0 m,%f m)); \n" % coord.COM[2])
                elif kk == 82:
                    file.write("MassModel1.setRadiusGyration(Vector3d(%f m,%f m,%f m)); \n" % (coord.RoG[0], coord.RoG[1], coord.RoG[2]))
                else:
                    file.write(lines[kk])

        if run == 1:
            # Run GeniE
            subprocess.run(
                self.GeniE_path + " " + self.GeniE_database + " /new " + self.GeniE_license + " /com=" + self.GeniE_JScommand + " /exit")

            os.remove(self.GeniE_database)

            # Run HydroD
            subprocess.run(
                 self.HydroD_path + " " + self.HydroD_database + " /new " + self.HydroD_license + " /com=" + self.HydroD_JScommand + " /exit")
            os.remove(self.HydroD_database)



