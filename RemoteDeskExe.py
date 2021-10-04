import numpy as np
import pickle
import sys
sys.path.append(r'C:\Users\Joar\Documents\1_Education\NTNU\OneDrive - NTNU\Thesis\Modelling\FD Model')
from FloaterParameters import FloaterParameters
from Environment import Environment
from Buoy import Buoy
from Mass import Mass
from Density import Density
from Area import Area
from GlobalCoordinateSystem import GlobalCoordinateSystem
from SystemMatrices import MatrixCalculation
from ComputeHydroCoefficients import CompHydroCoefficient
from plot_hydroD_results import plot_hydroD_results
from ReadWadamLis import ReadWadamLis

if __name__ == '__main__':
    x_space = np.arange(20,110,10)
    y_space = np.arange(20,110,10)
    column_diameter = [11,12,13,14]
    file_loc = r'C:\Users\Joar\Documents\1_Education\NTNU\pickle_files'
    for dd in column_diameter:
        for xx in x_space:
            for yy in y_space:
                file_name = "\sim_x_%d_y_%d_D%d" % (xx, yy, dd)
                file_path = file_loc + file_name
                mufp = FloaterParameters(xx, yy, dd)
                env = Environment()
                rho = Density()
                csa = Area(mufp)
                buoy = Buoy(mufp, csa, rho)
                mass = Mass(mufp, csa, buoy, rho)
                coord = GlobalCoordinateSystem(mufp, csa, mass, rho, buoy, env)
                hydrod_results = CompHydroCoefficient(coord, mass, mufp, run)
                sim = ReadWadamLis(hydrod_results.HydroD_result, 0, 1)
                matrix = MatrixCalculation(coord, mass, mufp, rho, env, csa, buoy)
                plot_hydroD_results(sim, mufp, env, 1)
                file_sim = open(file_path, 'wb')
                pickle.dump(sim, file_sim)
                file_sim.close()