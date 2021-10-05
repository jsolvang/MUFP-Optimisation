import pandas as pd
import numpy as np

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
    run = 1
    mufp = FloaterParameters(150, 50, 12)
    env = Environment()
    rho = Density()
    csa = Area(mufp)
    buoy = Buoy(mufp, csa, rho)
    mass = Mass(mufp, csa, buoy, rho)
    coord = GlobalCoordinateSystem(mufp, csa, mass, rho, buoy, env)
    hydrod_results = CompHydroCoefficient(coord, mass, mufp, run)
    sim = ReadWadamLis(hydrod_results.HydroD_result)
    matrix = MatrixCalculation(coord, mass, mufp, rho, env, csa, buoy)
    plot_hydroD_results(sim)
