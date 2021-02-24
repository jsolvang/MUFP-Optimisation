import pandas as pd
import numpy as np


class MatrixCalculation:
    def __init__(self, coord, mass, floater, rho, env, area, buoy):
        # Initialising Mass and Stiffness Matrix
        self.mass = np.zeros(shape=(6, 6))
        self.stiffness = np.zeros(shape=(6, 6))

        # Mass Matrix
        self.mass[0, 0] = mass.total
        self.mass[1, 1] = mass.total
        self.mass[2, 2] = mass.total
        self.mass[3, 3] = mass.total * np.square(coord.RoG[0])
        self.mass[4, 4] = mass.total * np.square(coord.RoG[1])
        self.mass[5, 5] = mass.total * np.square(coord.RoG[2])

        # Stiffness Matrix
        front_front = abs(coord.COB[0]) + floater.dia_column / 2
        front_back = abs(coord.COB[0]) - floater.dia_column / 2
        back_front = abs(coord.columnbackR[1] - coord.COB[0]) - floater.dia_column / 2
        back_back = abs(coord.columnbackR[1] - coord.COB[0]) + floater.dia_column / 2

        left_left = abs(coord.columnbackR[2]) + floater.dia_column / 2
        left_right = abs(coord.columnbackR[2]) - floater.dia_column / 2
        right_left = abs(coord.columnbackR[2]) - floater.dia_column / 2
        right_right = abs(coord.columnbackR[2]) + floater.dia_column / 2

        front_centre = abs(coord.COB[0])
        back_centre = abs(coord.columnbackR[1] - coord.COB[0])
        left_centre = abs(coord.columnbackL[2])
        right_centre = abs(coord.columnbackR[2])
        buoy_mass = coord.COM[1] - coord.COB[1]

        # Recording output for trouble shooting
        self.distances = [["Front Front", front_front],
                          ["Front Back", front_back],
                          ["Back Front", back_front],
                          ["Back Back", back_back],
                          ["Left Left", left_left],
                          ["Left Right", left_right],
                          ["Right Left", right_left],
                          ["Right Right", right_right],
                          ["Front Centre", front_centre],
                          ["Back Centre", back_centre],
                          ["Left Centre", left_centre],
                          ["Right Centre", right_centre],
                          ["Buoy Mass", buoy_mass]]

        # Calculating water plane moment of inertia using the parallel axis theorem
        I44_1 = (np.pi / 64) * np.power(floater.dia_column, 4) + np.divide(np.pi * np.square(floater.dia_column),
                                                                           4) * np.square(left_centre)
        I44_2 = ((np.pi / 64) * np.power(floater.dia_column, 4) + np.divide(np.pi * np.square(floater.dia_column),
                                                                            4) * np.square(right_centre))
        I44_3 = (np.pi / 64) * np.power(floater.dia_column, 4)

        self.I44 = I44_1 + I44_2 + I44_3

        I55_1 = (np.pi / 64) * np.power(floater.dia_column, 4) + np.divide(np.pi * np.square(floater.dia_column),
                                                                           4) * np.square(front_centre)
        I55_2 = 2 * ((np.pi / 64) * np.power(floater.dia_column, 4) + np.divide(np.pi * np.square(floater.dia_column),
                                                                                4) * np.square(back_centre))

        self.I55 = I55_1 + I55_2

        # Calculating metacentric height and adjusting for centre of buoyancy
        BM_44 = self.I44 / buoy.displaced_volume
        BM_55 = self.I55 / buoy.displaced_volume

        KM_44 = BM_44 + coord.COB[2]
        KM_55 = BM_55 + coord.COB[2]

        # Calculating stiffness terms
        self.stiffness[2, 2] = rho.water * 3 * area.column * env.g
        self.stiffness[3, 3] = rho.water * env.g * buoy.displaced_volume * KM_44
        self.stiffness[4, 4] = rho.water * env.g * buoy.displaced_volume * KM_55