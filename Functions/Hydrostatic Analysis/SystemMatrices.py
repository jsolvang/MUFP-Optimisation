import pandas as pd
import numpy as np


class MatrixCalculation:
    # Calculating mass matrix from system mass, moment of inertia and product of inertia
    # Calculating stiffness matrix from the metacentric height of the floater
    def __init__(self, coord, mass, floater, rho, env, area, buoy)
        self.mass = np.zeros(shape=(6, 6))
        self.stiffness = np.zeros(shape=(6, 6))

        # Mass Matrix
        self.mass[0, 0] = mass.total
        self.mass[1, 1] = mass.total
        self.mass[2, 2] = mass.total
        self.mass[3, 3] = mass.total * np.square(coord.RoG[0])
        self.mass[4, 4] = mass.total * np.square(coord.RoG[1])
        self.mass[5, 5] = mass.total * np.square(coord.RoG[2])
        self.mass[3, 5] = -mass.total * coord.PoI[1]
        self.mass[5, 3] = -mass.total * coord.PoI[1]
        self.mass[0, 4] = mass.total * coord.COM[2]
        self.mass[4, 0] = mass.total * coord.COM[2]
        self.mass[1, 3] = -mass.total * coord.COM[2]
        self.mass[3, 1] = -mass.total * coord.COM[2]

        # Stiffness Matrix
        # Calculating distances from COM and the nearest and furthest edges of each of the columns
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
        buoy_mass = coord.COM[2] - coord.COB[2]

        # Recording output for trouble shooting
        self.distances = {"Front Front": front_front,
                          "Front Back": front_back,
                          "Back Front": back_front,
                          "Back Back": back_back,
                          "Left Left": left_left,
                          "Left Right": left_right,
                          "Right Left": right_left,
                          "Right Right": right_right,
                          "Front Centre": front_centre,
                          "Back Centre": back_centre,
                          "Left Centre": left_centre,
                          "Right Centre": right_centre,
                          "CoB to CoM": buoy_mass}

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

        KM_44 = BM_44 + (coord.COB[2] - coord.COM[2])
        KM_55 = BM_55 + (coord.COB[2] - coord.COM[2])

        # Calculating stiffness terms
        self.stiffness[2, 2] = rho.water * 3 * area.column * env.g
        self.stiffness[3, 3] = rho.water * env.g * buoy.displaced_volume * KM_44
        self.stiffness[4, 4] = rho.water * env.g * buoy.displaced_volume * KM_55

        if self.stiffness[3, 3] < 5e6:
            self.stiffness[3, 3] = 1e6 - 5e5*((60 - floater.y_space)/30)
        if self.stiffness[4, 4] < 1e7:
            self.stiffness[4, 4] = 1e7 - 5e6*((60 - floater.x_space)/30)
