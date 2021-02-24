import pandas as pd
import numpy as np


class GeneralisedCoordinateSystem:
    def __init__(self, floater, area, mass, rho, buoy, env):
        self.columnfront = ["Front Column", 0, 0, (-floater.draft + floater.heave_height) / 2, mass.column]
        self.columnbackR = ["Back Right Column", -floater.x_space, floater.y_space / 2,
                            -(floater.draft - floater.heave_height) / 2, mass.column]
        self.columnbackL = ["Back Left Column", -floater.x_space, -floater.y_space / 2,
                            -(floater.draft - floater.heave_height) / 2, mass.column]

        self.heavefront = ["Front Heave Plate", 0, 0, (-floater.draft + floater.heave_height / 2), mass.heave]
        self.heavebackL = ["Back Left Heave Plate", -floater.x_space, -floater.y_space / 2,
                           (-floater.draft + floater.heave_height / 2), mass.heave]
        self.heavebackR = ["Back Right Heave Plate", -floater.x_space, floater.y_space / 2,
                           (-floater.draft + floater.heave_height / 2), mass.heave]

        self.hubL = ["Left RNA", -floater.x_space, -floater.hub_space / 2, floater.hub_height,
                     mass.hub + mass.rotor + mass.nacelle]
        self.hubR = ["Right RNA", -floater.x_space, floater.hub_space / 2, floater.hub_height,
                     mass.hub + mass.rotor + mass.nacelle]

        self.towerL = ["Left Tower", -floater.x_space, -np.divide(floater.y_space + floater.hub_space, 2) / 2,
                       floater.hub_height / 2.3, mass.tower]
        self.towerR = ["Right Tower", -floater.x_space, np.divide(floater.y_space + floater.hub_space, 2) / 2,
                       floater.hub_height / 2.3, mass.tower]

        self.buoy_columnfront = ["Front Column Buoy", 0, 0,
                                 - 1 / 2 * (floater.draft - floater.heave_height), buoy.column]
        self.buoy_columnbackR = ["Back Right Column Buoy", -floater.x_space, -floater.y_space / 2,
                                 - 1 / 2 * (floater.draft - floater.heave_height), buoy.column]
        self.buoy_columnbackL = ["Back Left Column Buoy", -floater.x_space, floater.y_space / 2,
                                 - 1 / 2 * (floater.draft - floater.heave_height), buoy.column]

        self.buoy_heavefront = ["Front Heave Buoy", 0, 0, (-floater.draft + floater.heave_height / 2), buoy.heave]
        self.buoy_heavebackR = ["Back Right Heave Buoy", -floater.x_space, -floater.y_space / 2,
                                (-floater.draft + floater.heave_height / 2), buoy.heave]
        self.buoy_heavebackL = ["Back Left Heave Buoy", -floater.x_space, floater.y_space / 2,
                                (-floater.draft + floater.heave_height / 2), buoy.heave]

        self.mass_df = pd.DataFrame([self.columnfront, self.columnbackL, self.columnbackR,
                                     self.heavefront, self.heavebackL, self.heavebackR,
                                     self.hubL, self.hubR,
                                     self.towerL, self.towerR])

        self.buoy_df = pd.DataFrame([self.buoy_columnfront, self.buoy_columnbackL, self.buoy_columnbackR,
                                     self.buoy_heavefront, self.buoy_heavebackL, self.buoy_heavebackR])
        self.mass_df.columns = ['Component', 'x', 'y', 'z', 'Mass [kg]']
        self.buoy_df.columns = ['Component', 'x', 'y', 'z', 'Buoy [kg]']

        self.mass_df.set_index('Component', inplace=True)
        self.buoy_df.set_index('Component', inplace=True)

        self.mass_df['weight_contribution'] = np.divide(self.mass_df['Mass [kg]'], np.sum(self.mass_df['Mass [kg]']))
        x = pd.DataFrame(self.mass_df['x'] * self.mass_df['weight_contribution'])
        x.columns = ['X_average']
        x['Y_average'] = pd.DataFrame(self.mass_df['y'] * self.mass_df['weight_contribution'])
        x['Z_average'] = pd.DataFrame(self.mass_df['z'] * self.mass_df['weight_contribution'])

        self.COM = [x['X_average'].sum(), x['Y_average'].sum(), x['Z_average'].sum()]
        self.buoy_df['buoy_contribution'] = np.divide(self.buoy_df['Buoy [kg]'], np.sum(self.buoy_df['Buoy [kg]']))
        x_buoy = pd.DataFrame(self.buoy_df['x'] * self.buoy_df['buoy_contribution'])
        x_buoy.columns = ['X_average']
        x_buoy['Y_average'] = self.buoy_df['y'] * self.buoy_df['buoy_contribution']
        x_buoy['Z_average'] = self.buoy_df['z'] * self.buoy_df['buoy_contribution']
        self.COB = [x_buoy['X_average'].sum(), x_buoy['Y_average'].sum(), x_buoy['Z_average'].sum()]

        [mass.ballast_front, mass.ballast_back] = _ballast_allocation(self, mass, floater)

        ballast_heightF = mass.ballast_front / (area.heave * rho.concrete)
        ballast_heightB = mass.ballast_back / (area.heave * rho.concrete)

        self.ballastfront = ["Front Ballast", 0, 0, -floater.draft + (ballast_heightF / 2), mass.ballast_front]
        self.ballastbackL = ["Back Left Ballast", -floater.x_space, -floater.y_space / 2,
                             -floater.draft + (ballast_heightB / 2), mass.ballast_back]
        self.ballastbackR = ["Back Right Ballast", -floater.x_space, floater.y_space / 2,
                             -floater.draft + (ballast_heightB / 2), mass.ballast_back]

        if ballast_heightF > floater.heave_height:
            excess_ballast_front = mass.ballast_front - area.heave * rho.concrete
            ballast_height_in_column = np.divide(excess_ballast_front, (area.column * rho.concrete))
            z_coord = np.divide(
                ((mass.ballast_front - excess_ballast_front) * (floater.heave_height / 2) + excess_ballast_front) * (
                    np.divide(ballast_height_in_column, 2)), mass.ballast_back)
            self.ballastfront = ["Front Ballast", 0, 0, -floater.draft + z_coord, mass.ballast_front]

        if ballast_heightB > floater.heave_height:
            excess_ballast_back = mass.ballast_back - area.heave * rho.concrete
            ballast_height_in_column = np.divide(excess_ballast_back, (area.column * rho.concrete))
            z_coord = np.divide(
                ((mass.ballast_back - excess_ballast_back) * (floater.heave_height / 2) + excess_ballast_back) * (
                    np.divide(ballast_height_in_column, 2)), mass.ballast_back)
            self.ballastbackL = ["Back Left Ballast", -floater.x_space, -floater.y_space / 2, -floater.draft + z_coord,
                                 mass.ballast_back]
            self.ballastbackR = ["Back Right Ballast", -floater.x_space, floater.y_space / 2, -floater.draft + z_coord,
                                 mass.ballast_back]

        self.mass_df = pd.DataFrame(
            [self.ballastfront, self.ballastbackL, self.ballastbackR, self.columnfront, self.columnbackL,
             self.columnbackR,
             self.heavefront, self.heavebackL, self.heavebackR,
             self.hubL, self.hubR,
             self.towerL, self.towerR])

        self.buoy_df = pd.DataFrame([self.buoy_columnfront, self.buoy_columnbackL, self.buoy_columnbackR,
                                     self.buoy_heavefront, self.buoy_heavebackL, self.buoy_heavebackR])
        self.mass_df.columns = ['Component', 'x', 'y', 'z', 'Mass [kg]']
        self.buoy_df.columns = ['Component', 'x', 'y', 'z', 'Buoy [kg]']

        self.mass_df.set_index('Component', inplace=True)
        self.buoy_df.set_index('Component', inplace=True)

        self.mass_df['weight_contribution'] = np.divide(self.mass_df['Mass [kg]'], np.sum(self.mass_df['Mass [kg]']))
        x = pd.DataFrame(self.mass_df['x'] * self.mass_df['weight_contribution'])
        x.columns = ['X_average']
        x['Y_average'] = pd.DataFrame(self.mass_df['y'] * self.mass_df['weight_contribution'])
        x['Z_average'] = pd.DataFrame(self.mass_df['z'] * self.mass_df['weight_contribution'])

        self.COM = [x['X_average'].sum(), x['Y_average'].sum(), x['Z_average'].sum()]
        self.buoy_df['buoy_contribution'] = np.divide(self.buoy_df['Buoy [kg]'], np.sum(self.buoy_df['Buoy [kg]']))
        x_buoy = pd.DataFrame(self.buoy_df['x'] * self.buoy_df['buoy_contribution'])
        x_buoy.columns = ['X_average']
        x_buoy['Y_average'] = self.buoy_df['y'] * self.buoy_df['buoy_contribution']
        x_buoy['Z_average'] = self.buoy_df['z'] * self.buoy_df['buoy_contribution']
        self.COB = [x_buoy['X_average'].sum(), x_buoy['Y_average'].sum(), x_buoy['Z_average'].sum()]

        # Radius of gyration
        self.rog_df = pd.DataFrame(self.mass_df['x'] - self.COM[0])
        self.rog_df['y'] = self.mass_df['y'] - self.COM[1]
        self.rog_df['z'] = self.mass_df['z'] - self.COM[2]
        self.rog_df['Mass [kg]'] = self.mass_df['Mass [kg]']
        self.rog_df['I_x'] = self.rog_df['Mass [kg]'] * np.square(
            np.sqrt(np.square(self.rog_df['x']) + np.square(self.rog_df['z'])))
        self.rog_df['I_y'] = self.rog_df['Mass [kg]'] * np.square(
            np.sqrt(np.square(self.rog_df['y']) + np.square(self.rog_df['z'])))
        self.rog_df['I_z'] = self.rog_df['Mass [kg]'] * np.square(
            np.sqrt(np.square(self.rog_df['x']) + np.square(self.rog_df['y'])))

        self.RoG = [np.sqrt(np.sum(self.rog_df['I_x']) / np.sum(self.mass_df['Mass [kg]'])),
                    np.sqrt(np.sum(self.rog_df['I_y']) / np.sum(self.mass_df['Mass [kg]'])),
                    np.sqrt(np.sum(self.rog_df['I_z']) / np.sum(self.mass_df['Mass [kg]']))]

        ### Mass Matrix

        floater.mass_matrix[0, 0] = mass.total
        floater.mass_matrix[1, 1] = mass.total
        floater.mass_matrix[2, 2] = mass.total
        floater.mass_matrix[3, 3] = mass.total * np.square(self.RoG[0])
        floater.mass_matrix[4, 4] = mass.total * np.square(self.RoG[1])
        floater.mass_matrix[5, 5] = mass.total * np.square(self.RoG[2])

        ### Stiffness Matrix

        front_front = abs(self.COB[0]) + floater.dia_column / 2
        front_back = abs(self.COB[0]) - floater.dia_column / 2
        back_front = abs(self.columnbackR[1] - self.COB[0]) - floater.dia_column / 2
        back_back = abs(self.columnbackR[1] - self.COB[0]) + floater.dia_column / 2

        left_left = abs(self.columnbackR[2]) + floater.dia_column / 2
        left_right = abs(self.columnbackR[2]) - floater.dia_column / 2
        right_left = abs(self.columnbackR[2]) - floater.dia_column / 2
        right_right = abs(self.columnbackR[2]) + floater.dia_column / 2

        front_centre = abs(self.COB[0])
        back_centre = abs(self.columnbackR[1] - self.COB[0])
        left_centre = abs(self.columnbackL[2])
        right_centre = abs(self.columnbackR[2])
        buoy_mass = self.COM[1] - self.COB[1]

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

        # Displaced volume pr unit angle, Assuming small angles
        volumedisp_front = area.column * front_back + 0.5 * (
                front_front - front_back) * area.column
        volumedisp_back = 2 * (area.column * back_front + 0.5 * (
                back_front - back_back) * area.column)
        volumedisp_left = area.column * left_right + 0.5 * (
                left_left - left_right) * area.column
        volumedisp_right = area.column * right_left + 0.5 * (
                right_right - right_left) * area.column

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

        BM_44 = self.I44 / buoy.displaced_volume
        BM_55 = self.I55 / buoy.displaced_volume

        KM_44 = BM_44 + self.COB[2]
        KM_55 = BM_55 + self.COB[2]

        floater.stiffness_matrix[2, 2] = rho.water * 3 * area.column * env.g

        floater.stiffness_matrix[3, 3] = rho.water*env.g*buoy.displaced_volume*KM_44

        floater.stiffness_matrix[4, 4] = rho.water*env.g*buoy.displaced_volume*KM_55

        def _ballast_allocation(coord, mass, floater):
            ratio = mass.ballast_total / np.sum(coord.mass_df['Mass [kg]'])
            COM_length = (coord.COB[0] - coord.COM[0])
            COM_ballast = coord.COB[0] + (COM_length / ratio)
            mass_ballast_back = ((-COM_ballast * mass.ballast_total) / floater.x_space) / 2
            mass_ballast_front = mass.ballast_total - 2 * mass_ballast_back
            return mass_ballast_front, mass_ballast_back