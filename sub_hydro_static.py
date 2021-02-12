import pandas as pd
import numpy as np


class FloaterParameters:
    def __init__(self, x, y, dia):
        # Free Variables
        self.x_space = x
        self.y_space = y

        # Currently fixed variables
        self.dia_column = dia
        self.dia_heave = 20
        self.heave_height = 2

        # Fixed variables
        self.hub_height = 85
        self.hub_space = 142.1
        self.draft = 14
        self.height = 25
        self.free_board = self.height - self.draft
        self.thickness = 0.07
        self.column_height = self.height - self.heave_height


class Density:
    def __init__(self):
        self.steel = 8500
        self.water = 1025
        self.concrete = 2500


class CrossSectionalArea:
    def __init__(self, floater):
        self.column = np.divide((np.pi * np.square(floater.dia_column)), 4)
        self.heave = np.divide((np.pi * np.square(floater.dia_heave)), 4)
        self.heave_top = self.heave - self.column

        self.cc_column = self.column - np.divide(
            (np.pi * np.square(floater.dia_column - 2 * floater.thickness)), 4)
        self.cc_heave = self.heave - np.divide(
            (np.pi * np.square(floater.dia_heave - np.multiply(2, floater.thickness))), 4)


class Buoy:
    def __init__(self, floater, area, rho):
        self.heave = area.heave * floater.heave_height * rho.water
        self.column = area.column * (floater.draft - floater.heave_height) * rho.water

        # Displaced Mass
        self.total = 3 * (self.heave + self.column)
        self.displaced_volume = self.total / rho.water


class Mass:
    def __init__(self, floater, area, buoy, rho):
        # System mass
        self.hub = 56780
        self.nacelle = 240000
        self.rotor = 110000

        t_bot = 0.04
        t_top = 0.02
        r_bot_outter = 3
        r_top_outter = r_bot_outter - t_bot
        r_out = 1.935
        r_in = r_out - t_top
        H = 87.6

        self.tower = 347460

        tower_volume_out = np.divide(np.pi, 3) * H * (np.square(r_bot_outter) + np.square(r_out) + r_bot_outter * r_out)
        tower_volume_in = np.divide(np.pi, 3) * H * (np.square(r_top_outter) + np.square(r_in) + r_top_outter * r_in)
        self.tower = (tower_volume_out - tower_volume_in) * rho.steel

        # Column mass
        self.column = area.cc_column * floater.column_height * rho.steel + area.column * floater.thickness * rho.steel

        # Heave Plate Mass
        self.heave = (area.cc_heave * floater.heave_height * rho.steel) + (
                area.heave_top * floater.thickness * rho.steel) + (area.heave * floater.thickness * rho.steel)
        float = 3 * (self.column + self.heave)
        self.turbines = 2 * (self.hub + self.nacelle + self.rotor + self.tower)
        back = 2 * (self.column + self.heave) + self.turbines
        front = (self.column + self.heave)
        total = self.turbines + float

        # Ballast weight
        self.ballast_total = buoy.total - total
        self.floater = float + self.ballast_total
        self.total = self.turbines + self.floater


class GeneralisedCoordinateSystem:
    def __init__(self, floater, area, mass, rho, buoy):
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
                       floater.hub_height / 1.3, mass.tower]
        self.towerR = ["Right Tower", -floater.x_space, np.divide(floater.y_space + floater.hub_space, 2) / 2,
                       floater.hub_height / 1.3, mass.tower]

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
        self.rog_df['I'] = self.rog_df['Mass [kg]'] * np.square(
            np.sqrt(np.square(self.rog_df['x']) + np.square(self.rog_df['y']) + np.square(self.rog_df['z'])))
        self.RoG = np.sqrt(np.sum(self.rog_df['I']) / np.sum(self.mass_df['Mass [kg]']))


def _ballast_allocation(coord, mass, floater):
    ratio = mass.ballast_total / np.sum(coord.mass_df['Mass [kg]'])
    COM_length = (coord.COB[0] - coord.COM[0])
    COM_ballast = coord.COB[0] + (COM_length / ratio)
    mass_ballast_back = ((-COM_ballast * mass.ballast_total) / floater.x_space) / 2
    mass_ballast_front = mass.ballast_total - 2 * mass_ballast_back
    return mass_ballast_front, mass_ballast_back


if __name__ == '__main__':
    mufp = FloaterParameters(150, 50)
    rho = Density()
    csa = CrossSectionalArea(mufp)
    buoy = Buoy(mufp, csa, rho)
    mass = Mass(mufp, csa, buoy, rho)
    coord = GeneralisedCoordinateSystem(mufp, csa, mass, rho, buoy)
