import pandas as pd
import numpy as np


class GlobalCoordinateSystem:
    # This class allocates the floaters mass and buoyancy into a global coordinate system (GCS).
    # The mass and buoyancy distribution are then used to calculate the moments of inertia and products of inertia
    # The GSC is centered around the front column center in the x/y axis and at the SWL in the z axis.
    def __init__(self, floater, area, mass, rho, buoy, env):
        self._determine_mass_coordinates(floater, area, mass, rho, buoy)
        self._determine_buoy_coordinates(floater, area, mass, rho, buoy)
        self._allocate_necessary_ballasting(floater, area, mass, rho)
        self._calculate_hydrod_inputs(floater, area, mass, rho, buoy, env)
        self._calculate_individual_column_mass_df()

    def _determine_mass_coordinates(self, floater, area, mass, rho, buoy):
        # Creating DataFrame of system mass coordinates using input variables in FloaterParameter object
        self.columnfront = ["Front Column", 0, 0, (-floater.draft + floater.heave_height) + floater.column_height / 2,
                            mass.column]
        self.columnbackR = ["Back Right Column", -floater.x_space, floater.y_space / 2,
                            (-floater.draft + floater.heave_height) + floater.column_height / 2, mass.column]
        self.columnbackL = ["Back Left Column", -floater.x_space, -floater.y_space / 2,
                            (-floater.draft + floater.heave_height) + floater.column_height / 2, mass.column]

        self.heavefront = ["Front Heave Plate", 0, 0, (-floater.draft + floater.heave_height / 2), mass.heave]
        self.heavebackL = ["Back Left Heave Plate", -floater.x_space, -floater.y_space / 2,
                           (-floater.draft + floater.heave_height / 2), mass.heave]
        self.heavebackR = ["Back Right Heave Plate", -floater.x_space, floater.y_space / 2,
                           (-floater.draft + floater.heave_height / 2), mass.heave]

        self.hubL = ["Left RNA", -floater.x_space, -floater.hub_space / 2, floater.hub_height,
                     mass.hub + mass.rotor + mass.nacelle]
        self.hubR = ["Right RNA", -floater.x_space, floater.hub_space / 2, floater.hub_height,
                     mass.hub + mass.rotor + mass.nacelle]

        self.towerL = ["Left Tower", -floater.x_space,
                       -floater.y_space / 2 - (np.subtract(floater.hub_space / 2, floater.y_space / 2) / 2.3),
                       floater.hub_height / 2.3, mass.tower]
        self.towerR = ["Right Tower", -floater.x_space,
                       floater.y_space / 2 + (np.subtract(floater.hub_space / 2, floater.y_space / 2) / 2.3),
                       floater.hub_height / 2.3, mass.tower]

        self.bracing_xyL = ["Left xy Bracing", -floater.x_space/2, -floater.y_space/4, 0, mass.xy_bracing]

        self.bracing_xyR = ["Right xy Bracing", -floater.x_space/2, floater.y_space/4, 0, mass.xy_bracing]

        self.bracing_y = ["Back xy Bracing", -floater.x_space, 0, 0, mass.y_bracing]

        self.unballasted_mass_df = pd.DataFrame([self.columnfront, self.columnbackL, self.columnbackR,
                                                 self.heavefront, self.heavebackL, self.heavebackR,
                                                 self.hubL, self.hubR,
                                                 self.towerL, self.towerR,
                                                 self.bracing_xyL, self.bracing_xyR, self.bracing_y])

        self.unballasted_mass_df.columns = ['Component', 'x', 'y', 'z', 'Mass [kg]']
        self.unballasted_mass_df.set_index('Component', inplace=True)

        # Calculating the un-ballasted COM by taking a weighted average of all mass components
        self.unballasted_mass_df['weight_contribution'] = np.divide(self.unballasted_mass_df['Mass [kg]'],
                                                                    np.sum(self.unballasted_mass_df['Mass [kg]']))
        x = pd.DataFrame(self.unballasted_mass_df['x'] * self.unballasted_mass_df['weight_contribution'])
        x.columns = ['X_average']
        x['Y_average'] = pd.DataFrame(self.unballasted_mass_df['y'] * self.unballasted_mass_df['weight_contribution'])
        x['Z_average'] = pd.DataFrame(self.unballasted_mass_df['z'] * self.unballasted_mass_df['weight_contribution'])
        self.COM_unballasted = [x['X_average'].sum(), x['Y_average'].sum(), x['Z_average'].sum()]

    def _determine_buoy_coordinates(self, floater, area, mass, rho, buoy):
        # Calculating the centre of buoyancy for each column and heave plate
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

        self.buoy_df = pd.DataFrame([self.buoy_columnfront, self.buoy_columnbackL, self.buoy_columnbackR,
                                     self.buoy_heavefront, self.buoy_heavebackL, self.buoy_heavebackR])

        self.buoy_df.columns = ['Component', 'x', 'y', 'z', 'Buoy [kg]']
        self.buoy_df.set_index('Component', inplace=True)

        # Calculating the COB by taking a weighted average of all buoyant components
        self.buoy_df['buoy_contribution'] = np.divide(self.buoy_df['Buoy [kg]'], np.sum(self.buoy_df['Buoy [kg]']))
        x_buoy = pd.DataFrame(self.buoy_df['x'] * self.buoy_df['buoy_contribution'])
        x_buoy.columns = ['X_average']
        x_buoy['Y_average'] = self.buoy_df['y'] * self.buoy_df['buoy_contribution']
        x_buoy['Z_average'] = self.buoy_df['z'] * self.buoy_df['buoy_contribution']
        self.COB = [x_buoy['X_average'].sum(), x_buoy['Y_average'].sum(), x_buoy['Z_average'].sum()]

    def _allocate_necessary_ballasting(self, floater, area, mass, rho):
        # Finding the required COM of the ballast to align final vertical centre of gravity (VCG) with the VCB
        ratio = mass.ballast_total / np.sum(self.unballasted_mass_df['Mass [kg]'])
        COM_length = (self.COB[0] - self.COM_unballasted[0])
        COM_ballast = self.COB[0] + (COM_length / ratio)

        mass.ballast_back = ((-COM_ballast * mass.ballast_total) / floater.x_space) / 2
        mass.ballast_front = mass.ballast_total - 2 * mass.ballast_back

        if mass.ballast_back >= 0:
            # Calculating the required height in heave plate for front and back ballasting
            ballast_heightF = mass.ballast_front / (area.heave * rho.concrete)
            ballast_heightB = mass.ballast_back / (area.heave * rho.concrete)
            self.ballast_heightF = ballast_heightF
            self.ballast_heightB = ballast_heightB

            # Assigning ballast weight and COM
            if ballast_heightF > floater.heave_height:
                excess_ballast_front = mass.ballast_front - area.heave * floater.heave_height * rho.concrete
                ballast_height_in_column = np.divide(excess_ballast_front, (area.column * rho.concrete))
                z_coord = np.divide(
                    ((mass.ballast_front - excess_ballast_front) * (floater.heave_height / 2) + excess_ballast_front) * (
                        np.divide(ballast_height_in_column, 2)), mass.ballast_back)
                self.ballastfront = ["Front Ballast", 0, 0, -floater.draft + z_coord, mass.ballast_front]
            else:
                self.ballastfront = ["Front Ballast", 0, 0, -floater.draft + (ballast_heightF / 2), mass.ballast_front]

            if ballast_heightB > floater.heave_height:
                excess_ballast_back = mass.ballast_back - area.heave * floater.heave_height * rho.concrete
                ballast_height_in_column = np.divide(excess_ballast_back, (area.column * rho.concrete))
                z_coord = np.divide(
                    ((mass.ballast_back - excess_ballast_back) * (floater.heave_height / 2) + excess_ballast_back) * (
                        np.divide(ballast_height_in_column, 2)), mass.ballast_back)
                self.ballastbackL = ["Back Left Ballast", -floater.x_space, -floater.y_space / 2, -floater.draft + z_coord,
                                     mass.ballast_back]
                self.ballastbackR = ["Back Right Ballast", -floater.x_space, floater.y_space / 2, -floater.draft + z_coord,
                                     mass.ballast_back]
            else:
                self.ballastbackL = ["Back Left Ballast", -floater.x_space, -floater.y_space / 2,
                                     -floater.draft + (ballast_heightB / 2), mass.ballast_back]
                self.ballastbackR = ["Back Right Ballast", -floater.x_space, floater.y_space / 2,
                                     -floater.draft + (ballast_heightB / 2), mass.ballast_back]
        else:
            mass.ballast_back = 0
            mass.ballast_front = mass.ballast_total

            ballast_heightF = mass.ballast_front / (area.heave * rho.concrete)
            ballast_heightB = mass.ballast_back / (area.heave * rho.concrete)
            self.ballast_heightF = ballast_heightF
            self.ballast_heightB = ballast_heightB

            if ballast_heightF > floater.heave_height:
                excess_ballast_front = mass.ballast_front - area.heave * floater.heave_height * rho.concrete
                ballast_height_in_column = np.divide(excess_ballast_front, (area.column * rho.concrete))
                z_coord = np.divide(
                    ((mass.ballast_front - excess_ballast_front) * (floater.heave_height / 2) + excess_ballast_front) * (
                        np.divide(ballast_height_in_column, 2)), mass.ballast_back)
                self.ballastfront = ["Front Ballast", 0, 0, -floater.draft + z_coord, mass.ballast_front]
            else:
                self.ballastfront = ["Front Ballast", 0, 0, -floater.draft + (ballast_heightF / 2), mass.ballast_front]

            self.ballastbackL = ["Back Left Ballast", -floater.x_space, -floater.y_space / 2,
                                 -floater.draft, mass.ballast_back]
            self.ballastbackR = ["Back Right Ballast", -floater.x_space, floater.y_space / 2,
                                 -floater.draft, mass.ballast_back]

        # Creating final mass DataFrame including ballasting
        self.mass_df = pd.DataFrame(
            [self.ballastfront, self.ballastbackL, self.ballastbackR, self.columnfront, self.columnbackL,
             self.columnbackR, self.heavefront, self.heavebackL, self.heavebackR, self.hubL, self.hubR, self.towerL,
             self.towerR, self.bracing_xyL, self.bracing_xyR, self.bracing_y])
        self.mass_df.columns = ['Component', 'x', 'y', 'z', 'Mass [kg]']
        self.mass_df.set_index('Component', inplace=True)
        self.mass_df['weight_contribution'] = np.divide(self.mass_df['Mass [kg]'], np.sum(self.mass_df['Mass [kg]']))


    def _calculate_hydrod_inputs(self, floater, area, mass, rho, buoy, env):

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
        self.rog_df['I_x'] = self.rog_df['Mass [kg]'] * (np.square(self.rog_df['x']) + np.square(self.rog_df['z']))
        self.rog_df['I_y'] = self.rog_df['Mass [kg]'] * (np.square(self.rog_df['y']) + np.square(self.rog_df['z']))
        self.rog_df['I_z'] = self.rog_df['Mass [kg]'] * (np.square(self.rog_df['x']) + np.square(self.rog_df['y']))

        self.RoG = [np.sqrt(np.sum(self.rog_df['I_x']) / np.sum(self.mass_df['Mass [kg]'])),
                    np.sqrt(np.sum(self.rog_df['I_y']) / np.sum(self.mass_df['Mass [kg]'])),
                    np.sqrt(np.sum(self.rog_df['I_z']) / np.sum(self.mass_df['Mass [kg]']))]

        self.PoI_df = pd.DataFrame(self.mass_df['x'] - self.COM[0])
        self.PoI_df['y'] = self.mass_df['y'] - self.COM[1]
        self.PoI_df['z'] = self.mass_df['z'] - self.COM[2]
        self.PoI_df['Mass [kg]'] = self.mass_df['Mass [kg]']
        self.PoI_df['Ixy'] = self.PoI_df['x'] * self.PoI_df['y'] * self.PoI_df['Mass [kg]']
        self.PoI_df['Iyz'] = self.PoI_df['y'] * self.PoI_df['z'] * self.PoI_df['Mass [kg]']
        self.PoI_df['Ixz'] = self.PoI_df['x'] * self.PoI_df['z'] * self.PoI_df['Mass [kg]']
        self.PoI = [sum(self.PoI_df['Ixy']) / mass.total, sum(self.PoI_df['Ixz']) / mass.total,
                    sum(self.PoI_df['Ixy']) / mass.total]

    def _calculate_individual_column_mass_df(self):
        # Creating mass dataframes for each column

        self.front_column_df = pd.DataFrame(
            [self.ballastfront, self.columnfront, self.heavefront])
        self.front_column_df.columns = ['Component', 'x', 'y', 'z', 'Mass [kg]']
        self.front_column_df.set_index('Component', inplace=True)
        self.front_column_df['weight_contribution'] = np.divide(self.front_column_df['Mass [kg]'], np.sum(self.front_column_df['Mass [kg]']))
        self.front_column_df['x'] += self.COM[0]
        self.front_column_df['y'] += self.COM[1]

        self.left_column_df = pd.DataFrame(
            [self.ballastbackL, self.columnbackL,self.heavebackL, self.hubL, self.towerL])
        self.left_column_df.columns = ['Component', 'x', 'y', 'z', 'Mass [kg]']
        self.left_column_df.set_index('Component', inplace=True)
        self.left_column_df['weight_contribution'] = np.divide(self.left_column_df['Mass [kg]'], np.sum(self.left_column_df['Mass [kg]']))
        self.left_column_df['x'] += self.COM[0]
        self.left_column_df['y'] += self.COM[1]

        self.right_column_df = pd.DataFrame(
            [self.ballastbackR, self.columnbackR, self.heavebackR, self.hubR, self.towerR])
        self.right_column_df.columns = ['Component', 'x', 'y', 'z', 'Mass [kg]']
        self.right_column_df.set_index('Component', inplace=True)
        self.right_column_df['weight_contribution'] = np.divide(self.right_column_df['Mass [kg]'], np.sum(self.right_column_df['Mass [kg]']))
        self.right_column_df['x'] += self.COM[0]
        self.right_column_df['y'] += self.COM[1]
