import numpy as np
from scipy import interpolate
import pandas as pd


class InternalForces:

    def __init__(self, column_pressure_panels, force_column, mass_column_df, RAO, sim, f_rad, area, rho, env, buoy, mufp, coord, mass):
        # Creating empty matrices
        self.internal_force = np.zeros(shape=(sim.numheadangles, len(RAO[0, 0, :, 0]), 6)) + 0j
        self.mass_matrix, self.COM = calc_column_mass_matrix(mass_column_df,area, rho, mufp, mass, coord)
        self.mass_column = self.mass_matrix[0,0]
        self.stiffness_matrix = calc_column_stiffness(rho, area, env, buoy, mufp, coord, self.COM)

        # Interpolating column force data to match RAO sampling frequency
        force_column_interp = np.zeros(shape=(sim.numheadangles, len(RAO[0, 0, :, 0]), 3)) + 0j
        for ii, _ in enumerate(sim.front_column_force[:, 0, 0]):
            for jj, _ in enumerate(sim.front_column_force[0, 0, :]):
                force_column_interp[ii, :, jj] = np.interp(f_rad, sim.wave_disc[:, 4], force_column[ii, :, jj])

        self.force_column = force_column_interp

        # Calculating internal forces transfer functions
        self.translation_internal_forces(RAO, sim)
        self.rotational_internal_forces(column_pressure_panels, RAO, sim)

    def translation_internal_forces(self, RAO, sim):
        # Code returning the internal force applied to a column in each of the translational degrees of freedom
        for ii in np.arange(0, sim.numheadangles, 1).astype(int):
            # RAO indexing [headangle,kinemetics,frequency,dof]
            # Internal Force Surge
            self.internal_force[ii, :, 0] = RAO[ii, 2, :, 0] * self.mass_column - self.force_column[ii, :, 0] + RAO[ii, 0, :, 0]*self.stiffness_matrix[0, 0]
            # Internal Force Sway
            self.internal_force[ii, :, 1] = RAO[ii, 2, :, 1] * self.mass_column - self.force_column[ii, :, 1] + RAO[ii, 0, :, 1]*self.stiffness_matrix[1, 1]
            # Internal Force Heave
            self.internal_force[ii, :, 2] = RAO[ii, 2, :, 2] * self.mass_column - self.force_column[ii, :, 2] + RAO[ii, 0, :, 2]*self.stiffness_matrix[2, 2]

    def rotational_internal_forces(self, column_pressure_panels, RAO, sim):
        # Rotational Forces
        # Internal forces roll
        self.rot_force_column = np.zeros(shape=(2, len(RAO[0, 0, :, 0]), 3)) * 1j
        self.rot_force_column[:, :, 0] = rot_calc(column_pressure_panels, sim, 'roll', self.COM)
        self.rot_force_column[:, :, 1] = rot_calc(column_pressure_panels, sim, 'pitch', self.COM)
        self.rot_force_column[:, :, 2] = rot_calc(column_pressure_panels, sim, 'yaw',  self.COM)

        # Finding internal forces
        for ii in np.arange(0, sim.numheadangles, 1).astype(int):
            self.internal_force[ii, :, 3] = RAO[ii, 2, :, 3] * self.mass_matrix[3,3] - self.rot_force_column[ii, :, 0] + RAO[ii, 0, :, 3]*self.stiffness_matrix[3, 3]
            self.internal_force[ii, :, 4] = RAO[ii, 2, :, 4] * self.mass_matrix[4,4] - self.rot_force_column[ii, :, 1] + RAO[ii, 0, :, 4]*self.stiffness_matrix[4, 4]
            self.internal_force[ii, :, 5] = RAO[ii, 2, :, 5] * self.mass_matrix[5,5] - self.rot_force_column[ii, :, 2] + RAO[ii, 0, :, 5]*self.stiffness_matrix[5, 5]

def rot_calc(column_pressure_panels, sim, dof, COM):
    # Rotational Forces
    # Creating empty force matrix
    force_rao = np.zeros(shape=(2, len(column_pressure_panels[0, :, 0, 0]))) * 1j
    if dof == 'roll':
        force1 = column_pressure_panels[:, :, :, 13] * -(column_pressure_panels[:, :, :, 6])
        force2 = column_pressure_panels[:, :, :, 14] * (column_pressure_panels[:, :, :, 5] - COM[1])

        for ii in np.arange(0, sim.numheadangles, 1).astype(int):
            for jj in np.arange(0, sim.numwavelengths, 1).astype(int):
                force_rao[ii, jj] = np.sum(force1[ii, jj, :] + force2[ii, jj, :])


    elif dof == 'pitch':
        force1 = column_pressure_panels[:, :, :, 12] * (column_pressure_panels[:, :, :, 6])
        force2 = column_pressure_panels[:, :, :, 14] * -(column_pressure_panels[:, :, :, 4] - COM[0])

        for ii in np.arange(0, sim.numheadangles, 1).astype(int):
            for jj in np.arange(0, sim.numwavelengths, 1).astype(int):
                force_rao[ii, jj] = np.sum(force1[ii, jj, :]) + np.sum(force2[ii, jj, :])

    elif dof == 'yaw':
        force1 = column_pressure_panels[:, :, :, 12] * -(column_pressure_panels[:, :, :, 5] - COM[1])
        force2 = column_pressure_panels[:, :, :, 13] * (column_pressure_panels[:, :, :, 4] - COM[0])

        for ii in np.arange(0, sim.numheadangles, 1).astype(int):
            for jj in np.arange(0, sim.numwavelengths, 1).astype(int):
                force_rao[ii, jj] = np.sum(force1[ii, jj, :]) + np.sum(force2[ii, jj, :])
    else:
        print('invalid degree of freedom')

    return force_rao


def calc_column_mass_matrix(mass_column_df, area, rho, mufp, mass, coord):
    mass_matrix = np.zeros(shape=(6, 6))
    mass_column = sum(mass_column_df['Mass [kg]'])
    COM = calc_centre_of_mass(mass_column_df)
    rog_df, RoG = calc_radius_of_gyration(mass_column_df, COM)
    PoI, PoI_df = calc_product_inertia(mass_column_df, COM, mass_column)

    # Mass Matrix
    mass_matrix[0, 0] = mass_column
    mass_matrix[1, 1] = mass_column
    mass_matrix[2, 2] = mass_column
    mass_matrix[3, 3] = mass_column * np.square(RoG[0]) + mass_column*COM[2]**2
    mass_matrix[4, 4] = mass_column * np.square(RoG[1]) + mass_column*COM[2]**2
    mass_matrix[3, 5] = -mass_column * np.square(PoI[1])
    mass_matrix[5, 3] = -mass_column * np.square(PoI[1])
    mass_matrix[0, 4] = mass_column * COM[2]
    mass_matrix[4, 0] = mass_column * COM[2]
    mass_matrix[1, 3] = -mass_column * COM[2]
    mass_matrix[3, 1] = -mass_column * COM[2]
    mass_matrix[5, 5] = 0.5 * area.column * mufp.thickness * rho.steel * (mufp.dia_column / 2) ** 2 \
                       + area.cc_column * mufp.column_height * rho.steel * (mufp.dia_column / 2) ** 2 \
                       + 0.5 * (area.heave_top * mufp.thickness * rho.steel) * ((mufp.dia_column / 2) ** 2 + (mufp.dia_heave / 2) ** 2) \
                       + area.cc_heave * mufp.heave_height * rho.steel * (mufp.dia_heave / 2) ** 2 \
                       + 0.5 * area.heave * mufp.thickness * rho.steel * (mufp.dia_heave / 2) ** 2 \
                       + 0.5 * mass.ballast_front * (mufp.dia_heave / 2) ** 2

    m1 = area.column * mufp.thickness * rho.steel
    m2 = area.cc_column * mufp.column_height * rho.steel
    m3 = area.heave_top * mufp.thickness * rho.steel
    m4 = area.cc_heave * mufp.heave_height * rho.steel
    m5 = mass.ballast_front
    m6 = area.heave * mufp.thickness * rho.steel

    I1 = 0.25 * m1 * (mufp.dia_column / 2) ** 2 + m1 * (mufp.free_board ** 2)
    I2 = (1 / 12) * m2 * (3 * ((mufp.dia_column / 2) ** 2 + ((mufp.dia_column / 2) - mufp.thickness) ** 2) + mufp.column_height ** 2) \
         + m2 * (mufp.free_board- mufp.column_height/2) ** 2
    I3 = 0.5 * m3 * ((mufp.dia_heave / 2) ** 2 + ((mufp.dia_column / 2) - mufp.thickness) ** 2) + m3 * 12 ** 2
    I4 = (1 / 12) * m4 * (3 * ((mufp.dia_heave / 2) ** 2 + ((mufp.dia_heave / 2) - mufp.thickness) ** 2) + mufp.heave_height ** 2) \
            + m4 * (mufp.draft -1) ** 2
    I5 = 0.25 * m5 * (mufp.dia_heave / 2) ** 2 + m5 * (mufp.draft - coord.ballast_heightF) ** 2
    I6 = 0.25 * m6 * (mufp.dia_heave / 2) ** 2 + m6 * (mufp.draft ** 2)

    mass_matrix[3, 3] = I1 + I2 + I3 + I4 + I5 + I6
    mass_matrix[4, 4] = I1 + I2 + I3 + I4 + I5 + I6

    return mass_matrix, COM

def calc_product_inertia(mass_column_df, COM, mass_column):
    PoI_df = pd.DataFrame(mass_column_df['x'] - COM[0])
    PoI_df['y'] = mass_column_df['y'] - COM[1]
    PoI_df['z'] = mass_column_df['z'] - COM[2]
    PoI_df['Mass [kg]'] = mass_column_df['Mass [kg]']
    PoI_df['Ixy'] = PoI_df['x'] * PoI_df['y'] * PoI_df['Mass [kg]']
    PoI_df['Iyz'] = PoI_df['y'] * PoI_df['z'] * PoI_df['Mass [kg]']
    PoI_df['Ixz'] = PoI_df['x'] * PoI_df['z'] * PoI_df['Mass [kg]']
    PoI = [sum(PoI_df['Ixy']) / mass_column, sum(PoI_df['Ixz']) / mass_column,
           sum(PoI_df['Ixy']) / mass_column]
    return PoI, PoI_df


def calc_centre_of_mass(mass_df):
    x = pd.DataFrame(mass_df['x'] * mass_df['weight_contribution'])
    x.columns = ['X_average']
    x['Y_average'] = pd.DataFrame(mass_df['y'] * mass_df['weight_contribution'])
    x['Z_average'] = pd.DataFrame(mass_df['z'] * mass_df['weight_contribution'])
    COM = [x['X_average'].sum(), x['Y_average'].sum(), x['Z_average'].sum()]
    return COM

def calc_radius_of_gyration(mass_df, COM):
    rog_df = pd.DataFrame(mass_df['x'] - COM[0])
    rog_df['y'] = mass_df['y'] - COM[1]
    rog_df['z'] = mass_df['z'] - COM[2]
    rog_df['Mass [kg]'] = mass_df['Mass [kg]']
    rog_df['I_x'] = rog_df['Mass [kg]'] * (np.square(rog_df['y']) + np.square(rog_df['z']))
    rog_df['I_y'] = rog_df['Mass [kg]'] * (np.square(rog_df['x']) + np.square(rog_df['z']))
    rog_df['I_z'] = rog_df['Mass [kg]'] * (np.square(rog_df['x']) + np.square(rog_df['y']))

    RoG = [np.sqrt(np.sum(rog_df['I_x']) / np.sum(mass_df['Mass [kg]'])),
                np.sqrt(np.sum(rog_df['I_y']) / np.sum(mass_df['Mass [kg]'])),
                np.sqrt(np.sum(rog_df['I_z']) / np.sum(mass_df['Mass [kg]']))]
    return rog_df, RoG

def calc_column_stiffness(rho, area, env, buoy, mufp, coord, COM):
    stiffness_matrix = np.zeros(shape=(6,6))
    stiffness_matrix[0, 0] = 1e6
    stiffness_matrix[1, 1] = 1e6
    stiffness_matrix[2, 2] = area.column * rho.water * env.g
    stiffness_matrix[5, 5] = 1e9

    I_wp = (np.pi / 64) * np.power(mufp.dia_column, 4)
    BM = I_wp / (buoy.heave + buoy.column)
    column_buoy = buoy.column + buoy.heave

    column_buoy_df = coord.buoy_df.loc[['Front Column Buoy', 'Front Heave Buoy']]


    COB = sum(column_buoy_df['z'] * column_buoy_df['Buoy [kg]']) / column_buoy

    KM = BM + (COB - COM[2])
    stiffness_matrix[3, 3] = KM * rho.water * env.g * (buoy.displaced_volume / 3)
    stiffness_matrix[4, 4] = KM * rho.water * env.g * (buoy.displaced_volume / 3)

    return stiffness_matrix

