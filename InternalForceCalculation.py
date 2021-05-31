import numpy as np
from scipy import interpolate
import pandas as pd


class InternalForces:

    def __init__(self, column_pressure_panels, force_column, mass_column_df, RAO, sim, f_rad, area, rho, env):
        # Creating empty matrices
        self.internal_force = np.zeros(shape=(sim.numheadangles, len(RAO[0, :, 0]), 6)) * 1j
        self.mass_column = sum(mass_column_df['Mass [kg]'])

        self.COM = calc_centre_of_mass(mass_column_df)
        self.rog_df, self.RoG = calc_radius_of_gyration(mass_column_df, self.COM)
        self.mom_of_inert = np.multiply(self.mass_column, np.square(self.RoG))

        # Interpolating column force data to match RAO sampling frequency
        force_column_interp = np.zeros(shape=(sim.numheadangles, len(RAO[0, :, 0]), 9))
        for ii, _ in enumerate(sim.front_column_force[:, 0, 0]):
            for jj, _ in enumerate(sim.front_column_force[0, 0, :]):
                force_column_interp[ii, :, jj] = np.interp(f_rad, sim.wave_disc[:, 4], force_column[ii, :, jj])

        self.force_column = force_column_interp

        # Calculating internal forces transfer functions
        self.translation_internal_forces(RAO, area, rho, env)
        self.rotational_internal_forces(column_pressure_panels, RAO, sim)

    def translation_internal_forces(self, RAO, area, rho, env):
        # Code returning the internal force applied to a column in each of the translational degrees of freedom
        # Internal Force Surge
        self.internal_force[0, :, 0] = (np.real(RAO[2, :, 0]) * self.mass_column - self.force_column[0, :, 0]) \
                                       + 1j * (np.imag(RAO[2, :, 0]) * self.mass_column - self.force_column[0, :, 1])
        # Internal Force Sway
        self.internal_force[0, :, 1] = (np.real(RAO[2, :, 1]) * self.mass_column - self.force_column[0, :, 3]) \
                                       + 1j * (np.imag(RAO[2, :, 1]) * self.mass_column - self.force_column[0, :, 4])
        # Internal Force Heave
        self.internal_force[0, :, 2] = (np.real(RAO[2, :, 2]) * self.mass_column - self.force_column[0, :, 6]) \
                                       + 1j * (np.imag(RAO[2, :, 2]) * self.mass_column - self.force_column[0, :, 7])

    def rotational_internal_forces(self, column_pressure_panels, RAO, sim):
        # Rotational Forces
        # Internal forces roll

        self.rot_force_column = np.zeros(shape=(2, len(RAO[0, :, 0]), 3)) * 1j
        self.rot_force_column[:, :, 0] = rot_calc(column_pressure_panels, RAO, sim, 'roll', self.COM)
        self.rot_force_column[:, :, 1] = rot_calc(column_pressure_panels, RAO, sim, 'pitch', self.COM)
        self.rot_force_column[:, :, 2] = rot_calc(column_pressure_panels, RAO, sim, 'yaw',  self.COM)

        # Finding internal forces
        self.internal_force[0, :, 3] = (np.real(RAO[2, :, 3]) * self.mom_of_inert[0] - np.real(self.rot_force_column[0, :, 0])) \
                                       + 1j * (np.imag(RAO[2, :, 3]) * self.mom_of_inert[0] - np.imag(self.rot_force_column[0, :, 0]))
        self.internal_force[0, :, 4] = (np.real(RAO[2, :, 4]) * self.mom_of_inert[1] - np.real(self.rot_force_column[0, :, 1])) \
                                       + 1j * (np.imag(RAO[2, :, 4]) * self.mom_of_inert[1] - np.imag(self.rot_force_column[0, :, 1]))
        self.internal_force[0, :, 5] = (np.real(RAO[2, :, 5]) * self.mom_of_inert[2] - np.real(self.rot_force_column[0, :, 2])) \
                                       + 1j * (np.imag(RAO[2, :, 5]) * self.mom_of_inert[2] - np.imag(self.rot_force_column[0, :, 2]))


def rot_calc(column_pressure_panels, RAO, sim, dof, COM):
    # Rotational Forces
    # Creating empty force matrix
    force = np.zeros(shape=(len(column_pressure_panels[:, 0, 0, 0]), len(column_pressure_panels[0, :, 0, 0]), 6))
    force_rao = np.zeros(shape=(2, len(RAO[0, :, 0]))) * 1j
    dof_f1 = np.zeros(shape=(
        len(sim.front_column[:, 0, 0, 0]), len(column_pressure_panels[0, :, 0, 0]), len(column_pressure_panels[0, 0, :, 4]), 3))
    dof_f2 = np.zeros(dof_f1.shape)
    if dof == 'roll':
        force_index = [15, 18]
        momentarm_index = [6, 5]
        com_index = [2,1]
    elif dof == 'pitch':
        force_index = [12, 18]
        momentarm_index = [6, 4]
        com_index = [2,0]
    elif dof == 'yaw':
        force_index = [12, 15]
        momentarm_index = [5, 4]
        com_index = [1,0]
    else:
        print('invalid degree of freedom')


    for jj, _ in enumerate(column_pressure_panels[:, 0, 0, 0]):
        for kk, _ in enumerate(column_pressure_panels[0, :, 0, 0]):
            for ii in np.arange(0, 3, 1):
                dof_f1[jj, kk, :, ii] = column_pressure_panels[jj, kk, :,
                                        force_index[0] + ii] * (column_pressure_panels[0, 0, :,
                                                               momentarm_index[0]] - COM[com_index[0]])
                dof_f2[jj, kk, :, ii] = column_pressure_panels[jj, kk, :,
                                        force_index[1] + ii] * (column_pressure_panels[0, 0, :,
                                                               momentarm_index[1]] - COM[com_index[1]])

    for jj, _ in enumerate(column_pressure_panels[:, 0, 0, 0]):
        for kk, _ in enumerate(column_pressure_panels[0, :, 0, 0]):
                force[jj, kk, 0:3] = np.sum(dof_f1[jj, kk, :, 0:3], axis=0)
                force[jj, kk, 3:6] = np.sum(dof_f2[jj, kk, :, 0:3], axis=0)
        force_rao[jj, kk] = (force[jj, kk, 0] + force[jj, kk, 3]) + 1j * (force[jj, kk, 1] + force[jj, kk, 4])

    return force_rao


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
    rog_df['I_x'] = rog_df['Mass [kg]'] * np.square(
        np.sqrt(np.square(rog_df['y']) + np.square(rog_df['z'])))
    rog_df['I_y'] = rog_df['Mass [kg]'] * np.square(
        np.sqrt(np.square(rog_df['x']) + np.square(rog_df['z'])))
    rog_df['I_z'] = rog_df['Mass [kg]'] * np.square(
        np.sqrt(np.square(rog_df['x']) + np.square(rog_df['y'])))

    RoG = [np.sqrt(np.sum(rog_df['I_x']) / np.sum(mass_df['Mass [kg]'])),
                np.sqrt(np.sum(rog_df['I_y']) / np.sum(mass_df['Mass [kg]'])),
                np.sqrt(np.sum(rog_df['I_z']) / np.sum(mass_df['Mass [kg]']))]
    return rog_df, RoG