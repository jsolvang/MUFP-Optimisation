import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle
import sys
from InternalForceCalculation import InternalForces, rot_calc
from RAO_Calculation import calulate_RAOs, calculate_column_raos
from FloaterParameters import FloaterParameters
from Environment import Environment
from Buoy import Buoy
from Mass import Mass
from Density import Density
from Area import Area
from SystemMatrices import MatrixCalculation
from GlobalCoordinateSystem import GlobalCoordinateSystem

file_loc = r'C:\Users\Joar\Documents\1_Education\NTNU'
function_path = file_loc + r'\Functions'
sys.path.append(function_path)

if __name__ == "__main__":
    # Switch to 1 if the answers plotted
    plot = 1
    # Switch to see the plots for the second
    headangle = 0

    # Manually defining floater parameters - These functions are only used to define distances while calculating
    # individual column and RNA raos.
    pull = np.array([60, 80, 12, 14])
    # Hand Calculations
    mufp = FloaterParameters(pull[0], pull[1], pull[2], pull[3])
    env = Environment()
    rho = Density()
    csa = Area(mufp)
    buoy = Buoy(mufp, csa, rho)
    mass = Mass(mufp, csa, buoy, rho)
    coord = GlobalCoordinateSystem(mufp, csa, mass, rho, buoy, env)
    matrix = MatrixCalculation(coord, mass, mufp, rho, env, csa, buoy)
    matrix.stiffness[0, 0] = 1e6
    matrix.stiffness[1, 1] = 1e6
    matrix.stiffness[5, 5] = 1e9

    ####### Comparing Dispersion Pressure Data to HydroD Wave excitation Data ###################
    #############################################################################################

    # Importing HydroD Simulation where model is constrained by high stiffnesses in all DOF
    file_name = r"\pickle_files\fixed_motion"
    file_path = file_loc + file_name
    infile = open(file_path, 'rb')
    fixed_motion = pickle.load(infile)
    infile.close()

    # Calculating forces acting of floater due to panel pressure
    full_body_force = np.zeros(shape=(2, 84, 3)) + 0j
    for ii in np.arange(0, fixed_motion.numheadangles, 1).astype(int):
        for jj in np.arange(0, fixed_motion.numwavelengths, 1).astype(int):
            for kk in [0, 1, 2]:
                full_body_force[ii, jj, kk] = np.sum(fixed_motion.panel_pressure[ii, jj, :, 12 + kk])

    # Calculating moments acting on floater
    rot_force_column = np.zeros(shape=(2, 84, 3)) + 0j
    rot_force_column[:, :, 0] = rot_calc(fixed_motion.panel_pressure, fixed_motion, 'roll', [0, 0, 0])
    rot_force_column[:, :, 1] = rot_calc(fixed_motion.panel_pressure, fixed_motion, 'pitch', [0, 0, 0])
    rot_force_column[:, :, 2] = rot_calc(fixed_motion.panel_pressure, fixed_motion, 'yaw', [0, 0, 0])

    # Code plotting forces due to panel pressure compared to wave excitation extracted directly from HydroD
    if plot == 1:
        ylab = np.array([['Surge [N/m]', 'Sway [N/m]', 'Heave[N/m]'],
                         ['Roll [Nm/m]', 'Pitch [Nm/m]', 'Yaw [Nm/m]'], ])
        plt.rcParams["figure.figsize"] = (15, 15)
        plt.rcParams.update({'font.size': 15})
        fig, axs1 = plt.subplots(3, 2)
        for jj in [0, 1, 2]:
            axs1[jj, 0].plot(fixed_motion.wave_disc[:, 4], np.imag(full_body_force[headangle, :, jj]), '-xb',
                             label='Real - Manual Panel Pressure', linewidth=2)
            axs1[jj, 0].scatter(fixed_motion.wave_disc[:, 4], fixed_motion.WAVEEX[headangle, :, jj, 1],
                                label='Real - HydroD Excitation', s=100, facecolors='none', edgecolors='r')
            axs1[jj, 0].set(ylabel=ylab[0, jj])
            axs1[jj, 0].set(xlim=(0, 4))
            axs1[2, 0].set(xlabel=('Frequency [rad/s]'))
            axs1[jj, 0].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        for jj in [3, 4, 5]:
            axs1[jj - 3, 1].plot(fixed_motion.wave_disc[:, 4], np.imag(rot_force_column[headangle, :, jj - 3]), '-xb',
                                 linewidth=2)
            axs1[jj - 3, 1].scatter(fixed_motion.wave_disc[:, 4], fixed_motion.WAVEEX[headangle, :, jj, 1], s=100,
                                    facecolors='none', edgecolors='r')
            axs1[2, 1].set(xlabel=('Frequency [rad/s]'))
            axs1[jj - 3, 1].set(ylabel=ylab[1, jj - 3])
            axs1[jj - 3, 1].set(xlim=(0, 4))
            axs1[jj - 3, 1].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        plt.figlegend(['$\sum$Pressure Panel Forces', 'HydroD Excitation'], bbox_to_anchor=(.82, .8))

    ####### Comparing Manual Internal Force Calculation to HydroD Sectional Loads ###############
    #############################################################################################

    # Importing "Free motion" simulation from pickle file
    file_name = r"\pickle_files\free_motion"
    file_path = file_loc + file_name
    infile = open(file_path, 'rb')
    hydrod_loads = pickle.load(infile)
    infile.close()

    # Sign convention for internal loads is opposite
    hydrod_loads.secloads[:, :, :, 0:6, :] *= -1

    # Calculating rigid body RAO.
    RAO_rigid, _, _, _ = calulate_RAOs(hydrod_loads.wave_disc, hydrod_loads.mass_matrix, hydrod_loads.ADDEDMASS,
                                       hydrod_loads.DAMPING, hydrod_loads.stiffness_matrix, hydrod_loads.WAVEEX, 0,
                                       mufp, coord)

    # Mapping rigid body RAO to single column RAOs
    RAO = np.zeros(shape=(3, 2, 3, 84, 8)) + 0j
    RAO[0, :, :, :, :] = calculate_column_raos(RAO_rigid, matrix, 'front')
    RAO[1, :, :, :, :] = calculate_column_raos(RAO_rigid, matrix, 'left')
    RAO[2, :, :, :, :] = calculate_column_raos(RAO_rigid, matrix, 'right')

    # Code plotting RAOs
    if plot == 1:
        sec = 0  # sec = sectional_load -- 0 is front column, 1 is left, 2 is right
        Y_lab = ["RAO Surge, [m/m]",
                 "RAO Sway, [m/m]",
                 "RAO Heave, [m/m]",
                 "RAO Roll, [m/m]",
                 "RAO Pitch, [m/m]",
                 "RAO Yaw, [m/m]"]
        plt.rcParams["figure.figsize"] = (15, 15)
        plt.rcParams.update({'font.size': 15})
        plt.rcParams['font.family'] = 'serif'
        plot_dof = [0, 1, 2, 3, 4, 5]
        fig, axs = plt.subplots(len(plot_dof))
        for x, jj in enumerate(plot_dof):
            axs[x].plot(hydrod_loads.wave_disc[:, 4], np.real(RAO[sec, headangle, 0, :, jj]), '-xb',
                        label="Real - Hand Calc", linewidth=2, markersize=4)
            axs[x].plot(hydrod_loads.wave_disc[:, 4], np.imag(RAO[sec, headangle, 0, :, jj]), '--Dr',
                        label="Imag - Hand Calc", linewidth=2, markersize=4)
            axs[x].plot(hydrod_loads.wave_disc[:, 4], np.absolute(RAO[sec, headangle, 0, :, jj]), ':.g',
                        label="Absolute - Hand Calc", linewidth=2, markersize=4)
            axs[x].set(xlim=[0, 5])
            if x < 3:
                axs[x].set(ylim=[-5, 5])
            else:
                axs[x].set(ylim=[-0.1, 0.1])
            axs[x].set(xlabel='Frequency, rad/s', ylabel=Y_lab[jj])
            axs[x].legend(loc="upper right", bbox_to_anchor=(1.3, 1.01))
            axs[x].grid(b=True, which='both', axis='both')
        plt.rcParams["figure.figsize"] = (15, 5)

    # Calculating internal forces in the system
    IF_manual = []
    IF_manual.append(InternalForces(hydrod_loads.front_column, hydrod_loads.front_column_force, coord.front_column_df,
                                    RAO[0, :, :, :, :], hydrod_loads, hydrod_loads.wave_disc[:, 4], csa, rho, env, buoy,
                                    mufp, coord, mass))
    IF_manual.append(InternalForces(hydrod_loads.left_column, hydrod_loads.left_column_force, coord.left_column_df,
                                    RAO[1, :, :, :, :], hydrod_loads, hydrod_loads.wave_disc[:, 4], csa, rho, env, buoy,
                                    mufp, coord, mass))
    IF_manual.append(InternalForces(hydrod_loads.right_column, hydrod_loads.right_column_force, coord.right_column_df,
                                    RAO[2, :, :, :, :], hydrod_loads, hydrod_loads.wave_disc[:, 4], csa, rho, env, buoy,
                                    mufp, coord, mass))

    # Code plotting Internal forces compared to sectional loads
    # Currently only plotting internal forces, however uncomment lines if you want inertial/pressure panels/stiffness plots
    if plot == 1:
        sec = 0  # sec = sectional_load -- 0 is front column, 1 is left, 2 is right
        y_lab = ["Fx, [N/m]", "Fy, [N/m]", "Fz, [N/m]", "Mx, [Nm/m]", "My, [Nm/m]", "Mz, [Nm/m]"]
        plt.rcParams["figure.figsize"] = (15, 20)
        plt.rcParams.update({'font.size': 15})
        fig, axs1 = plt.subplots(6)
        for jj in [0, 1, 2, 3, 4, 5]:
            axs1[jj].scatter(hydrod_loads.wave_disc[:, 4], np.real(hydrod_loads.secloads[headangle, :, sec, jj]),
                             label='HydroD: Internal Force', s=100, facecolors='none', edgecolors='r')
            axs1[jj].plot(hydrod_loads.wave_disc[:, 4], np.real(IF_manual[sec].internal_force[headangle, :, jj]), '-*b',
                          label='Manual Calculation: Internal Force', linewidth=3)
            # if jj < 3:
            #    axs1[jj].plot(hydrod_loads.wave_disc[:,4], np.real(RAO[sec, 2, :, jj]) *IF_manual[sec].mass_column, 'b', label='Acceleration')
            #    axs1[jj].plot(hydrod_loads.wave_disc[:,4], -np.real(IF_manual[sec].force_column[0, :, jj]), 'g', label='Pressure Panels')
            # else:
            #    axs1[jj].plot(hydrod_loads.wave_disc[:,4], np.real(RAO[sec,headangle, 2, :, jj]) * IF_manual[sec].mass_matrix[jj,jj], 'b', label='Acceleration')
            #    axs1[jj].plot(hydrod_loads.wave_disc[:,4], -np.real(IF_manual[sec].rot_force_column[0, :, jj-3]), 'g', label='Pressure Panels')
            # axs1[jj].plot(hydrod_loads.wave_disc[:,4], np.real(RAO[sec,0, :, jj])*IF_manual[sec].stiffness_matrix[jj,jj], 'r', label='Stiffness')
            axs1[jj].set(ylabel=y_lab[jj])
            axs1[jj].set(xlim=(0, 4))
            axs1[2].set(xlabel=('Frequency [rad/s]'))
            axs1[jj].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        # plt.figlegend(['Manual Calculation','Inertial Load','$\sum$ Pressure Panels','Stiffness', 'HydroD: Internal Force'], bbox_to_anchor=(.82, .77))
        plt.figlegend(['HydroD - Sectional Loads', 'Internal Force Calculation'], bbox_to_anchor=(0.95, 0.97))
        axs1[0].set(ylim=(-1e6, 1e6))
        axs1[1].set(ylim=(-5e6, 5e6))
        axs1[2].set(ylim=(-6e5, 6e5))
        axs1[3].set(ylim=(-4e7, 4e7))
        axs1[4].set(ylim=(-4e7, 3e8))
        axs1[5].set(ylim=(-2e8, 2e8))
        plt.tight_layout()

    plt.show()
