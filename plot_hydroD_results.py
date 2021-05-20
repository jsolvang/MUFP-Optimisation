import numpy as np
import subprocess
import os
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from kSolve import ksolve

def plot_hydroD_results(results, floater, env, plot_sync):
    SMALL_SIZE = 8
    MEDIUM_SIZE = 10
    BIGGER_SIZE = 12

    if plot_sync == 1:
        alpha_plot = 0.6
        L_x = np.zeros(21)
        k_x = np.zeros(21)

        L_y = np.zeros(21)
        k_y = np.zeros(21)

        omega_list = np.arange(0.01, 5, 0.01)
        k_list = np.zeros(len(omega_list))
        L_list = np.zeros(len(omega_list))

        phi_x = np.zeros(len(omega_list))
        phi_y = np.zeros(len(omega_list))


        for ii in np.linspace(0, len(omega_list) - 1, len(omega_list)).astype(int):
            [k_list[ii], L_list[ii]] = ksolve(omega_list[ii], 50, 9.81)
            phi_x[ii] = (floater.x_space / L_list[ii]) - np.floor(floater.x_space / L_list[ii])
            phi_y[ii] = (floater.y_space / L_list[ii]) - np.floor(floater.y_space / L_list[ii])

        phi_x = phi_x * 360 - 180
        phi_y = phi_y * 360 - 180

        for i in np.linspace(1, 20, 20).astype(int):
           L_x[i] = floater.x_space / i
           k_x[i] = 2 * np.pi / L_x[i]
        omega_x = np.sqrt(env.g * k_x * np.tanh(k_x * env.h))
        omega_x = np.delete(omega_x, np.where(omega_x < 0.5))
        omega_x = np.delete(omega_x, np.where(omega_x > 2.5))

        for i in np.linspace(1, 20, 20).astype(int):
            L_y[i] = floater.y_space / i
            k_y[i] = 2 * np.pi / L_y[i]
        omega_y = np.sqrt(env.g * k_y * np.tanh(k_y * env.h))
        omega_y = np.delete(omega_y, np.where(omega_y < 0.5))
        omega_y = np.delete(omega_y, np.where(omega_y > 2.5))

    plt.rcParams.update({'font.size': 20})

    plt.rcParams["figure.figsize"] = (20, 20)
    styles = [['-*b', '-*g', '-*r', '-*b', '-*g', '-*r'],
              ['--b', '--g', '--r', '--b', '--g', '--r'],
              [':b', ':g', ':r', ':b', ':g', ':r'],
              ['-b', '-g', '-r', '-b', '-g', '-r'],
              ['--b', '--g', '--r', '--b', '--g', '--r'],
              [':b', ':g', ':r', ':b', ':g', ':r']]

    legend_A = [['A11', 'A12', 'A13', 'A14', 'A15', 'A16'],
                ['A21', 'A22', 'A23', 'A24', 'A25', 'A26'],
                ['A31', 'A32', 'A33', 'A34', 'A35', 'A36'],
                ['A41', 'A42', 'A43', 'A44', 'A45', 'A46'],
                ['A51', 'A52', 'A53', 'A54', 'A55', 'A56'],
                ['A61', 'A62', 'A63', 'A64', 'A65', 'A66']]

    legend_D = [['D11', 'D12', 'D13', 'D14', 'D15', 'D16'],
                ['D21', 'D22', 'D23', 'D24', 'D25', 'D26'],
                ['D31', 'D32', 'D33', 'D34', 'D35', 'D36'],
                ['D41', 'D42', 'D43', 'D44', 'D45', 'D46'],
                ['D51', 'D52', 'D53', 'D54', 'D55', 'D56'],
                ['D61', 'D62', 'D63', 'D64', 'D65', 'D66']]

    legend_F = ['|$X_1$|', '|$X_2$|', '|$X_3$|', '|$X_4$|', '|$X_5$|', '|$X_6$|']

    legend_T = [r'$\theta_1$', r'$\theta_2$', r'$\theta_3$', r'$\theta_4$', r'$\theta_5$', r'$\theta_6$']

    fig, axs1 = plt.subplots(2, 2)
    for ii in np.linspace(0, 2, 3).astype(int):
        for jj in np.linspace(0, 2, 3).astype(int):
            axs1[0, 0].plot(results.wave_disc[:, 4], results.ADDEDMASS[:, ii, jj], styles[ii][jj], label=legend_A[ii][jj],
                           linewidth=2, markersize=5)
    if plot_sync == 1:
        for xc in omega_x:
            axs1[0, 0].axvline(x=xc, linestyle='--', alpha=alpha_plot, linewidth=2, color='c')
        for xx in omega_y:
            axs1[0, 0].axvline(x=xx, color='m', linestyle='--', alpha=alpha_plot, linewidth=2)
    axs1[0, 0].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    axs1[0, 0].legend(loc="upper right")
    axs1[0, 0].set(xlabel='omega, rad/s', ylabel='Added Mass [kg]')
    axs1[0, 0].grid(b=True, which='both', axis='both')


    for ii in np.linspace(3, 5, 3).astype(int):
        for jj in np.linspace(0, 2, 3).astype(int):
            axs1[1, 0].plot(results.wave_disc[:, 4], results.ADDEDMASS[:, ii, jj], styles[ii][jj], label=legend_A[ii][jj],
                           linewidth=2, markersize=5)
    if plot_sync == 1:
        for xc in omega_x:
            axs1[1, 0].axvline(x=xc, linestyle='--', alpha=alpha_plot, linewidth=2, color='c')
        for xx in omega_y:
            axs1[1, 0].axvline(x=xx, color='m', linestyle='--', alpha=alpha_plot, linewidth=2)
    axs1[1, 0].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    axs1[1, 0].legend(loc="upper right")
    axs1[1, 0].set(xlabel='omega, rad/s', ylabel='Added Mass [kg-m^2]')
    axs1[1, 0].grid(b=True, which='both', axis='both')

    for ii in np.linspace(0, 2, 3).astype(int):
        for jj in np.linspace(3, 5, 3).astype(int):
            axs1[1, 1].plot(results.wave_disc[:, 4], results.ADDEDMASS[:, ii, jj], styles[ii][jj], label=legend_A[ii][jj],
                           linewidth=2, markersize=5)
    if plot_sync == 1:
        for xc in omega_x:
            axs1[1, 1].axvline(x=xc, linestyle='--', alpha=alpha_plot, linewidth=2, color='c')
        for xx in omega_y:
            axs1[1, 1].axvline(x=xx, color='m', linestyle='--', alpha=alpha_plot, linewidth=2)
    axs1[1, 1].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    axs1[1, 1].legend(loc="upper right")
    axs1[1, 1].set(xlabel='omega, rad/s', ylabel='Added Mass [kg-m]')
    axs1[1, 1].grid(b=True, which='both', axis='both')

    for ii in np.linspace(3, 5, 3).astype(int):
        for jj in np.linspace(3, 5, 3).astype(int):
            axs1[0, 1].plot(results.wave_disc[:, 4], results.ADDEDMASS[:, ii, jj], styles[ii][jj], label=legend_A[ii][jj],
                           linewidth=2, markersize=5)
    if plot_sync == 1:
        for xc in omega_x:
            axs1[0, 1].axvline(x=xc, linestyle='--', alpha=alpha_plot, linewidth=2, color='c')
        for xx in omega_y:
            axs1[0, 1].axvline(x=xx, color='m', linestyle='--', alpha=alpha_plot, linewidth=2)
    axs1[0, 1].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    axs1[0, 1].legend(loc="upper right")
    axs1[0, 1].set(xlabel='omega, rad/s', ylabel='Added Mass [kg-m]')
    axs1[0, 1].grid(b=True, which='both', axis='both')

    fig, axs2 = plt.subplots(2, 2)
    for ii in np.linspace(0, 2, 3).astype(int):
        for jj in np.linspace(0, 2, 3).astype(int):
            axs2[0, 0].plot(results.wave_disc[:, 4], results.DAMPING[:, ii, jj], styles[ii][jj], label=legend_D[ii][jj],
                           linewidth=2, markersize=2)
    if plot_sync == 1:
        for xc in omega_x:
            axs2[0, 0].axvline(x=xc, linestyle='--', alpha=alpha_plot, linewidth=2, color='c')
        for xx in omega_y:
            axs2[0, 0].axvline(x=xx, color='m', linestyle='--', alpha=alpha_plot, linewidth=2)
    axs2[0, 0].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    axs2[0, 0].legend(loc="upper right")
    axs2[0, 0].set(xlabel='omega, rad/s', ylabel='Damping [kg/s]')
    axs2[0, 0].grid(b=True, which='both', axis='both')

    for ii in np.linspace(0, 2, 3).astype(int):
        for jj in np.linspace(3, 5, 3).astype(int):
            axs2[1, 1].plot(results.wave_disc[:, 4], results.DAMPING[:, ii, jj], styles[ii][jj], label=legend_D[ii][jj],
                           linewidth=2, markersize=2)
    if plot_sync == 1:
        for xc in omega_x:
            axs2[1, 1].axvline(x=xc, linestyle='--', alpha=alpha_plot, linewidth=2, color='c')
        for xx in omega_y:
            axs2[1, 1].axvline(x=xx, color='m', linestyle='--', alpha=alpha_plot, linewidth=2)
    axs2[1, 1].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    axs2[1, 1].legend(loc="upper right")
    axs2[1, 1].set(xlabel='omega, rad/s', ylabel='Damping [kg-m/s^2]')
    axs2[1, 1].grid(b=True, which='both', axis='both')

    for ii in np.linspace(3, 5, 3).astype(int):
        for jj in np.linspace(0, 2, 3).astype(int):
            axs2[1, 0].plot(results.wave_disc[:, 4], results.DAMPING[:, ii, jj], styles[ii][jj], label=legend_D[ii][jj],
                           linewidth=2, markersize=2)
    if plot_sync == 1:
        for xc in omega_x:
            axs2[1, 0].axvline(x=xc, linestyle='--', alpha=alpha_plot, linewidth=2, color='c')
        for xx in omega_y:
            axs2[1, 0].axvline(x=xx, color='m', linestyle='--', alpha=alpha_plot, linewidth=2)
    axs2[1, 0].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    axs2[1, 0].legend(loc="upper right")
    axs2[1, 0].set(xlabel='omega, rad/s', ylabel='Damping [kg-m/s]')
    axs2[1, 0].grid(b=True, which='both', axis='both')

    for ii in np.linspace(3, 5, 3).astype(int):
        for jj in np.linspace(3, 5, 3).astype(int):
            axs2[0, 1].plot(results.wave_disc[:, 4], results.DAMPING[:, ii, jj], styles[ii][jj], label=legend_D[ii][jj],
                           linewidth=2, markersize=2)
    if plot_sync == 1:
        for xc in omega_x:
            axs2[0, 1].axvline(x=xc, linestyle='--', alpha=alpha_plot, linewidth=2, color='c')
        for xx in omega_y:
            axs2[0, 1].axvline(x=xx, color='m', linestyle='--', alpha=alpha_plot, linewidth=2)
    axs2[0, 1].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    axs2[0, 1].legend(loc="upper right")
    axs2[0, 1].set(xlabel='omega, rad/s', ylabel='Damping [kg-m/s]')
    axs2[0, 1].grid(b=True, which='both', axis='both')

    fig, axs3 = plt.subplots(2, 2)
    for jj in [0,1,2]:
        axs3[0, 0].plot(results.wave_disc[:, 4], results.WAVEEX[0, :, jj, 2], styles[0][jj], label=legend_F[jj], linewidth=2, markersize=5)
    if plot_sync == 1:
        for xc in omega_x:
            axs3[0, 0].axvline(x=xc, linestyle='--', alpha=alpha_plot, linewidth=2, color='c')
        for xx in omega_y:
            axs3[0, 0].axvline(x=xx, color='m', linestyle='--', alpha=alpha_plot, linewidth=2)
    axs3[0, 0].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    axs3[0, 0].legend(loc="upper right")
    axs3[0, 0].set(xlabel='omega, rad/s', ylabel='0 Deg Exciting Force Translation [N]')
    axs3[0, 0].grid(b=True, which='both', axis='both')

    for jj in [3,4,5]:
        axs3[1, 0].plot(results.wave_disc[:, 4], results.WAVEEX[0, :, jj, 2], styles[0][jj], label=legend_F[jj], linewidth=2, markersize=5)
    if plot_sync == 1:
        for xc in omega_x:
            axs3[1, 0].axvline(x=xc, linestyle='--', alpha=alpha_plot, linewidth=2, color='c')
        for xx in omega_y:
            axs3[1, 0].axvline(x=xx, color='m', linestyle='--', alpha=alpha_plot, linewidth=2)
    axs3[1, 0].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    axs3[1, 0].legend(loc="upper right")
    axs3[1, 0].set(xlabel='omega, rad/s', ylabel='0 Deg Exciting Force Rotation [Nm]')
    axs3[1, 0].grid(b=True, which='both', axis='both')

    for jj in [0,2]:
        axs3[0, 1].plot(results.wave_disc[:, 4], results.WAVEEX[0, :, jj, 3], '*', label=legend_F[jj], linewidth=2, markersize=5)
    # if plot_sync == 1:
        # axs[0, 1].plot(omega_list, phi_x, linestyle='-', alpha=0.3, linewidth=2)
    axs3[0, 1].legend(loc="upper right")
    axs3[0, 1].set(xlabel='omega, rad/s', ylabel='0 Deg - Phase Shift - Translation [N]')
    axs3[0, 1].grid(b=True, which='both', axis='both')

    for jj in [4]:
        axs3[1, 1].plot(results.wave_disc[:, 4], results.WAVEEX[0, :, jj, 3], '*g', label=legend_F[jj], linewidth=2, markersize=5)
    # if plot_sync == 1:
        # axs[1, 1].plot(omega_list,phi_x, linestyle='-', alpha=0.3, linewidth=2)
    axs3[1, 1].legend(loc="upper right")
    axs3[1, 1].set(xlabel='omega, rad/s', ylabel='0 Deg - Phase Shift - Rotation [Nm]')
    axs3[1, 1].grid(b=True, which='both', axis='both')

#
#    fig, axs4 = plt.subplots(2, 2)
#    for jj in [0,2]:
#        axs4[0, 0].plot(results.wave_disc[:, 4], results.WAVEEX[1, :, jj, 2], label=legend_F[jj], linewidth=2, markersize=2)
#        if plot_sync == 1:
#            for xc in omega_y:
#                axs4[0, 0].axvline(x=xc, linestyle='-', alpha=0.3, linewidth=2)
#        axs4[0, 0].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
#        axs4[0, 0].legend(loc="upper right")
#        axs4[0, 0].set(xlabel='omega, rad/s', ylabel='90 Deg Exciting Force Translation [N]')
#        axs4[0, 0].grid(b=True, which='both', axis='both')
#
#    for jj in [4]:
#        axs4[1, 0].plot(results.wave_disc[:, 4], results.WAVEEX[1, :, jj, 2], '-g', label=legend_T[jj], linewidth=2, markersize=2)
#        if plot_sync == 1:
#            for xc in omega_y:
#                axs4[1, 0].axvline(x=xc, linestyle='-', alpha=0.3, linewidth=2)
#        axs4[1, 0].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
#        axs4[1, 0].legend(loc="upper right")
#        axs4[1, 0].set(xlabel='omega, rad/s', ylabel='90 Deg Exciting Force Rotation [N]')
#        axs4[1, 0].grid(b=True, which='both', axis='both')
#
#    for jj in [0,2]:
#        axs4[0, 1].plot(results.wave_disc[:, 4], results.WAVEEX[1, :, jj, 3], '*', label=legend_F[jj], linewidth=2, markersize=5)
#        # if plot_sync == 1:
#            # axs[0, 1].plot(omega_list, phi_y, linestyle='-', alpha=0.3, linewidth=2)
#        axs4[0, 1].legend(loc="upper right")
#        axs4[0, 1].set(xlabel='omega, rad/s', ylabel='90 Deg - Phase Shift -  Translation [N]')
#        axs4[0, 1].grid(b=True, which='both', axis='both')

#    for jj in [4]:
#        axs4[1, 1].plot(results.wave_disc[:, 4], results.WAVEEX[1, :, jj, 3], '*g', label=legend_F[jj], linewidth=2, markersize=5)
#        # if plot_sync == 1:
#            # axs[1, 1].plot(omega_list, phi_y, linestyle='-', alpha=0.3, linewidth=2)
#        axs4[1, 1].legend(loc="upper right")
#        axs4[1, 1].set(xlabel='omega, rad/s', ylabel='90 Deg - Phase Shift - Rotation [N]')
#        axs4[1, 1].grid(b=True, which='both', axis='both')

