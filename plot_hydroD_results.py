import numpy as np
import subprocess
import os
import matplotlib.pyplot as plt

def plot_hydroD_results(results):
    SMALL_SIZE = 8
    MEDIUM_SIZE = 10
    BIGGER_SIZE = 12

    plt.rc('font', size=BIGGER_SIZE)  # controls default text sizes
    plt.rc('axes', titlesize=BIGGER_SIZE)  # fontsize of the axes title
    plt.rc('axes', labelsize=BIGGER_SIZE)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=BIGGER_SIZE)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=BIGGER_SIZE)  # fontsize of the tick labels
    plt.rc('legend', fontsize=BIGGER_SIZE)  # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    plt.rcParams["figure.figsize"] = (20, 20)
    styles = [['-b', '-g', '-r', '-b', '-g', '-r'],
              ['--b', '--g', '--r', '--b', '--g', '--r'],
              [':b', ':g', ':r', ':b', ':g', ':r'],
              ['-.b', '-.g', '-.r', '-.b', '-.g', '-.r'],
              ['bo--', 'go--', 'ro--', 'bo--', 'go--', 'ro--'],
              ['b+--', 'g+--', 'r+--', 'b+--', 'g+--', 'r+--']]

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

    fig, axs = plt.subplots(2, 2)
    for ii in np.linspace(0, 2, 3).astype(int):
        for jj in np.linspace(0, 2, 3).astype(int):
            axs[0, 0].plot(results.wave_disc[:, 4], results.ADDEDMASS[:, ii, jj], styles[ii][jj], label=legend_A[ii][jj],
                           linewidth=2, markersize=5)
    axs[0, 0].legend(loc="upper right")
    axs[0, 0].set(ylabel='Added Mass [kg]')
    axs[0, 0].grid(b=True, which='both', axis='both')

    for ii in np.linspace(3, 5, 3).astype(int):
        for jj in np.linspace(0, 2, 3).astype(int):
            axs[1, 0].plot(results.wave_disc[:, 4], results.ADDEDMASS[:, ii, jj], styles[ii][jj], label=legend_A[ii][jj],
                           linewidth=2, markersize=5)
    axs[1, 0].legend(loc="upper right")
    axs[1, 0].set(ylabel='Added Mass [kg-m^2]')
    axs[1, 0].grid(b=True, which='both', axis='both')

    for ii in np.linspace(0, 2, 3).astype(int):
        for jj in np.linspace(3, 5, 3).astype(int):
            axs[1, 1].plot(results.wave_disc[:, 4], results.ADDEDMASS[:, ii, jj], styles[ii][jj], label=legend_A[ii][jj],
                           linewidth=2, markersize=5)
    axs[1, 1].legend(loc="upper right")
    axs[1, 1].set(xlabel='omega, rad/s', ylabel='Added Mass [kg-m]')
    axs[1, 1].grid(b=True, which='both', axis='both')

    for ii in np.linspace(3, 5, 3).astype(int):
        for jj in np.linspace(3, 5, 3).astype(int):
            axs[0, 1].plot(results.wave_disc[:, 4], results.ADDEDMASS[:, ii, jj], styles[ii][jj], label=legend_A[ii][jj],
                           linewidth=2, markersize=5)
    axs[0, 1].legend(loc="upper right")
    axs[0, 1].set(xlabel='omega, rad/s', ylabel='Added Mass [kg-m]')
    axs[0, 1].grid(b=True, which='both', axis='both')

    fig, axs = plt.subplots(2, 2)
    for ii in np.linspace(0, 2, 3).astype(int):
        for jj in np.linspace(0, 2, 3).astype(int):
            axs[0, 0].plot(results.wave_disc[:, 4], results.DAMPING[:, ii, jj], styles[ii][jj], label=legend_D[ii][jj],
                           linewidth=2, markersize=2)
    axs[0, 0].legend(loc="upper right")
    axs[0, 0].set(ylabel='Damping [kg/s]')
    axs[0, 0].grid(b=True, which='both', axis='both')

    for ii in np.linspace(0, 2, 3).astype(int):
        for jj in np.linspace(3, 5, 3).astype(int):
            axs[1, 1].plot(results.wave_disc[:, 4], results.DAMPING[:, ii, jj], styles[ii][jj], label=legend_D[ii][jj],
                           linewidth=2, markersize=2)
    axs[1, 1].legend(loc="upper right")
    axs[1, 1].set(ylabel='Damping [kg-m/s^2]')
    axs[1, 1].grid(b=True, which='both', axis='both')

    for ii in np.linspace(3, 5, 3).astype(int):
        for jj in np.linspace(0, 2, 3).astype(int):
            axs[1, 0].plot(results.wave_disc[:, 4], results.DAMPING[:, ii, jj], styles[ii][jj], label=legend_D[ii][jj],
                           linewidth=2, markersize=2)
    axs[1, 0].legend(loc="upper right")
    axs[1, 0].set(xlabel='omega, rad/s', ylabel='Damping [kg-m/s]')
    axs[1, 0].grid(b=True, which='both', axis='both')

    for ii in np.linspace(3, 5, 3).astype(int):
        for jj in np.linspace(3, 5, 3).astype(int):
            axs[0, 1].plot(results.wave_disc[:, 4], results.DAMPING[:, ii, jj], styles[ii][jj], label=legend_D[ii][jj],
                           linewidth=2, markersize=2)
    axs[0, 1].legend(loc="upper right")
    axs[0, 1].set(xlabel='omega, rad/s', ylabel='Damping [kg-m/s]')
    axs[0, 1].grid(b=True, which='both', axis='both')

    fig, axs = plt.subplots(2, 2)
    for jj in np.linspace(0, 2, 3).astype(int):
        axs[0, 0].plot(results.wave_disc[:, 4], results.WAVEEX[0, :, jj, 2], label=legend_F[jj], linewidth=2, markersize=2)
        axs[0, 0].legend(loc="upper right")
        axs[0, 0].set(xlabel='omega, rad/s', ylabel='0 Deg Exciting Force [N]')
        axs[0, 0].grid(b=True, which='both', axis='both')

    for jj in np.linspace(3, 5, 3).astype(int):
        axs[0, 1].plot(results.wave_disc[:, 4], results.WAVEEX[0, :, jj, 2], label=legend_T[jj], linewidth=2, markersize=2)
        axs[0, 1].legend(loc="upper right")
        axs[0, 1].set(xlabel='omega, rad/s', ylabel='0 Deg Exciting Force [N]')
        axs[0, 1].grid(b=True, which='both', axis='both')

    for jj in np.linspace(0, 2, 3).astype(int):
        axs[1, 0].plot(results.wave_disc[:, 4], results.WAVEEX[1, :, jj, 2], label=legend_F[jj], linewidth=2, markersize=2)
        axs[1, 0].legend(loc="upper right")
        axs[1, 0].set(xlabel='omega, rad/s', ylabel='90 Deg Exciting Force [N]')
        axs[1, 0].grid(b=True, which='both', axis='both')

    for jj in np.linspace(3, 5, 3).astype(int):
        axs[1, 1].plot(results.wave_disc[:, 4], results.WAVEEX[1, :, jj, 2], label=legend_T[jj], linewidth=2, markersize=2)
        axs[1, 1].legend(loc="upper right")
        axs[1, 1].set(xlabel='omega, rad/s', ylabel='90 Deg Exciting Force [N]')
        axs[1, 1].grid(b=True, which='both', axis='both')

