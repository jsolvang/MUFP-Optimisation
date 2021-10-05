import os
import matplotlib.pyplot as plt
import numpy as np


def sub_plots(x, y, y_lab, plot, y_lim, sci):
    plt.rcParams.update({'font.size': 20})
    plt.rcParams["figure.figsize"] = (15, 5 * len(plot))
    fig, axs = plt.subplots(len(plot))
    for ii,jj in enumerate(plot):
        axs[ii].plot(x, np.real(y[:, jj]), '-r', label='Real')
        axs[ii].plot(x, np.imag(y[:, jj]), '--b', label='Imaginary')
        axs[ii].plot(x, np.abs(y[:, jj]), ':g', label='Absolute', linewidth=3)
        if sci == 1:
            axs[ii].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        axs[ii].set(xlabel='Frequency, [rad/s]', ylabel=y_lab[jj])
        axs[ii].legend(loc="upper right")
        axs[ii].set(xlim=(0, 5))
        if y_lim != [0,0]:
            axs[ii].set(ylim=y_lim)
    plt.tight_layout()
    plt.rcParams["figure.figsize"] = (15, 5)

def sub_plots_spectrums(x, y, y_lab, plot, y_lim, sci):
    plt.rcParams.update({'font.size': 20})
    plt.rcParams["figure.figsize"] = (15, 5 * len(plot))
    fig, axs = plt.subplots(len(plot))
    for ii,jj in enumerate(plot):
        axs[ii].plot(x, y[:, jj], '-b', linewidth=3)
        if sci == 1:
            axs[ii].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        axs[ii].set(xlabel='Frequency, [rad/s]', ylabel=y_lab[jj])
        axs[ii].set(xlim=(0, 2))
        if y_lim != [0,0]:
            axs[ii].set(ylim=y_lim)
    plt.tight_layout()
    plt.rcParams["figure.figsize"] = (15, 5)