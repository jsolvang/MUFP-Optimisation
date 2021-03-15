import numpy as np
import matplotlib.pyplot as plt

def calulate_RAOs(wave_disc, mass, ADDEDMASS, DAMPING, stiffness, WAVEEX, MOTIONS, plot):
    Y = np.zeros(shape=(len(wave_disc[:, 4]), 6, 6)) * 1j
    H = np.zeros(shape=(len(wave_disc[:, 4]), 6)) * 1j
    X = np.zeros(shape=(3, len(wave_disc[:, 4]), 6)) * 1j

    Y_lab = ["RAO Surge, [m/m]",
             "RAO Sway, [m/m]",
             "RAO Heave, [m/m]",
             "RAO Roll, [rad/m]",
             "RAO Pitch, [rad/m]",
             "RAO Yaw, [rad/m]"]

    omega_index = np.linspace(0, len(wave_disc[:, 4]) - 1, len(wave_disc[:, 4])).astype(int)

    for jj in omega_index:
        Y[jj, :, :] = -np.square(wave_disc[jj, 4]) * (ADDEDMASS[jj, :, :] + mass) + 1j * wave_disc[
            jj, 4] * DAMPING[jj, :, :] + stiffness

    for jj in omega_index:
        H[jj, :] = WAVEEX[0, jj, :, 0] + 1j * WAVEEX[0, jj, :, 1]

    for jj in omega_index:
        X[0, jj, :] = np.dot(np.linalg.inv(Y[jj, :, :]), H[jj, :])
        X[1, :, :] = X[0, :, :] * wave_disc[jj, 4]
        X[2, :, :] = X[0, :, :] * np.square(wave_disc[jj, 4])

    if plot == 1:
        plt.rcParams["figure.figsize"] = (10, 20)
        fig, axs = plt.subplots(6)
        for jj in np.linspace(0, 5, 6).astype(int):
            axs[jj].plot(wave_disc[:, 4], np.absolute(X[0, :, jj]), '-b', label="Transfer Function", linewidth=2,
                         markersize=5)
            axs[jj].plot(wave_disc[:, 4], MOTIONS[0, :, jj, 2], '--r', label="HydroD Output", linewidth=2,
                         markersize=5)
            axs[jj].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
            axs[jj].set(xlabel='omega, rad/s', ylabel=Y_lab[jj], title="Displacement")
            axs[jj].legend(loc="upper right")
            axs[jj].grid(b=True, which='both', axis='both')
        plt.tight_layout()

        fig, axs = plt.subplots(6)
        for jj in np.linspace(0, 5, 6).astype(int):
            axs[jj].plot(wave_disc[:, 4], np.absolute(X[1, :, jj]), '-b', label="Transfer Function", linewidth=2,
                         markersize=5)
            axs[jj].set(xlabel='omega, rad/s', ylabel=Y_lab[jj], title="Velcoity")
            axs[jj].legend(loc="upper right")
            axs[jj].grid(b=True, which='both', axis='both')
        plt.tight_layout()

        fig, axs = plt.subplots(6)
        for jj in np.linspace(0, 5, 6).astype(int):
            axs[jj].plot(wave_disc[:, 4], np.absolute(X[2, :, jj]), '-b', label="Transfer Function", linewidth=2,
                         markersize=5)
            axs[jj].set(xlabel='omega, rad/s', ylabel=Y_lab[jj], title="Acceleration")
            axs[jj].legend(loc="upper right")
            axs[jj].grid(b=True, which='both', axis='both')
        plt.tight_layout()

    return X
