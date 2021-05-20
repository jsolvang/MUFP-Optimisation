import numpy as np
import matplotlib.pyplot as plt


def calculate_column_raos(RAO, matrix, column):
    column_rao = RAO

    if column == 'front':
        column_rao[:, :, 2] = RAO[:, :, 2] + RAO[:, :, 4] * matrix.distances['Front Centre']
    elif column == 'left':
        column_rao[:, :, 2] = RAO[:, :, 2] - RAO[:, :, 3] * matrix.distances['Left Centre'] - RAO[:, :, 4] * matrix.distances['Back Centre']
    elif column == 'right':
        column_rao[:, :, 2] = RAO[:, :, 2] + RAO[:, :, 3] * matrix.distances['Right Centre'] - RAO[:, :, 4] * matrix.distances['Back Centre']
    else:
        print('Invalid Column Input')

    return column_rao


def calulate_RAOs(wave_disc, mass, ADDEDMASS, DAMPING, stiffness, WAVEEX, plot, floater, coord):
    Y = np.zeros(shape=(len(wave_disc[:, 4]), 6, 6)) + 0j
    Y_inv = np.zeros(shape=(len(wave_disc[:, 4]), 6, 6)) + 0j
    H = np.zeros(shape=(len(wave_disc[:, 4]), 6)) + 0j
    X = np.zeros(shape=(3, len(wave_disc[:, 4]), 8)) + 0j

    omega_index = np.linspace(0, len(wave_disc[:, 4]) - 1, len(wave_disc[:, 4])).astype(int)

    for jj in omega_index:
        Y[jj, :, :] = -np.square(wave_disc[jj, 4]) * (ADDEDMASS[jj, :, :] + mass) + 1j * wave_disc[
            jj, 4] * DAMPING[jj, :, :] + stiffness

    for jj in omega_index:
        Y_inv[jj, :, :] = np.linalg.inv(Y[jj, :, :])

    for jj in omega_index:
        H[jj, :] = WAVEEX[0, jj, :, 0] + 1j * WAVEEX[0, jj, :, 1]

    for jj in omega_index:
        X[0, jj, 0:6] = np.dot(np.linalg.inv(Y[jj, :, :]), H[jj, :])
        X[1, jj, 0:6] = np.dot(np.linalg.inv(Y[jj, :, :]), H[jj, :]) * wave_disc[jj, 4]
        X[2, jj, 0:6] = np.dot(np.linalg.inv(Y[jj, :, :]), H[jj, :]) * np.square(wave_disc[jj, 4])

    X[0, :, 6] = X[0, :, 4] * np.sqrt(floater.hub_height ** 2 + (floater.x_space - coord.COM[0]) ** 2)
    X[0, :, 7] = X[0, :, 3] * np.sqrt(floater.hub_height ** 2 + (floater.y_space / 2) ** 2)

    X[1, :, 6] = X[1, :, 4] * np.sqrt(floater.hub_height ** 2 + (floater.x_space - coord.COM[0]) ** 2)
    X[1, :, 7] = X[1, :, 3] * np.sqrt(floater.hub_height ** 2 + (floater.y_space / 2) ** 2)

    X[2, :, 6] = X[2, :, 4] * np.sqrt(floater.hub_height ** 2 + (floater.x_space - coord.COM[0]) ** 2)
    X[2, :, 7] = X[2, :, 3] * np.sqrt(floater.hub_height ** 2 + (floater.y_space / 2) ** 2)

    if plot == 1:

        Y_lab = ["RAO Surge, [m/m]",
                 "RAO Sway, [m/m]",
                 "RAO Heave, [m/m]",
                 "RAO Roll, [rad/m]",
                 "RAO Pitch, [rad/m]",
                 "RAO Yaw, [rad/m]",
                 "RAO Nacelle For-aft [m/m]",
                 "RAO Nacelle Side-Side [m/m]"]

        plt.rcParams["figure.figsize"] = (15, 40)

        fig, axs = plt.subplots(8)
        for jj in np.linspace(0, 7, 8).astype(int):
            axs[jj].plot(wave_disc[:, 4], np.real(X[0, :, jj]), '-r', label="Real", linewidth=2,
                         markersize=5)
            axs[jj].plot(wave_disc[:, 4], np.imag(X[1, :, jj]), '--b', label="Imaginary", linewidth=2,
                         markersize=5)
            axs[jj].plot(wave_disc[:, 4], np.absolute(X[2, :, jj]), ':g', label="Absolute", linewidth=2,
                         markersize=5)
            axs[jj].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
            axs[jj].set(xlabel='omega, rad/s', ylabel=Y_lab[jj])
            axs[jj].set(xlim=[0, 2])
            axs[jj].legend(loc="upper right")
            axs[jj].grid(b=True, which='both', axis='both')

        plt.tight_layout()
    plt.rcParams["figure.figsize"] = (15, 5)
    return X, Y, H, Y_inv
