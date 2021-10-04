import struct
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import signal
from scipy import interpolate


def WindForceSpectrum():

    path1 = r"C:\Users\Joar\Documents\1_Education\NTNU\OneDrive - NTNU\Thesis\Modelling\SIMA\FixedRotor5MW\conditionSet\296-20210425142045\conditionSet_1"
    path2 = r"C:\Users\Joar\Documents\1_Education\NTNU\OneDrive - NTNU\Thesis\Modelling\SIMA\FixedRotor5MW\conditionSet\296-20210425142045\conditionSet_2"

    simtime = 600.
    dt = 0.1
    transient = 100

    # ---------------------------------------------Wind Turbine One ---------------------------------------------#
    with open(path1 + r'\sima_witurb.bin', mode='rb') as file:
        CC = file.read()

    cols = len(np.asarray(struct.unpack("f" * (len(CC) // 4), CC))) / (simtime / dt)
    CC = np.transpose(
        np.reshape(np.asarray(struct.unpack("f" * (len(CC) // 4), CC)), (int(cols), int(simtime / dt)), order='F'))

    time = CC[:, 1]
    omega = CC[:, 2] * np.pi / 180  # convert from deg/s to rad/s
    HubWindX1 = CC[:, 7]
    AeroForceX1 = CC[:, 10] * 1e3

    # ---------------------------------------------Wind Turbine Two ---------------------------------------------#

    with open(path2 + r'\sima_witurb.bin', mode='rb') as file:
        CC = file.read()

    cols = len(np.asarray(struct.unpack("f" * (len(CC) // 4), CC))) / (simtime / dt)
    CC = np.transpose(
        np.reshape(np.asarray(struct.unpack("f" * (len(CC) // 4), CC)), (int(cols), int(simtime / dt)), order='F'))

    time = CC[:, 1]
    omega = CC[:, 2] * np.pi / 180  # convert from deg/s to rad/s
    HubWindX2 = CC[:, 7]
    AeroForceX2 = CC[:, 10] * 1e3

    # ---------------------------------------------Aerodynamic  Moment  ---------------------------------------------#

    mom_wind = AeroForceX1[int(transient / dt):] * 71.05 + AeroForceX2[int(transient / dt):] * -71.05

    # ----------------------------------------- Power Spectral Density-----------------------------------------------#

    f_xx, sx_wind = signal.welch(AeroForceX1[int(transient / dt):] + AeroForceX2[int(transient / dt):], 1 / dt)
    f_yy, sy_wind = signal.welch((AeroForceX1[int(transient / dt):] + AeroForceX2[int(transient / dt):]) * (3 * np.pi / 180), 1 / dt)
    f_mz, smz_wind = signal.welch(mom_wind, 1 / dt)

    f_xx = f_xx * 2 * np.pi
    sx_wind = sx_wind / (2 * np.pi)

    f_yy = f_yy * 2 * np.pi
    sy_wind = sy_wind / (2 * np.pi)

    f_mz = f_mz * 2 * np.pi
    smz_wind = smz_wind / (2 * np.pi)

    return f_xx, sx_wind, sy_wind, smz_wind















