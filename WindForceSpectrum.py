import struct
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import signal
from scipy import interpolate


def WindForceSpectrum():

    path = r"C:\Users\Joar\Documents\1_Education\NTNU\OneDrive - NTNU\Thesis\Modelling\SIMA\FixedRotor5MW\Initial\147-20210416154156"

    # ---------------------------------------------Wind Turbine Results---------------------------------------------#
    with open(path + r'\sima_witurb.bin', mode='rb') as file:
        CC = file.read()

    # Total simulation time and storage time step
    simtime = 600.
    dt = 0.1
    transient = 100

    # Element number in storage list ()
    n_beam = 1

    cols = len(np.asarray(struct.unpack("f" * (len(CC) // 4), CC))) / (simtime / dt)
    CC = np.transpose(
        np.reshape(np.asarray(struct.unpack("f" * (len(CC) // 4), CC)), (int(cols), int(simtime / dt)), order='F'))

    time = CC[:, 1]
    omega = CC[:, 2] * np.pi / 180  # convert from deg/s to rad/s
    genTq = CC[:, 4]
    genPwr = CC[:, 5]
    azimuth = CC[:, 6]
    HubWindX = CC[:, 7]
    HubWindY = CC[:, 8]
    HubWindZ = CC[:, 9]
    AeroForceX = CC[:, 10]*1e3
    AeroMomX = CC[:, 13]*1e3
    Bl1Pitch = CC[:, 16]
    Bl2Pitch = CC[:, 17]
    Bl3Pitch = CC[:, 18]

    f_xx, s_wind = signal.welch(AeroForceX[int(transient / dt):] * 2, 1 / dt)

    return f_xx, s_wind, time, AeroForceX















