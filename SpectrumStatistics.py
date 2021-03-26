import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle
import sys
import seaborn as sb
import math
from numpy import random
from scipy.signal import find_peaks

def spectrum_response(s, Tr, RAO, sim, f, df_rad):
    RAO_interp = np.zeros(shape=(3, len(f), 6))
    response = np.zeros(shape=(3, len(f), 6))
    m0 = np.zeros(shape=(3, 6))
    m1 = np.zeros(shape=(3, 6))
    m2 = np.zeros(shape=(3, 6))
    Tz = np.zeros(shape=(3, 6))
    Significant_Amplitude = np.zeros(shape=(3, 6))
    N_mpm = np.zeros(shape=(3, 6))
    MPM = np.zeros(shape=(3, 6))
    DOF = [0, 2, 4]


    for jj in np.linspace(0, 2, 3).astype(int):
        for ii in np.linspace(0, len(RAO[1, 1, :]) - 1, len(RAO[1, 1, :])).astype(int):
            response[jj, :, ii] = np.square(RAO_interp[jj, :, ii]) * s

    for jj in np.linspace(0, 2, 3).astype(int):
        for ii in np.linspace(0, len(RAO[1, 1, :]) - 1, len(RAO[1, 1, :])).astype(int):
            m0[jj, ii] = sum(response[jj, :, ii] * df_rad)
            m1[jj, ii] = sum(response[jj, :, ii] * f * df_rad)
            m2[jj, ii] = sum(response[jj, :, ii] * np.square(f) * df_rad)

    for jj in np.linspace(0, 2, 3).astype(int):
        for ii in np.linspace(0, 2, 3).astype(int):
            Tz[jj, DOF[ii]] = 2 * np.pi * np.sqrt(np.divide(m0[jj, DOF[ii]], m2[jj, DOF[ii]]))
            Significant_Amplitude[jj, DOF[ii]] = 4 * np.sqrt(m0[jj, DOF[ii]])
            N_mpm[jj, DOF[ii]] = np.divide(Tr, Tz[jj, DOF[ii]])
            MPM[jj, DOF[ii]] = np.sqrt(2 * m0[jj, DOF[ii]] * np.log(N_mpm[jj, DOF[ii]]))

    return Tz, Significant_Amplitude, N_mpm, MPM, response