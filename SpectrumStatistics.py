import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle
import sys
import seaborn as sb
import math
from numpy import random
from scipy.signal import find_peaks


def spectrum_response(RAO, s):
    response = np.zeros(shape=(len(RAO[:,0,0]), len(s), len(RAO[0,0,:])))
    for jj in np.linspace(0, len(RAO[:,0,0])-1, len(RAO[:,0,0])).astype(int):
        for ii in np.linspace(0, len(RAO[1, 1, :]) - 1, len(RAO[1, 1, :])).astype(int):
            response[jj, :, ii] = np.square(np.absolute(RAO[jj, :, ii])) * s[:, ii]
    return response


def spectrum_response_statistics_single(Tr, response, f, df_rad):
    m0 = sum(np.absolute(response) * df_rad)
    m1 = sum(np.absolute(response) * f * df_rad)
    m2 = sum(np.absolute(response) * np.square(f) * df_rad)

    Tz = 2 * np.pi * np.sqrt(np.divide(m0, m2))
    Significant_Amplitude= 2 * np.sqrt(m0)
    N_mpm = np.divide(Tr, Tz)
    MPM = np.sqrt(2 * m0 * np.log(N_mpm))

    return Tz, Significant_Amplitude, N_mpm, MPM, m0, m1, m2


def spectrum_response_statistics(Tr, response, f, df_rad):
    m0 = np.zeros(shape=(3, 8))
    m1 = np.zeros(shape=(3, 8))
    m2 = np.zeros(shape=(3, 8))
    Tz = np.zeros(shape=(3, 8))
    Significant_Amplitude = np.zeros(shape=(3, 8))
    N_mpm = np.zeros(shape=(3, 8))
    MPM = np.zeros(shape=(3, 8))
    DOF = [0, 1, 2, 3, 4, 5, 6, 7]

    for jj in np.linspace(0, len(response[:,0,0])-1, len(response[:,0,0])).astype(int):
        for ii in np.linspace(0, len(response[0, 0, :]) - 1, len(response[0, 0, :])).astype(int):
            m0[jj, ii] = sum(np.absolute(response[jj, :, ii]) * df_rad)
            m1[jj, ii] = sum(np.absolute(response[jj, :, ii]) * f * df_rad)
            m2[jj, ii] = sum(np.absolute(response[jj, :, ii]) * np.square(f) * df_rad)

    for jj in np.linspace(0, len(response[:,0,0])-1, len(response[:,0,0])).astype(int):
        for ii in np.linspace(0, len(N_mpm[0,:])-1, len(N_mpm[0,:])).astype(int):
            Tz[jj, DOF[ii]] = 2 * np.pi * np.sqrt(np.divide(m0[jj, DOF[ii]], m2[jj, DOF[ii]]))
            Significant_Amplitude[jj, DOF[ii]] = 2 * np.sqrt(m0[jj, DOF[ii]])
            N_mpm[jj, DOF[ii]] = np.divide(Tr, Tz[jj, DOF[ii]])
            if N_mpm[jj, DOF[ii]] <= 1:
                MPM[jj, DOF[ii]] = Significant_Amplitude[jj, DOF[ii]]
            else:
                MPM[jj, DOF[ii]] = np.sqrt(2 * m0[jj, DOF[ii]] * np.log(N_mpm[jj, DOF[ii]]))

    return Tz, Significant_Amplitude, N_mpm, MPM, m0, m1, m2

