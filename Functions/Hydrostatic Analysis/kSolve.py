from scipy.optimize import fsolve
import pandas as pd
import numpy as np


def ksolve(omega, h, g):
    k = fsolve(dispersion_equation, 1, args=(g,h,omega), maxfev=1000)
    k = k[0]
    L = 2*np.pi/k
    return k, L


def dispersion_equation(k, g, h, omega):
    fun = g*k*np.tanh(k*h)-omega**2
    return fun
