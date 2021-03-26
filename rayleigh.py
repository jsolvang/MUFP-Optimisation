import numpy as np

def Rayleigh(x,m0):
    f = (x/m0)*np.exp(-np.divide(x**2,2*m0))
    return f