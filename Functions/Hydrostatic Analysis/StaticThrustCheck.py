import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle
import sys
import seaborn as sb
import math
from numpy import random
from scipy.signal import find_peaks

def StaticThrustCheck(stiffness_matrix,floater):
    # Easy deflection check against maximum static loading
    rated_thrust = 820000
    deflection = (2*rated_thrust*floater.hub_height)/stiffness_matrix[4, 4]
    return deflection*(180/np.pi)