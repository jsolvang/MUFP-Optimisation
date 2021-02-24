import pandas as pd
import numpy as np


class CrossSectionalArea:
    def __init__(self, floater):
        self.column = np.divide((np.pi * np.square(floater.dia_column)), 4)
        self.heave = np.divide((np.pi * np.square(floater.dia_heave)), 4)
        self.heave_top = self.heave - self.column

        self.cc_column = self.column - np.divide(
            (np.pi * np.square(floater.dia_column - 2 * floater.thickness)), 4)
        self.cc_heave = self.heave - np.divide(
            (np.pi * np.square(floater.dia_heave - np.multiply(2, floater.thickness))), 4)
