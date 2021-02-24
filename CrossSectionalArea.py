import pandas as pd
import numpy as np


class CrossSectionalArea:
    def __init__(self, floater):
        # Calculating the cross-sectional area of the columns and heave plates
        self.column = np.divide((np.pi * np.square(floater.dia_column)), 4)
        self.heave = np.divide((np.pi * np.square(floater.dia_heave)), 4)

        # Area of steel between heave plate and column
        self.heave_top = self.heave - self.column

        # Cross sectional area of steel for column and heave plate
        self.cc_column = self.column - np.divide(
            (np.pi * np.square(floater.dia_column - 2 * floater.thickness)), 4)
        self.cc_heave = self.heave - np.divide(
            (np.pi * np.square(floater.dia_heave - np.multiply(2, floater.thickness))), 4)
