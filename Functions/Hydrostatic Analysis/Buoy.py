import pandas as pd
import numpy as np


class Buoy:
    def __init__(self, floater, area, rho):
        # Calculating the buoyancy of the columns and heave plates
        self.heave = area.heave * floater.heave_height * rho.water
        self.column = area.column * (floater.draft - floater.heave_height) * rho.water

        # Displaced Mass
        self.total = 3 * (self.heave + self.column)
        self.displaced_volume = self.total / rho.water

