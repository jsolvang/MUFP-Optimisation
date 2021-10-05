import pandas as pd
import numpy as np


class FloaterParameters:
        # This class defines the floater dimensions
        # Input in meters
    def __init__(self, x, y, dia, draft):
        # Free Variables
        self.x_space = x
        self.y_space = y

        # Currently fixed variables
        self.dia_column = dia
        self.dia_heave = 20
        self.heave_height = 2

        # Fixed variables
        self.hub_height = 85
        self.hub_space = 142.1
        self.draft = draft
        self.free_board = 11
        self.height = self.draft + self.free_board
        self.thickness = 0.07
        self.column_height = self.height - self.heave_height


