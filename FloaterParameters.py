import pandas as pd
import numpy as np


class FloaterParameters:
    def __init__(self, x, y, dia):
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
        self.draft = 14
        self.height = 25
        self.thickness = 0.07


        self.free_board = self.height - self.draft
        self.column_height = self.height - self.heave_height


