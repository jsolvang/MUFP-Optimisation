import pandas as pd
import numpy as np


class Environment:
    # The class defines environmental parameters
    def __init__(self):
        self.g = 9.80665 #m/s^2
        self.h = 300 # Still Water Level
