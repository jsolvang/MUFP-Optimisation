import pandas as pd
import numpy as np


class Mass:
    def __init__(self, floater, area, buoy, rho):
        # Turbine masses taken from the NREL 5MW reference turbine
        self.hub = 56780
        self.nacelle = 240000
        self.rotor = 110000
        self._turbine_tower_calculation(floater, rho)

        # Column mass
        self.column = area.cc_column * floater.column_height * rho.steel + area.column * floater.thickness * rho.steel

        # Heave Plate Mass
        self.heave = (area.cc_heave * floater.heave_height * rho.steel) \
                     + (area.heave_top * floater.thickness * rho.steel) \
                     + (area.heave * floater.thickness * rho.steel)

        # Total un-ballasted mass of platform
        float = 3 * (self.column + self.heave)
        self.turbines = 2 * (self.hub + self.nacelle + self.rotor + self.tower)
        total = self.turbines + float

        # Calculating required ballasting mass
        self.ballast_total = buoy.total - total
        self.floater = float + self.ballast_total
        self.total = self.turbines + self.floater

    def _turbine_tower_calculation(self, floater, rho):
        # Predetermined turbine parameters - Change these when structural components are taken into acccount
        t_bot = 0.04
        t_top = 0.02
        r_bot_outer = 3
        r_top_outer = r_bot_outer - t_bot
        r_out = 1.935
        r_in = r_out - t_top

        # Pythagoras to find the tower length
        tower_length = np.sqrt(
            np.square(np.divide(floater.hub_space / 2 - floater.x_space, 2)) + np.square(floater.hub_height))

        # Tower area found by calculating the area of a tapered column. Subtracting the "in" (internal) volume.
        tower_volume_out = np.divide(np.pi, 3) * tower_length * \
                           (np.square(r_bot_outer) + np.square(r_out) + r_bot_outer * r_out)
        tower_volume_in = np.divide(np.pi, 3) * tower_length * \
                          (np.square(r_top_outer) + np.square(r_in) + r_top_outer * r_in)
        self.tower = (tower_volume_out - tower_volume_in) * rho.steel
