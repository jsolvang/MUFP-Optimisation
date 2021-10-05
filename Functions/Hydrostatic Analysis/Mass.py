import pandas as pd
import numpy as np


class Mass:
    def __init__(self, floater, area, buoy, rho):
        # Turbine masses taken from the NREL 5MW reference turbine
        self.hub = 56780
        self.nacelle = 240000
        self.rotor = 110000
        self._turbine_tower_calculation(floater, rho)
        self._calculate_bracing_mass(floater, rho)

        # Column mass
        self.column = area.cc_column * floater.column_height * rho.steel + area.column * floater.thickness * rho.steel

        # Heave Plate Mass
        self.heave = (area.cc_heave * floater.heave_height * rho.steel) \
                     + (area.heave_top * floater.thickness * rho.steel) \
                     + (area.heave * floater.thickness * rho.steel)

        # Total un-ballasted mass of platform
        float = 3 * (self.column + self.heave)
        self.turbines = 2 * (self.hub + self.nacelle + self.rotor + self.tower)
        self.total_unballasted = self.turbines + float + self.bracing

        # Calculating required ballasting mass
        self.ballast_total = buoy.total - self.total_unballasted
        if self.ballast_total < 0:
            self.ballast_total = 1
        self.floater = float + self.ballast_total
        self.total = self.turbines + self.floater + self.bracing
        self.steel_mass = 3 * (self.column + self.heave) + self.turbines + self.bracing

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
            np.square(np.divide(floater.hub_space / 2 - floater.y_space, 2)) + np.square(floater.hub_height))

        # Tower area found by calculating the area of a tapered column. Subtracting the "in" (internal) volume.
        tower_volume_out = np.divide(np.pi, 3) * tower_length * \
                           (np.square(r_bot_outer) + np.square(r_out) + r_bot_outer * r_out)
        tower_volume_in = np.divide(np.pi, 3) * tower_length * \
                          (np.square(r_top_outer) + np.square(r_in) + r_top_outer * r_in)
        self.tower = (tower_volume_out - tower_volume_in) * rho.steel

    def _calculate_bracing_mass(self, floater, rho):
        # Calculating length of bracing between columns
        l_xy = np.sqrt(floater.x_space ** 2 + (floater.y_space/2) ** 2) - floater.dia_column
        l_y = floater.y_space - floater.dia_column
        h = 14

        # Calculating number of bays necessary to obtain 30 deg weld angle
        nxy_bays, lxy_bay = _calc_nbays(l_xy, h)
        ny_bays, ly_bay = _calc_nbays(l_y, h)

        # Calculating length of diagonal bracing members
        lxy_diag = np.sqrt(lxy_bay ** 2 + h ** 2)
        ly_diag = np.sqrt(ly_bay ** 2 + h ** 2)

        # Calculating diameter based on slenderness rule of thumb
        dxy_diag = 0.029 * lxy_diag
        dy_diag = 0.029 * ly_diag
        dxy_hori = 0.023 * l_xy
        dy_hori = 0.023 * l_y

        # Calculating thickness of each member
        txy_hori = dxy_hori / 30
        ty_hori = dy_hori / 30
        txy_diag = dxy_diag / 45
        ty_diag = dy_diag / 45

        self.xy_bracing = 2 * (1 / 4) * np.pi * (dxy_hori ** 2 - (dxy_hori - 2 * txy_hori) ** 2) * rho.steel * l_xy \
                        + nxy_bays * 2 * (1 / 4) * np.pi * (dxy_diag ** 2 - (dxy_diag - 2 * txy_diag) ** 2) * rho.steel * lxy_diag \

        self.y_bracing = 2 * (1 / 4) * np.pi * (dy_hori ** 2 - (dy_hori - 2 * ty_hori) ** 2) * rho.steel * l_y \
                       + ny_bays * 2 * (1 / 4) * np.pi * (dy_diag ** 2 - (dy_diag - 2 * ty_diag) ** 2) * rho.steel * ly_diag

        self.bracing = 2*self.xy_bracing + self.y_bracing


def _calc_nbays(L, h):
    found = 0
    nbays = 1
    while found == 0:
        lbay = (L / 2) / nbays
        theta = np.arctan(h / lbay) * (180 / np.pi)
        if nbays > 10:
            print('Too many bays')
        if theta > 30:
            found == 0
            break
        else:
            nbays += 1
    return nbays, lbay