from matplotlib import cm
import numpy as np
import ccmaps


class VariableValues:
    def __init__(self, contourf_levels=None, contour_levels=None, ticks=None,
                 cmap=cm.viridis, extend="max"):
        self.contourf_levels = contourf_levels
        self.contour_levels = contour_levels
        self.ticks = ticks
        self.cmap = cmap
        self.extend = extend

variable_cmaps = {
    "CAPE": VariableValues(
        contourf_levels=np.arange(0, 1001, 50),
        ticks=np.arange(0, 1001, 200),
        cmap=ccmaps.get_cmap("severe"),
    ),
    "CIN": VariableValues(
        contourf_levels=np.concatenate([
            np.arange(-400, -200, 50), np.arange(-200, -100, 25),
            np.arange(-100, -50, 12.5), np.arange(-50, 0.1, 6.25)]),
        ticks=[-400, -200, -100, -50, -25, 0],
        cmap=ccmaps.get_cmap("severe_r"),
        extend="min",
    ),
    "CLOUDCOVER": VariableValues(
        contourf_levels=np.arange(0, 101, 5),
        ticks=np.arange(0, 101, 20),
        cmap=cm.Greys_r,
        extend="both",
    ),
    "CTT": VariableValues(
        contourf_levels=np.concatenate(
            [np.arange(-90, -20, 1), np.arange(-20, 41, 2)]),
        ticks=[-90, -80, -70, -60, -50, -40, -30, -20, 0, 20, 40],
        cmap=ccmaps.get_cmap("ir"),
        extend="both",
    ),
    "DBZ": VariableValues(
        contourf_levels=np.arange(17, 96, 1),
        ticks=np.arange(20, 91, 10),
        cmap=ccmaps.get_cmap("radarscope_trunc"),
        extend="both",
    ),
    "GPH250": VariableValues(
        contour_levels=np.arange(954, 1150, 6),
    ),
    "GPH500": VariableValues(
        contour_levels=np.arange(500, 601, 3)
    ),
    "GPH700": VariableValues(
        contour_levels=np.arange(270, 331, 3)
    ),
    "GPH850": VariableValues(
        contour_levels=np.arange(120, 181, 3)
    ),
    "GPH925": VariableValues(
        contour_levels=np.arange(60, 91, 3)
    ),
    "IVT": VariableValues(
        contourf_levels=np.arange(250, 1501, 50),
        ticks=np.arange(250, 1501, 250),
        cmap=cm.YlOrRd,
    ),
    "LR": VariableValues(
        contourf_levels=np.concatenate(
            [np.arange(0, 6, 1.5), np.arange(6, 10.1, 0.25)]),
        ticks=[0, 6, 7, 8, 9, 10],
        cmap=ccmaps.get_cmap("severe"),
    ),
    "OLR": VariableValues(
        contourf_levels=np.arange(150, 301, 5),
        ticks=np.arange(150, 301, 30),
        cmap=cm.Greys,
        extend="both",
    ),
    "OMEGA": VariableValues(
        contourf_levels=np.concatenate(
            [np.arange(-50, 0, 5), np.arange(5, 51, 5)]),
        ticks=np.arange(-50, 51, 25),
        cmap=cm.PiYG_r,
        extend="both",
    ),
    "PWAT": VariableValues(
        contourf_levels=np.arange(0, 61, 2),
        ticks=np.arange(0, 61, 10),
        cmap=ccmaps.get_cmap("tpw"),
        extend="both",
    ),
    "QPF": VariableValues(
        contourf_levels=np.concatenate([
            np.arange(0, 3, 0.5), np.arange(3, 10, 1), np.arange(10, 25, 1.5),
            np.arange(25, 50, 2.5), np.arange(50, 100, 5),
            np.arange(100, 200, 10), np.arange(200, 401, 20)]),
        ticks=[0, 3,10, 25, 50, 100, 200, 400],
        cmap=ccmaps.get_cmap("qpf"),
    ),
    "RH": VariableValues(
        contourf_levels=np.arange(0, 101, 5),
        ticks=np.arange(0, 101, 20),
        cmap=cm.BrBG,
    ),
    "SLP": VariableValues(
        contour_levels=np.arange(800, 1200, 2),
    ),
    "SNOW": VariableValues(
        contourf_levels=np.concatenate([
            np.arange(0, 3, 0.5), np.arange(3, 10, 1), np.arange(10, 25, 1.5),
            np.arange(25, 50, 2.5), np.arange(50, 100, 5),
            np.arange(100, 200, 10), np.arange(200, 401, 20)]),
        ticks=[0, 3,10, 25, 50, 100, 200, 400],
        cmap=ccmaps.get_cmap("snow"),
    ),
    "SWDOWN": VariableValues(
        contourf_levels=np.arange(0, 1201, 50),
        ticks=np.arange(0, 1201, 200),
        cmap=cm.magma,
    ),
    "TMP": VariableValues(
        contourf_levels=np.arange(-20, 40.1, 0.5),
        ticks=np.arange(-20, 40.1, 10),
        cmap=ccmaps.get_cmap("temperature_fine"),
        extend="both",
    ),
    "VORT": VariableValues(
        contourf_levels=np.arange(0, 51, 2.5),
        ticks=np.arange(0, 51, 10),
        cmap=cm.hot_r,
        extend="both",
    ),
    "W": VariableValues(
        contourf_levels=np.concatenate(
            [np.arange(-1, 0, 0.1), np.arange(0.1, 1.1, 0.1)]),
        ticks=np.arange(-1, 1.1, 0.5),
        cmap=cm.PiYG,
        extend="both",
    ),
    "WIND10M": VariableValues(
        contourf_levels=np.arange(5, 26, 1),
        ticks=np.arange(5, 26, 5),
        cmap=ccmaps.get_cmap("wind"),
    ),
    "WIND250": VariableValues(
        contourf_levels=np.arange(25, 106, 4),
        ticks=np.arange(25, 106, 20),
        cmap=ccmaps.get_cmap("wind"),
    ),
    "WIND500": VariableValues(
        contourf_levels=np.arange(10, 51, 2),
        ticks=np.arange(10, 51, 10),
        cmap=ccmaps.get_cmap("wind"),
    ),
    "WIND700": VariableValues(
        contourf_levels=np.arange(10, 51, 2),
        ticks=np.arange(10, 51, 10),
        cmap=ccmaps.get_cmap("wind"),
    ),
    "WIND850": VariableValues(
        contourf_levels=np.arange(10, 51, 2),
        ticks=np.arange(10, 51, 10),
        cmap=ccmaps.get_cmap("wind"),
    ),
    "WIND925": VariableValues(
        contourf_levels=np.arange(10, 51, 2),
        ticks=np.arange(10, 51, 10),
        cmap=ccmaps.get_cmap("wind"),
    ),
}
