import config
from variable_config import variable_cmaps
from scipy.ndimage import gaussian_filter
import matplotlib.colors as mpcolors
import numpy as np
from metpy.units import units
import wrf_calc


class ProductDefault:
    def __init__(self, name,
            contourf_var, contourf_kwargs={}, cmap_label="",
            barb_u=False, barb_v=False, barb_kwargs={},
            barb_u_2=False, barb_v_2=False, barb_kwargs_2={},
            barb_u_3=False, barb_v_3=False, barb_kwargs_3={},
            quiver_u=False, quiver_v=False, quiver_kwargs={},
            contour_var=False, contour_kwargs={},
            mask_var=False,
            title="",
            shapefile_kwargs={},
            filename=False,
            stations=True, station_kwargs={},
            hilo_var=False, hilo_kwargs={},
            scales=["regional"]):
        self.name = name
        self.contourf_var = contourf_var
        self.contourf_kwargs = {
            "levels": variable_cmaps[name].contourf_levels,
            "ticks": variable_cmaps[name].ticks,
            "cmap": variable_cmaps[name].cmap,
            "extend": variable_cmaps[name].extend,
            "norm": mpcolors.BoundaryNorm(
                variable_cmaps[name].contourf_levels,
                variable_cmaps[name].cmap.N, clip=True),
            **contourf_kwargs}
        self.cmap_label = cmap_label
        self.mask_var = mask_var
        self.barb_u = barb_u
        self.barb_v = barb_v
        self.barb_kwargs = {
            "length": 5,
            "linewidth": 0.5,
            "pivot": "middle",
            **barb_kwargs}
        self.barb_u_2 = barb_u_2
        self.barb_v_2 = barb_v_2
        self.barb_kwargs_2 = {
            "length": 5,
            "linewidth": 0.5,
            "pivot": "middle",
            **barb_kwargs_2}
        self.barb_u_3 = barb_u_3
        self.barb_v_3 = barb_v_3
        self.barb_kwargs_3 = {
            "length": 5,
            "linewidth": 0.5,
            "pivot": "middle",
            **barb_kwargs_3}
        self.quiver_u = quiver_u
        self.quiver_v = quiver_v
        self.quiver_kwargs = {
            "pivot": "middle",
            **quiver_kwargs}
        self.contour_var = contour_var
        self.contour_kwargs = {
            "colors": "black",
            "linewidths": 1,
            **contour_kwargs} 
        self.title = title
        self.shapefile_kwargs = {
            "edgecolor": "black",
            **shapefile_kwargs}
        self.filename = filename if filename else name
        self.stations = stations
        self.station_kwargs = {
            "decimal": 0,
            **station_kwargs}
        self.hilo_var = hilo_var
        self.hilo_kwargs = {
            "gaussian_sigma": config.mpas_slp_smooth_sigma,
            **hilo_kwargs}
        self.scales = scales
        

mpas_products = {
    "CAPE": ProductDefault(
        name="CAPE",
        contourf_var=lambda ds: ds["cape"],
        contourf_kwargs={"alpha": 0.5},
        cmap_label="[J kg$^{-1}$]",
        barb_u=lambda ds: ds["u10"].where(ds["cape"] > 50),
        barb_v=lambda ds: ds["v10"].where(ds["cape"] > 50),
        barb_kwargs={"pivot": "tip"},
        barb_u_2=lambda ds: ds["uzonal_850hPa"].where(ds["cape"] > 50),
        barb_v_2=lambda ds: ds["umeridional_850hPa"].where(ds["cape"] > 50),
        barb_kwargs_2={"barbcolor": "red", "pivot": "tip"},
        barb_u_3=lambda ds: ds["uzonal_500hPa"].where(ds["cape"] > 50),
        barb_v_3=lambda ds: ds["umeridional_500hPa"].where(ds["cape"] > 50),
        barb_kwargs_3={"barbcolor": "green", "pivot": "tip"},
        title="CAPE [J kg$^{-1}$] | 10m/850mb/500mb Wind (black/red/green)",
        stations=False,
    ),
    "CIN": ProductDefault(
        name="CIN",
        contourf_var=lambda ds: -ds["cin"],
        cmap_label="[J kg$^{-1}$]",
        title="CIN [J kg$^{-1}$]",
        stations=False,
    ),
    "CLOUDHIGH": ProductDefault(
        name="CLOUDCOVER",
        contourf_var=lambda ds: ds["cldfrac_high_UPP"] * 100,
        cmap_label="[%]",
        title="High Cloud Cover [%]",
        shapefile_kwargs={"color": "lime"},
        stations=False,
        filename="cloudhigh",
    ),
    "CLOUDLOW": ProductDefault(
        name="CLOUDCOVER",
        contourf_var=lambda ds: ds["cldfrac_low_UPP"] * 100,
        cmap_label="[%]",
        title="Low Cloud Cover [%]",
        shapefile_kwargs={"color": "lime"},
        stations=False,
        filename="cloudlow",
    ),
    "CLOUDMID": ProductDefault(
        name="CLOUDCOVER",
        contourf_var=lambda ds: ds["cldfrac_mid_UPP"] * 100,
        cmap_label="[%]",
        title="Mid Cloud Cover [%]",
        shapefile_kwargs={"color": "lime"},
        stations=False,
        filename="cloudmid",
    ),
    "CLOUDTOTAL": ProductDefault(
        name="CLOUDCOVER",
        contourf_var=lambda ds: ds["cldfrac_tot_UPP"] * 100,
        cmap_label="[%]",
        title="Total Cloud Cover [%]",
        shapefile_kwargs={"color": "lime"},
        stations=False,
        filename="cloudtotal",
    ),
    "DBZ": ProductDefault(
        name="DBZ",
        contourf_var=lambda ds: ds["refl10cm_max"],
        cmap_label="[dBZ]",
        contour_var=lambda ds: gaussian_filter(
            ds["mslp"] / 100, config.mpas_slp_smooth_sigma),
        contour_kwargs={"levels": variable_cmaps["SLP"].contour_levels},
        barb_u=lambda ds: ds["u10"],
        barb_v=lambda ds: ds["v10"],
        title="Composite Reflectivity [dBZ] | 10m Wind [m s$^{-1}$] | MSLP [mb]",
        stations=False,
        hilo_var=lambda ds: ds["mslp"] / 100,
    ),
    "DBZ1KM": ProductDefault(
        name="DBZ",
        contourf_var=lambda ds: ds["refl10cm_1km"],
        cmap_label="[dBZ]",
        contour_var=lambda ds: gaussian_filter(
            ds["mslp"] / 100, config.mpas_slp_smooth_sigma),
        contour_kwargs={"levels": variable_cmaps["SLP"].contour_levels},
        barb_u=lambda ds: ds["u10"],
        barb_v=lambda ds: ds["v10"],
        title="1km Reflectivity [dBZ] | 10m Wind [m s$^{-1}$] | MSLP [mb]",
        stations=False,
        hilo_var=lambda ds: ds["mslp"] / 100,
        filename="dbz1km",
    ),
    "LR700500": ProductDefault(
        name="LR",
        contourf_var=lambda ds: (
            ds["temperature_500hPa"] - ds["temperature_700hPa"]) /
            (ds["height_700hPa"] - ds["height_500hPa"]) * 1000,
        cmap_label="[°C km$^{-1}$]",
        title="700-500mb Lapse Rate [°C km$^{-1}$]",
        stations=False,
    ),
    "OLR": ProductDefault(
        name="OLR",
        contourf_var=lambda ds: ds["olrtoa"],
        cmap_label="[W m$^{-2}$]",
        title="Top of Atmosphere Outgoing Longwave Radiation [W m$^{-2}$]",
        stations=False,
    ),
    "PWAT": ProductDefault(
        name="PWAT",
        contourf_var=lambda ds: ds["precipw"],
        cmap_label="[kg m$^{-2}$]",
        barb_u=lambda ds: ds["uzonal_850hPa"],
        barb_v=lambda ds: ds["umeridional_850hPa"],
        contour_var=lambda ds: gaussian_filter(
            ds["height_850hPa"].where(ds["surface_pressure"] > 85000) / 10,
            config.mpas_height_smooth_sigma),
        contour_kwargs={"levels": variable_cmaps["GPH850"].contour_levels},
        title="Precipitable Water [kg m$^{-2}$] | 850mb Height [dam] | Wind [m s$^{-1}$]",
        stations=False,
    ),
    "QPF3H": ProductDefault(
        name="QPF",
        contourf_var=lambda ds: ds["raincv"] + ds["rainncv"],
        contourf_kwargs={"norm": None},
        cmap_label="[mm]",
        title="3hr Precipitation [mm]",
        stations=True,
        station_kwargs={"station_condition": lambda x: x > 0.5},
        filename="qpf3h",
    ),
    "QPFT": ProductDefault(
        name="QPF",
        contourf_var=lambda ds: ds["rainc"] + ds["rainnc"],
        contourf_kwargs={"norm": None},
        cmap_label="[mm]",
        title="Total Accumulated Precipitation [mm]",
        stations=True,
        station_kwargs={"station_condition": lambda x: x > 0.5},
        filename="qpft",
    ),
    "RH250": ProductDefault(
        name="RH",
        mask_var=lambda ds: ds["surface_pressure"].where(
            ds["surface_pressure"] < 25000),
        contourf_var=lambda ds: ds["relhum_250hPa"],
        contourf_kwargs={},
        cmap_label="[%]",
        barb_u=lambda ds: ds["uzonal_250hPa"],
        barb_v=lambda ds: ds["umeridional_250hPa"],
        contour_var=lambda ds: gaussian_filter(
            ds["height_250hPa"].where(ds["surface_pressure"] > 25000) / 10,
            config.mpas_height_smooth_sigma),
        contour_kwargs={"levels": variable_cmaps["GPH250"].contour_levels},
        title="250mb Relative Humidity [%] | Height [dam] | Wind [m s$^{-1}$]",
        stations=False,
        filename="rh250",
    ),
    "RH500": ProductDefault(
        name="RH",
        mask_var=lambda ds: ds["surface_pressure"].where(
            ds["surface_pressure"] < 50000),
        contourf_var=lambda ds: ds["relhum_500hPa"],
        contourf_kwargs={},
        cmap_label="[%]",
        barb_u=lambda ds: ds["uzonal_500hPa"],
        barb_v=lambda ds: ds["umeridional_500hPa"],
        contour_var=lambda ds: gaussian_filter(
            ds["height_500hPa"].where(ds["surface_pressure"] > 50000) / 10,
            config.mpas_height_smooth_sigma),
        contour_kwargs={"levels": variable_cmaps["GPH500"].contour_levels},
        title="500mb Relative Humidity [%] | Height [dam] | Wind [m s$^{-1}$]",
        stations=False,
        filename="rh500",
    ),
    "RH700": ProductDefault(
        name="RH",
        mask_var=lambda ds: ds["surface_pressure"].where(
            ds["surface_pressure"] < 70000),
        contourf_var=lambda ds: ds["relhum_700hPa"],
        contourf_kwargs={},
        cmap_label="[%]",
        barb_u=lambda ds: ds["uzonal_700hPa"],
        barb_v=lambda ds: ds["umeridional_700hPa"],
        contour_var=lambda ds: gaussian_filter(
            ds["height_700hPa"].where(ds["surface_pressure"] > 70000) / 10,
            config.mpas_height_smooth_sigma),
        contour_kwargs={"levels": variable_cmaps["GPH700"].contour_levels},
        title="700mb Relative Humidity [%] | Height [dam] | Wind [m s$^{-1}$]",
        stations=False,
        filename="rh700",
    ),
    "RH850": ProductDefault(
        name="RH",
        mask_var=lambda ds: ds["surface_pressure"].where(
            ds["surface_pressure"] < 85000),
        contourf_var=lambda ds: ds["relhum_850hPa"],
        contourf_kwargs={},
        cmap_label="[%]",
        barb_u=lambda ds: ds["uzonal_850hPa"],
        barb_v=lambda ds: ds["umeridional_850hPa"],
        contour_var=lambda ds: gaussian_filter(
            ds["height_850hPa"].where(ds["surface_pressure"] > 85000) / 10,
            config.mpas_height_smooth_sigma),
        contour_kwargs={"levels": variable_cmaps["GPH850"].contour_levels},
        title="850mb Relative Humidity [%] | Height [dam] | Wind [m s$^{-1}$]",
        stations=False,
        filename="rh850",
    ),
    "RH925": ProductDefault(
        name="RH",
        mask_var=lambda ds: ds["surface_pressure"].where(
            ds["surface_pressure"] < 92500),
        contourf_var=lambda ds: ds["relhum_925hPa"],
        contourf_kwargs={},
        cmap_label="[%]",
        barb_u=lambda ds: ds["uzonal_925hPa"],
        barb_v=lambda ds: ds["umeridional_925hPa"],
        contour_var=lambda ds: gaussian_filter(
            ds["height_925hPa"].where(ds["surface_pressure"] > 92500) / 10,
            config.mpas_height_smooth_sigma),
        contour_kwargs={"levels": variable_cmaps["GPH925"].contour_levels},
        title="925mb Relative Humidity [%] | Height [dam] | Wind [m s$^{-1}$]",
        stations=False,
        filename="rh925",
    ),
    "SNOW3H": ProductDefault(
        name="SNOW",
        contourf_var=lambda ds: ds["snowncv"] * 10,
        contourf_kwargs={"norm": None},
        cmap_label="[cm]",
        title="3hr Accumulated Snow [cm]",
        stations=True,
        station_kwargs={"station_condition": lambda x: x > 0.5},
        filename="snow3h",
    ),
    "SNOWDEPTH": ProductDefault(
        name="SNOW",
        contourf_var=lambda ds: ds["snowh"] * 10,
        contourf_kwargs={"norm": None},
        cmap_label="[cm]",
        title="Snow Depth [cm]",
        stations=True,
        station_kwargs={"station_condition": lambda x: x > 0.5},
        filename="snowdepth",
    ),
    "SNOWT": ProductDefault(
        name="SNOW",
        contourf_var=lambda ds: ds["acsnow"] * 10,
        contourf_kwargs={"norm": None},
        cmap_label="[cm]",
        title="Total Accumulated Snow [cm]",
        stations=True,
        station_kwargs={"station_condition": lambda x: x > 0.5},
        filename="snowt",
    ),
    "SWDOWN": ProductDefault(
        name="SWDOWN",
        contourf_var=lambda ds: ds["swdnb"],
        cmap_label="[W m$^{-2}$]",
        title="Surface Downward Shortwave Radiation Flux [W m$^{-2}$]",
        shapefile_kwargs={"color": "lime"},
        stations=False,
    ),
    "T2M": ProductDefault(
        name="TMP",
        contourf_var=lambda ds: ds["t2m"] - 273.15,
        cmap_label="[°C]",
        barb_u=lambda ds: ds["u10"],
        barb_v=lambda ds: ds["v10"],
        title="Temperature [°C] | 10m Wind [m s$^{-1}$",
        stations=True,
        filename="t2m",
    ),
    "T500": ProductDefault(
        name="TMP",
        mask_var=lambda ds: ds["surface_pressure"].where(
            ds["surface_pressure"] < 50000),
        contourf_var=lambda ds: ds["temperature_500hPa"] - 273.15,
        cmap_label="[°C]",
        barb_u=lambda ds: ds["uzonal_500hPa"],
        barb_v=lambda ds: ds["umeridional_500hPa"],
        contour_var=lambda ds: gaussian_filter(
            ds["height_500hPa"].where(ds["surface_pressure"] > 50000) / 10,
                config.mpas_height_smooth_sigma),
        contour_kwargs={"levels": variable_cmaps["GPH500"].contour_levels},
        title="500mb Temperature [°C] | Height [dam] | Wind [m s$^{-1}$]",
        stations=False,
        filename="t500",
    ),
    "T700": ProductDefault(
        name="TMP",
        mask_var=lambda ds: ds["surface_pressure"].where(
            ds["surface_pressure"] < 70000),
        contourf_var=lambda ds: ds["temperature_700hPa"] - 273.15,
        cmap_label="[°C]",
        barb_u=lambda ds: ds["uzonal_700hPa"],
        barb_v=lambda ds: ds["umeridional_700hPa"],
        contour_var=lambda ds: gaussian_filter(
            ds["height_700hPa"].where(ds["surface_pressure"] > 70000) / 10,
                config.mpas_height_smooth_sigma),
        contour_kwargs={"levels": variable_cmaps["GPH700"].contour_levels},
        title="700mb Temperature [°C] | Height [dam] | Wind [m s$^{-1}$]",
        stations=False,
        filename="t700",
    ),
    "T850": ProductDefault(
        name="TMP",
        mask_var=lambda ds: ds["surface_pressure"].where(
            ds["surface_pressure"] < 85000),
        contourf_var=lambda ds: ds["temperature_850hPa"] - 273.15,
        cmap_label="[°C]",
        barb_u=lambda ds: ds["uzonal_850hPa"],
        barb_v=lambda ds: ds["umeridional_850hPa"],
        contour_var=lambda ds: gaussian_filter(
            ds["height_850hPa"].where(ds["surface_pressure"] > 85000) / 10,
                config.mpas_height_smooth_sigma),
        contour_kwargs={"levels": variable_cmaps["GPH850"].contour_levels},
        title="850mb Temperature [°C] | Height [dam] | Wind [m s$^{-1}$]",
        stations=False,
        filename="t850",
    ),
    "T925": ProductDefault(
        name="TMP",
        mask_var=lambda ds: ds["surface_pressure"].where(
            ds["surface_pressure"] < 92500),
        contourf_var=lambda ds: ds["temperature_925hPa"] - 273.15,
        cmap_label="[°C]",
        barb_u=lambda ds: ds["uzonal_925hPa"],
        barb_v=lambda ds: ds["umeridional_925hPa"],
        contour_var=lambda ds: gaussian_filter(
            ds["height_925hPa"].where(ds["surface_pressure"] > 92500) / 10,
                config.mpas_height_smooth_sigma),
        contour_kwargs={"levels": variable_cmaps["GPH925"].contour_levels},
        title="925mb Temperature [°C] | Height [dam] | Wind [m s$^{-1}$]",
        stations=False,
        filename="t925",
    ),
    "VORT500": ProductDefault(
        name="VORT",
        mask_var=lambda ds: ds["surface_pressure"].where(
            ds["surface_pressure"] < 50000),
        contourf_var=lambda ds: ds["vorticity_500hPa"] * 10**5,
        cmap_label="[x10$^{-5}$ s$^{-1}$]",
        barb_u=lambda ds: ds["uzonal_500hPa"],
        barb_v=lambda ds: ds["umeridional_500hPa"],
        contour_var=lambda ds: gaussian_filter(
            ds["height_500hPa"].where(ds["surface_pressure"] > 50000) / 10,
                config.mpas_height_smooth_sigma),
        contour_kwargs={"levels": variable_cmaps["GPH500"].contour_levels},
        title="500mb Vorticity [x10$^{-5}$ s$^{-1}$] | Height [dam] | Wind [m s$^{-1}$]",
        stations=False,
        filename="vort500",
    ),
    "VORT700": ProductDefault(
        name="VORT",
        mask_var=lambda ds: ds["surface_pressure"].where(
            ds["surface_pressure"] < 70000),
        contourf_var=lambda ds: ds["vorticity_700hPa"] * 10**5,
        cmap_label="[x10$^{-5}$ s$^{-1}$]",
        barb_u=lambda ds: ds["uzonal_700hPa"],
        barb_v=lambda ds: ds["umeridional_700hPa"],
        contour_var=lambda ds: gaussian_filter(
            ds["height_700hPa"].where(ds["surface_pressure"] > 70000) / 10,
                config.mpas_height_smooth_sigma),
        contour_kwargs={"levels": variable_cmaps["GPH700"].contour_levels},
        title="700mb Vorticity [x10$^{-5}$ s$^{-1}$] | Height [dam] | Wind [m s$^{-1}$]",
        stations=False,
        filename="vort700",
    ),
    "VORT850": ProductDefault(
        name="VORT",
        mask_var=lambda ds: ds["surface_pressure"].where(
            ds["surface_pressure"] < 85000),
        contourf_var=lambda ds: ds["vorticity_850hPa"] * 10**5,
        cmap_label="[x10$^{-5}$ s$^{-1}$]",
        barb_u=lambda ds: ds["uzonal_850hPa"],
        barb_v=lambda ds: ds["umeridional_850hPa"],
        contour_var=lambda ds: gaussian_filter(
            ds["height_850hPa"].where(ds["surface_pressure"] > 85000) / 10,
                config.mpas_height_smooth_sigma),
        contour_kwargs={"levels": variable_cmaps["GPH850"].contour_levels},
        title="850mb Vorticity [x10$^{-5}$ s$^{-1}$] | Height [dam] | Wind [m s$^{-1}$]",
        stations=False,
        filename="vort850",
    ),
    "W500": ProductDefault(
        name="W",
        mask_var=lambda ds: ds["surface_pressure"].where(
            ds["surface_pressure"] < 50000),
        contourf_var=lambda ds: -ds["w_500hPa"],
        cmap_label="[m s$^{-1}$]",
        barb_u=lambda ds: ds["u10"],
        barb_v=lambda ds: ds["v10"],
        contour_var=lambda ds: gaussian_filter(
            ds["height_500hPa"].where(ds["surface_pressure"] > 50000) / 10,
            config.mpas_height_smooth_sigma),
        contour_kwargs={"levels": variable_cmaps["GPH500"].contour_levels},
        title="500mb Vertical Velocity [m s$^{-1}$] | Height [dam] | Wind [m s$^{-1}$]",
        stations=False,
        filename="w500",
    ),
    "W700": ProductDefault(
        name="W",
        mask_var=lambda ds: ds["surface_pressure"].where(
            ds["surface_pressure"] < 70000),
        contourf_var=lambda ds: -ds["w_700hPa"],
        cmap_label="[m s$^{-1}$]",
        barb_u=lambda ds: ds["u10"],
        barb_v=lambda ds: ds["v10"],
        contour_var=lambda ds: gaussian_filter(
            ds["height_700hPa"].where(ds["surface_pressure"] > 70000) / 10,
            config.mpas_height_smooth_sigma),
        contour_kwargs={"levels": variable_cmaps["GPH700"].contour_levels},
        title="700mb Vertical Velocity [m s$^{-1}$] | Height [dam] | Wind [m s$^{-1}$]",
        stations=False,
        filename="w700",
    ),
    "W850": ProductDefault(
        name="W",
        mask_var=lambda ds: ds["surface_pressure"].where(
            ds["surface_pressure"] < 85000),
        contourf_var=lambda ds: -ds["w_850hPa"],
        cmap_label="[m s$^{-1}$]",
        barb_u=lambda ds: ds["u10"],
        barb_v=lambda ds: ds["v10"],
        contour_var=lambda ds: gaussian_filter(
            ds["height_850hPa"].where(ds["surface_pressure"] > 85000) / 10,
            config.mpas_height_smooth_sigma),
        contour_kwargs={"levels": variable_cmaps["GPH850"].contour_levels},
        title="850mb Vertical Velocity [m s$^{-1}$] | Height [dam] | Wind [m s$^{-1}$]",
        stations=False,
        filename="w850",
    ),
    "WIND10M": ProductDefault(
        name="WIND10M",
        contourf_var=lambda ds: np.sqrt(ds["u10"]**2 + ds["v10"]**2),
        cmap_label="[m s$^{-1}$]",
        barb_u=lambda ds: ds["u10"],
        barb_v=lambda ds: ds["v10"],
        contour_var=lambda ds: gaussian_filter(
            ds["mslp"] / 100, config.mpas_slp_smooth_sigma),
        contour_kwargs={"levels": variable_cmaps["SLP"].contour_levels},
        title="10m Wind [m s$^{-1}$] | MSLP [mb]",
        stations=True,
        hilo_var=lambda ds: ds["mslp"] / 100,
    ),
    "WIND250": ProductDefault(
        name="WIND250",
        mask_var=lambda ds: ds["surface_pressure"].where(
            ds["surface_pressure"] < 25000),
        contourf_var=lambda ds: np.sqrt(
            ds["uzonal_250hPa"]**2 + ds["umeridional_250hPa"]**2),
        cmap_label="[m s$^{-1}$]",
        barb_u=lambda ds: ds["uzonal_250hPa"],
        barb_v=lambda ds: ds["umeridional_250hPa"],
        contour_var=lambda ds: gaussian_filter(
            ds["height_250hPa"].where(ds["surface_pressure"] > 25000) / 10,
            config.mpas_height_smooth_sigma),
        contour_kwargs={"levels": variable_cmaps["GPH250"].contour_levels},
        title="250mb Wind [m s$^{-1}$] | Height [dam]",
        stations=False,
    ),
    "WIND500": ProductDefault(
        name="WIND500",
        mask_var=lambda ds: ds["surface_pressure"].where(
            ds["surface_pressure"] < 50000),
        contourf_var=lambda ds: np.sqrt(
            ds["uzonal_500hPa"]**2 + ds["umeridional_500hPa"]**2),
        cmap_label="[m s$^{-1}$]",
        barb_u=lambda ds: ds["uzonal_500hPa"],
        barb_v=lambda ds: ds["umeridional_500hPa"],
        contour_var=lambda ds: gaussian_filter(
            ds["height_500hPa"].where(ds["surface_pressure"] > 50000) / 10,
            config.mpas_height_smooth_sigma),
        contour_kwargs={"levels": variable_cmaps["GPH500"].contour_levels},
        title="500mb Wind [m s$^{-1}$] | Height [dam]",
        stations=False,
    ),
    "WIND700": ProductDefault(
        name="WIND700",
        mask_var=lambda ds: ds["surface_pressure"].where(
            ds["surface_pressure"] < 70000),
        contourf_var=lambda ds: np.sqrt(
            ds["uzonal_700hPa"]**2 + ds["umeridional_700hPa"]**2),
        cmap_label="[m s$^{-1}$]",
        barb_u=lambda ds: ds["uzonal_700hPa"],
        barb_v=lambda ds: ds["umeridional_700hPa"],
        contour_var=lambda ds: gaussian_filter(
            ds["height_700hPa"].where(ds["surface_pressure"] > 70000) / 10,
            config.mpas_height_smooth_sigma),
        contour_kwargs={"levels": variable_cmaps["GPH700"].contour_levels},
        title="700mb Wind [m s$^{-1}$] | Height [dam]",
        stations=False,
    ),
    "WIND850": ProductDefault(
        name="WIND850",
        mask_var=lambda ds: ds["surface_pressure"].where(
            ds["surface_pressure"] < 85000),
        contourf_var=lambda ds: np.sqrt(
            ds["uzonal_850hPa"]**2 + ds["umeridional_850hPa"]**2),
        cmap_label="[m s$^{-1}$]",
        barb_u=lambda ds: ds["uzonal_850hPa"],
        barb_v=lambda ds: ds["umeridional_850hPa"],
        contour_var=lambda ds: gaussian_filter(
            ds["height_850hPa"].where(ds["surface_pressure"] > 85000) / 10,
            config.mpas_height_smooth_sigma),
        contour_kwargs={"levels": variable_cmaps["GPH850"].contour_levels},
        title="850mb Wind [m s$^{-1}$] | Height [dam]",
        stations=False,
    ),
    "WIND925": ProductDefault(
        name="WIND925",
        mask_var=lambda ds: ds["surface_pressure"].where(
            ds["surface_pressure"] < 92500),
        contourf_var=lambda ds: np.sqrt(
            ds["uzonal_925hPa"]**2 + ds["umeridional_925hPa"]**2),
        cmap_label="[m s$^{-1}$]",
        barb_u=lambda ds: ds["uzonal_925hPa"],
        barb_v=lambda ds: ds["umeridional_925hPa"],
        contour_var=lambda ds: gaussian_filter(
            ds["height_925hPa"].where(ds["surface_pressure"] > 92500) / 10,
            config.mpas_height_smooth_sigma),
        contour_kwargs={"levels": variable_cmaps["GPH925"].contour_levels},
        title="925mb Wind [m s$^{-1}$] | Height [dam]",
        stations=False,
    ),
}


wrf_products = lambda domain: {
    "CAPEMU" :ProductDefault(
        name="CAPE",
        contourf_var=lambda ds: ds.ds["capemu"],
        cmap_label="[J kg$^{-1}$]",
        barb_u=lambda ds: ds.ds["umet10"].where(ds.ds["capemu"] > 50),
        barb_v=lambda ds: ds.ds["vmet10"].where(ds.ds["capemu"] > 50),
        barb_kwargs={"pivot": "tip"},
        barb_u_2=lambda ds: wrf_calc.interp_ds(ds, "umet", "prs", 85000).where(
            ds.ds["capemu"] > 50),
        barb_v_2=lambda ds: wrf_calc.interp_ds(ds, "vmet", "prs", 85000).where(
            ds.ds["capemu"] > 50),
        barb_kwargs_2={"barbcolor": "red", "pivot": "tip"},
        barb_u_3=lambda ds: wrf_calc.interp_ds(ds, "umet", "prs", 50000).where(
            ds.ds["capemu"] > 50),
        barb_v_3=lambda ds: wrf_calc.interp_ds(ds, "umet", "prs", 50000).where(
            ds.ds["capemu"] > 50),
        barb_kwargs_3={"barbcolor": "green", "pivot": "tip"},
        title="Most Unstable CAPE [J kg$^{-1}$] | 10m/850mb/500mb Wind (black/red/green)",
        filename="capemu",
        stations=False,
        scales=99
    ),
    "CIN": ProductDefault(
        name="CIN",
        contourf_var=lambda ds: -ds.ds["cin"],
        title="CIN [J kg$^{-1}$]",
        stations=False,
        scales=99
    ),
    "CLOUDCOVERHIGH": ProductDefault(
        name="CLOUDCOVER",
        contourf_var=lambda ds: ds.ds["cloud_cover_high"] * 100,
        title="High Cloud Cover [%]",
        shapefile_kwargs={"color": "lime"},
        filename="cloudcoverhigh",
        stations=False,
        scales=99
    ),
    "CLOUDCOVERLOW": ProductDefault(
        name="CLOUDCOVER",
        contourf_var=lambda ds: ds.ds["cloud_cover_low"] * 100,
        title="Low Cloud Cover [%]",
        shapefile_kwargs={"color": "lime"},
        filename="cloudcoverlow",
        stations=False,
        scales=99
    ),
    "CLOUDCOVERMID": ProductDefault(
        name="CLOUDCOVER",
        contourf_var=lambda ds: ds.ds["cloud_cover_mid"] * 100,
        title="Mid Cloud Cover [%]",
        shapefile_kwargs={"color": "lime"},
        filename="cloudcovermid",
        stations=False,
        scales=99
    ),
    "CTT": ProductDefault(
        name="CTT",
        contourf_var=lambda ds: ds.ds["ctt"],
        contourf_kwargs={"norm": None},
        title="Cloud Top Temperature [°C]",
        stations=False,
        scales=99
    ),
    "DBZ": ProductDefault(
        name="DBZ",
        contourf_var=lambda ds: ds.ds["dbz"].max("bottom_top"),
        cmap_label="[dBZ]",
        barb_u=lambda ds: ds.ds["umet10"],
        barb_v=lambda ds: ds.ds["vmet10"],
        contour_var=lambda ds: gaussian_filter(
            ds.ds["slp"], config.wrf_slp_smooth_sigma[domain]),
        contour_kwargs={"levels": variable_cmaps["SLP"].contour_levels},
        title="Composite Reflectivity [dBZ] | 10m Wind [m s$^{-1}$] | MSLP [mb]",
        filename="dbzcomp",
        stations=False,
        hilo_var=lambda ds: ds.ds["slp"],
        hilo_kwargs={"gaussian_sigma": config.wrf_slp_smooth_sigma[domain]},
        scales=99
    ),
    "DBZ1KM": ProductDefault(
        name="DBZ",
        contourf_var=lambda ds: wrf_calc.interplevel(
            ds.ds["dbz"], ds.ds["hgt"], 1000),
        cmap_label="[dBZ]",
        barb_u=lambda ds: ds.ds["umet10"],
        barb_v=lambda ds: ds.ds["vmet10"],
        contour_var=lambda ds: gaussian_filter(
            ds.ds["slp"], config.wrf_slp_smooth_sigma[domain]),
        contour_kwargs={"levels": variable_cmaps["SLP"].contour_levels},
        title="1km Reflectivity [dBZ] | 10m Wind [m s$^{-1}$] | MSLP [mb]",
        filename="dbz1km",
        stations=False,
        hilo_var=lambda ds: ds.ds["slp"],
        hilo_kwargs={"gaussian_sigma": config.wrf_slp_smooth_sigma[domain]},
        scales=99
    ),
    "DPT2M": ProductDefault(
        name="TMP",
        contourf_var=lambda ds: ds.ds["dpt2m"],
        cmap_label="[°C]",
        barb_u=lambda ds: ds.ds["umet10"],
        barb_v=lambda ds: ds.ds["vmet10"],
        contour_var=lambda ds: gaussian_filter(
            ds.ds["slp"], config.wrf_slp_smooth_sigma[domain]),
        contour_kwargs={"levels": variable_cmaps["SLP"].contour_levels},
        title="2m Dewpoint [°C] | 10m Wind [m s$^{-1}$] | MSLP [mb]",
        filename="dpt2m",
        hilo_var=lambda ds: ds.ds["slp"],
        hilo_kwargs={"gaussian_sigma": config.wrf_slp_smooth_sigma[domain]},
        scales=99
    ),
    "IVT": ProductDefault(
        name="IVT",
        contourf_var=lambda ds: ds.ds["ivt"],
        cmap_label="[kg m$^{-1}$ s$^{-1}$]",
        quiver_u=lambda ds: ds.ds["ivt_u"].where(ds.ds["ivt"] > 250),
        quiver_v=lambda ds: ds.ds["ivt_v"].where(ds.ds["ivt"] > 250),
        quiver_kwargs={"scale": 50000},
        contour_var=lambda ds: gaussian_filter(
            ds.ds["slp"], config.wrf_slp_smooth_sigma[domain]),
        contour_kwargs={"levels": variable_cmaps["SLP"].contour_levels},
        title="Integrated Vapor Transport [kg m$^{-1}$ s$^{-1}$] | MSLP [mb]",
        stations=False,
        hilo_var=lambda ds: ds.ds["slp"],
        hilo_kwargs={"gaussian_sigma": config.wrf_slp_smooth_sigma[domain]},
        scales=2
    ),
    "LR700500": ProductDefault(
        name="LR",
        contourf_var=lambda ds: -(wrf_calc.interp_ds(
            ds, "tmp", "prs", 70000) - wrf_calc.interp_ds(
            ds, "tmp", "prs", 50000)) / \
        (wrf_calc.interp_ds(ds, "hgt", "prs", 70000) - \
            wrf_calc.interp_ds(ds, "hgt", "prs", 50000)) * 1000,
        title="700-500mb Lapse Rate [°C km$^{-1}$]",
        filename="lr700500",
        stations=False,
        scales=2
    ),
    "OLR": ProductDefault(
        name="OLR",
        contourf_var=lambda ds: ds.ds["OLR"],
        title="Top of Atmosphere Outgoing Longwave Radiation [W m$^{-2}$]",
        shapefile_kwargs={"color": "lime"},
        stations=False,
        scales=99
    ),
    "OMEGA500": ProductDefault(
        name="OMEGA",
        contourf_var=lambda ds: wrf_calc.interp_ds(
            ds, "omega", "prs", 50000) * 10,
        cmap_label="[dPa s$^{-1}$]",
        barb_u=lambda ds: wrf_calc.interp_ds(ds, "umet", "prs", 50000),
        barb_v=lambda ds: wrf_calc.interp_ds(ds, "vmet", "prs", 50000),
        contour_var=lambda ds: gaussian_filter(
            wrf_calc.interp_ds(ds, "hgt", "prs", 50000),
            config.wrf_height_smooth_sigma[domain]) / 10,
        mask_var=lambda ds: ds.ds["PSFC"].where(ds.ds["PSFC"] < 50000),
        title="500mb Omega [dPa s$^{-1}$] | Wind [m s$^{-1}$] | Height [dam]",
        filename="omega500",
        stations=False,
        scales=2
    ),
    "OMEGA700": ProductDefault(
        name="OMEGA",
        contourf_var=lambda ds: wrf_calc.interp_ds(
            ds, "omega", "prs", 70000) * 10,
        cmap_label="[dPa s$^{-1}$]",
        barb_u=lambda ds: wrf_calc.interp_ds(ds, "umet", "prs", 70000),
        barb_v=lambda ds: wrf_calc.interp_ds(ds, "vmet", "prs", 70000),
        contour_var=lambda ds: gaussian_filter(
            wrf_calc.interp_ds(ds, "hgt", "prs", 70000),
            config.wrf_height_smooth_sigma[domain]) / 10,
        mask_var=lambda ds: ds.ds["PSFC"].where(ds.ds["PSFC"] < 70000),
        title="700mb Omega [dPa s$^{-1}$] | Wind [m s$^{-1}$] | Height [dam]",
        filename="omega700",
        stations=False,
        scales=2
    ),
    "OMEGA850": ProductDefault(
        name="OMEGA",
        contourf_var=lambda ds: wrf_calc.interp_ds(
            ds, "omega", "prs", 85000) * 10,
        cmap_label="[dPa s$^{-1}$]",
        barb_u=lambda ds: wrf_calc.interp_ds(ds, "umet", "prs", 85000),
        barb_v=lambda ds: wrf_calc.interp_ds(ds, "vmet", "prs", 85000),
        contour_var=lambda ds: gaussian_filter(
            wrf_calc.interp_ds(ds, "hgt", "prs", 85000),
            config.wrf_height_smooth_sigma[domain]) / 10,
        mask_var=lambda ds: ds.ds["PSFC"].where(ds.ds["PSFC"] < 85000),
        title="850mb Omega [dPa s$^{-1}$] | Wind [m s$^{-1}$] | Height [dam]",
        filename="omega850",
        stations=False,
        scales=2
    ),
    "PWAT": ProductDefault(
        name="PWAT",
        contourf_var=lambda ds: ds.ds["pwat"],
        cmap_label="[mm]",
        barb_u=lambda ds: wrf_calc.interp_ds(ds, "umet", "prs", 85000),
        barb_v=lambda ds: wrf_calc.interp_ds(ds, "vmet", "prs", 85000),
        contour_var=lambda ds: gaussian_filter(
            wrf_calc.interp_ds(ds, "hgt", "prs", 85000),
            config.wrf_height_smooth_sigma[domain]) / 10,
        title="Precipitable Water [mm] | 850mb Wind [m s$^{-1}$] | Height [dam]",
        stations=False,
        scales=2
    ),
    "QPF1H": ProductDefault(
        name="QPF",
        contourf_var=lambda ds: ds.ds["PREC_ACC_C"] + ds.ds["PREC_ACC_NC"],
        contourf_kwargs={"norm": None},
        title="1hr Precipitation [mm]",
        filename="qpf1h",
        station_kwargs={"station_condition": lambda x: x > 0.5},
        scales=99
    ),
    "QPFT": ProductDefault(
        name="QPF",
        contourf_var=lambda ds: ds.ds["RAINNC"] + ds.ds["RAINC"],
        contourf_kwargs={"norm": None},
        title="Total Precipitation [mm]",
        filename="qpft",
        station_kwargs={"station_condition": lambda x: x > 0.5},
        scales=99
    ),
    "RH2M": ProductDefault(
        name="RH",
        contourf_var=lambda ds: ds.ds["rh2m"],
        cmap_label="[%]",
        barb_u=lambda ds: ds.ds["umet10"],
        barb_v=lambda ds: ds.ds["vmet10"],
        contour_var=lambda ds: gaussian_filter(
            ds.ds["slp"], config.wrf_slp_smooth_sigma[domain]),
        contour_kwargs={"levels": variable_cmaps["SLP"].contour_levels},
        title="2m Relative Humidity [%] | 10m Wind [m s$^{-1}}$] | MSLP [mb]",
        filename="rh2m",
        hilo_var=lambda ds: ds.ds["slp"],
        hilo_kwargs={"gaussian_sigma": config.wrf_slp_smooth_sigma[domain]},
        scales=99
    ),
    "RH250": ProductDefault(
        name="RH",
        contourf_var=lambda ds: wrf_calc.interp_ds(
            ds, "rh", "prs", 25000) * 100,
        cmap_label="[%]",
        barb_u=lambda ds: wrf_calc.interp_ds(ds, "umet", "prs", 25000),
        barb_v=lambda ds: wrf_calc.interp_ds(ds, "vmet", "prs", 25000),
        contour_var=lambda ds: gaussian_filter(
            wrf_calc.interp_ds(ds, "hgt", "prs", 25000),
            config.wrf_height_smooth_sigma[domain]) / 10,
        mask_var=lambda ds: ds.ds["PSFC"].where(ds.ds["PSFC"] < 25000),
        title="250mb Relative Humidity [%] | Wind [m s$^{-1}$] | Height [dam]",
        filename="rh250",
        stations=False,
        scales=2
    ),
    "RH500": ProductDefault(
        name="RH",
        contourf_var=lambda ds: wrf_calc.interp_ds(
            ds, "rh", "prs", 50000) * 100,
        cmap_label="[%]",
        barb_u=lambda ds: wrf_calc.interp_ds(ds, "umet", "prs", 50000),
        barb_v=lambda ds: wrf_calc.interp_ds(ds, "vmet", "prs", 50000),
        contour_var=lambda ds: gaussian_filter(
            wrf_calc.interp_ds(ds, "hgt", "prs", 50000),
            config.wrf_height_smooth_sigma[domain]) / 10,
        mask_var=lambda ds: ds.ds["PSFC"].where(ds.ds["PSFC"] < 50000),
        title="500mb Relative Humidity [%] | Wind [m s$^{-1}$] | Height [dam]",
        filename="rh500",
        stations=False,
        scales=2
    ),
    "RH700": ProductDefault(
        name="RH",
        contourf_var=lambda ds: wrf_calc.interp_ds(
            ds, "rh", "prs", 70000) * 100,
        cmap_label="[%]",
        barb_u=lambda ds: wrf_calc.interp_ds(ds, "umet", "prs", 70000),
        barb_v=lambda ds: wrf_calc.interp_ds(ds, "vmet", "prs", 70000),
        contour_var=lambda ds: gaussian_filter(
            wrf_calc.interp_ds(ds, "hgt", "prs", 70000),
            config.wrf_height_smooth_sigma[domain]) / 10,
        mask_var=lambda ds: ds.ds["PSFC"].where(ds.ds["PSFC"] < 70000),
        title="700mb Relative Humidity [%] | Wind [m s$^{-1}$] | Height [dam]",
        filename="rh700",
        stations=False,
        scales=2
    ),
    "RH850": ProductDefault(
        name="RH",
        contourf_var=lambda ds: wrf_calc.interp_ds(
            ds, "rh", "prs", 85000) * 100,
        cmap_label="[%]",
        barb_u=lambda ds: wrf_calc.interp_ds(ds, "umet", "prs", 85000),
        barb_v=lambda ds: wrf_calc.interp_ds(ds, "vmet", "prs", 85000),
        contour_var=lambda ds: gaussian_filter(
            wrf_calc.interp_ds(ds, "hgt", "prs", 85000),
            config.wrf_height_smooth_sigma[domain]) / 10,
        mask_var=lambda ds: ds.ds["PSFC"].where(ds.ds["PSFC"] < 85000),
        title="850mb Relative Humidity [%] | Wind [m s$^{-1}$] | Height [dam]",
        filename="rh850",
        stations=False,
        scales=2
    ),
    "RH925": ProductDefault(
        name="RH",
        contourf_var=lambda ds: wrf_calc.interp_ds(
            ds, "rh", "prs", 92500) * 100,
        cmap_label="[%]",
        barb_u=lambda ds: wrf_calc.interp_ds(ds, "umet", "prs", 92500),
        barb_v=lambda ds: wrf_calc.interp_ds(ds, "vmet", "prs", 92500),
        contour_var=lambda ds: gaussian_filter(
            wrf_calc.interp_ds(ds, "hgt", "prs", 92500),
            config.wrf_height_smooth_sigma[domain]) / 10,
        mask_var=lambda ds: ds.ds["PSFC"].where(ds.ds["PSFC"] < 92500),
        title="925mb Relative Humidity [%] | Wind [m s$^{-1}$] | Height [dam]",
        filename="rh925",
        stations=False,
        scales=2
    ),
    "SNOW1H": ProductDefault(
        name="SNOW",
        contourf_var=lambda ds: ds.ds["SNOW_ACC_NC"] * 10,
        contourf_kwargs={"norm": None},
        title="1hr Snow 10:1 [cm]",
        filename="snow1h",
        station_kwargs={"station_condition": lambda x: x > 0.5},
        scales=99
    ),
    "SNOWT": ProductDefault(
        name="SNOW",
        contourf_var=lambda ds: ds.ds["ACSNOW"] * 10,
        contourf_kwargs={"norm": None},
        title="Total Snow 10:1 [cm]",
        filename="snowt",
        station_kwargs={"station_condition": lambda x: x > 0.5},
        scales=99
    ),
    "SNOWDEPTH": ProductDefault(
        name="SNOW",
        contourf_var=lambda ds: ds.ds["SNOWH"] * 10,
        contourf_kwargs={"norm": None},
        title="Snow Depth [cm]",
        station_kwargs={"station_condition": lambda x: x > 0.5},
        scales=99
    ),
    "SWDOWN": ProductDefault(
        name="SWDOWN",
        contourf_var=lambda ds: ds.ds["SWDOWN"],
        title="Surface Downward Shortwave Radiation [W m$^{-2}$]",
        stations=False,
        scales=99
    ),
    "T2M": ProductDefault(
        name="TMP",
        contourf_var=lambda ds: ds.ds["T2"] - 273.15,
        cmap_label="[°C]",
        barb_u=lambda ds: ds.ds["umet10"],
        barb_v=lambda ds: ds.ds["vmet10"],
        contour_var=lambda ds: gaussian_filter(
            ds.ds["slp"], config.wrf_slp_smooth_sigma[domain]),
        contour_kwargs={"levels": variable_cmaps["SLP"].contour_levels},
        title="2m Temperature [°C] | 10m Wind [m s$^{-1}$] | MSLP [mb]",
        filename="t2m",
        hilo_var=lambda ds: ds.ds["slp"],
        hilo_kwargs={"gaussian_sigma": config.wrf_slp_smooth_sigma[domain]},
        scales=99
    ),
    "TA2M": ProductDefault(
        name="TMP",
        contourf_var=lambda ds: ds.ds["ta2"] - 273.15,
        cmap_label="[°C]",
        barb_u=lambda ds: ds.ds["umet10"],
        barb_v=lambda ds: ds.ds["vmet10"],
        contour_var=lambda ds: gaussian_filter(
            ds.ds["slp"], config.wrf_slp_smooth_sigma[domain]),
        contour_kwargs={"levels": variable_cmaps["SLP"].contour_levels},
        title="2m Apparent Temperature [°C] | 10m Wind [m s$^{-1}$] | MSLP [mb]",
        filename="ta2m",
        hilo_var=lambda ds: ds.ds["slp"],
        hilo_kwargs={"gaussian_sigma": config.wrf_slp_smooth_sigma[domain]},
        scales=99
    ),
    "T500": ProductDefault(
        name="TMP",
        contourf_var=lambda ds: wrf_calc.interp_ds(
            ds, "tmp", "prs", 50000) - 273.15,
        cmap_label="[°C]",
        barb_u=lambda ds: wrf_calc.interp_ds(ds, "umet", "prs", 50000),
        barb_v=lambda ds: wrf_calc.interp_ds(ds, "vmet", "prs", 50000),
        contour_var=lambda ds: gaussian_filter(
            wrf_calc.interp_ds(ds, "hgt", "prs", 50000),
            config.wrf_height_smooth_sigma[domain]) / 10,
        mask_var=lambda ds: ds.ds["PSFC"].where(ds.ds["PSFC"] < 50000),
        title="500mb Temperature [°C] | Wind [m s$^{-1}$] | Height [dam]",
        filename="t500",
        stations=False,
        scales=2
    ),
    "T700": ProductDefault(
        name="TMP",
        contourf_var=lambda ds: wrf_calc.interp_ds(
            ds, "tmp", "prs", 70000) - 273.15,
        cmap_label="[°C]",
        barb_u=lambda ds: wrf_calc.interp_ds(ds, "umet", "prs", 70000),
        barb_v=lambda ds: wrf_calc.interp_ds(ds, "vmet", "prs", 70000),
        contour_var=lambda ds: gaussian_filter(
            wrf_calc.interp_ds(ds, "hgt", "prs", 70000),
            config.wrf_height_smooth_sigma[domain]) / 10,
        mask_var=lambda ds: ds.ds["PSFC"].where(ds.ds["PSFC"] < 70000),
        title="700mb Temperature [°C] | Wind [m s$^{-1}$] | Height [dam]",
        filename="t700",
        stations=False,
        scales=2
    ),
    "T850": ProductDefault(
        name="TMP",
        contourf_var=lambda ds: wrf_calc.interp_ds(
            ds, "tmp", "prs", 85000) - 273.15,
        cmap_label="[°C]",
        barb_u=lambda ds: wrf_calc.interp_ds(ds, "umet", "prs", 85000),
        barb_v=lambda ds: wrf_calc.interp_ds(ds, "vmet", "prs", 85000),
        contour_var=lambda ds: gaussian_filter(
            wrf_calc.interp_ds(ds, "hgt", "prs", 85000),
            config.wrf_height_smooth_sigma[domain]) / 10,
        mask_var=lambda ds: ds.ds["PSFC"].where(ds.ds["PSFC"] < 85000),
        title="850mb Temperature [°C] | Wind [m s$^{-1}$] | Height [dam]",
        filename="t850",
        stations=False,
        scales=2
    ),
    "T925": ProductDefault(
        name="TMP",
        contourf_var=lambda ds: wrf_calc.interp_ds(
            ds, "tmp", "prs", 92500) - 273.15,
        cmap_label="[°C]",
        barb_u=lambda ds: wrf_calc.interp_ds(ds, "umet", "prs", 92500),
        barb_v=lambda ds: wrf_calc.interp_ds(ds, "vmet", "prs", 92500),
        contour_var=lambda ds: gaussian_filter(
            wrf_calc.interp_ds(ds, "hgt", "prs", 92500),
            config.wrf_height_smooth_sigma[domain]) / 10,
        mask_var=lambda ds: ds.ds["PSFC"].where(ds.ds["PSFC"] < 92500),
        title="925mb Temperature [°C] | Wind [m s$^{-1}$] | Height [dam]",
        filename="t925",
        stations=False,
        scales=2
    ),
    "VORT500": ProductDefault(
        name="VORT",
        contourf_var=lambda ds: wrf_calc.interp_ds(
            ds, "vort", "prs", 50000) * 1e5,
        cmap_label="[x10$^{-5}$ s$^{-1}$]",
        barb_u=lambda ds: wrf_calc.interp_ds(ds, "umet", "prs", 50000),
        barb_v=lambda ds: wrf_calc.interp_ds(ds, "vmet", "prs", 50000),
        contour_var=lambda ds: gaussian_filter(
            wrf_calc.interp_ds(ds, "hgt", "prs", 50000),
            config.wrf_height_smooth_sigma[domain]) / 10,
        mask_var=lambda ds: ds.ds["PSFC"].where(ds.ds["PSFC"] < 50000),
        title="500mb Vorticity [x10$^{-5}$ s$^{-1}$] | Wind [m s$^{-1}$] | Height [dam]",
        filename="vort500",
        stations=False,
        scales=2
    ),
    "VORT700": ProductDefault(
        name="VORT",
        contourf_var=lambda ds: wrf_calc.interp_ds(
            ds, "vort", "prs", 70000) * 1e5,
        cmap_label="[x10$^{-5}$ s$^{-1}$]",
        barb_u=lambda ds: wrf_calc.interp_ds(ds, "umet", "prs", 70000),
        barb_v=lambda ds: wrf_calc.interp_ds(ds, "vmet", "prs", 70000),
        contour_var=lambda ds: gaussian_filter(
            wrf_calc.interp_ds(ds, "hgt", "prs", 70000),
            config.wrf_height_smooth_sigma[domain]) / 10,
        mask_var=lambda ds: ds.ds["PSFC"].where(ds.ds["PSFC"] < 70000),
        title="700mb Vorticity [x10$^{-5}$ s$^{-1}$] | Wind [m s$^{-1}$] | Height [dam]",
        filename="vort700",
        stations=False,
        scales=2
    ),
    "VORT850": ProductDefault(
        name="VORT",
        contourf_var=lambda ds: wrf_calc.interp_ds(
            ds, "vort", "prs", 85000) * 1e5,
        cmap_label="[x10$^{-5}$ s$^{-1}$]",
        barb_u=lambda ds: wrf_calc.interp_ds(ds, "umet", "prs", 85000),
        barb_v=lambda ds: wrf_calc.interp_ds(ds, "vmet", "prs", 85000),
        contour_var=lambda ds: gaussian_filter(
            wrf_calc.interp_ds(ds, "hgt", "prs", 85000),
            config.wrf_height_smooth_sigma[domain]) / 10,
        mask_var=lambda ds: ds.ds["PSFC"].where(ds.ds["PSFC"] < 85000),
        title="850mb Vorticity [x10$^{-5}$ s$^{-1}$] | Wind [m s$^{-1}$] | Height [dam]",
        filename="vort850",
        stations=False,
        scales=2
    ),
    "WIND10M": ProductDefault(
        name="WIND10M",
        contourf_var=lambda ds: wrf_calc.wind_speed(
            ds.ds["umet10"], ds.ds["vmet10"]),
        cmap_label="[m s$^{-1}$]",
        barb_u=lambda ds: ds.ds["umet10"],
        barb_v=lambda ds: ds.ds["vmet10"],
        contour_var=lambda ds: gaussian_filter(
            ds.ds["slp"], config.wrf_slp_smooth_sigma[domain]),
        contour_kwargs={"levels": variable_cmaps["SLP"].contour_levels},
        title="10m Wind [m s$^{-1}$] | MSLP [mb]",
        filename="wind10m",
        hilo_var=lambda ds: ds.ds["slp"],
        hilo_kwargs={"gaussian_sigma": config.wrf_slp_smooth_sigma[domain]},
        scales=99
    ),
    "WIND250": ProductDefault(
        name="WIND250",
        contourf_var=lambda ds: wrf_calc.wind_speed(
            wrf_calc.interp_ds(ds, "umet", "prs", 25000),
            wrf_calc.interp_ds(ds, "vmet", "prs", 25000)),
        cmap_label="[m s$^{-1}$]",
        barb_u=lambda ds: wrf_calc.interp_ds(ds, "umet", "prs", 25000),
        barb_v=lambda ds: wrf_calc.interp_ds(ds, "vmet", "prs", 25000),
        contour_var=lambda ds: gaussian_filter(
            wrf_calc.interp_ds(ds, "hgt", "prs", 25000),
            config.wrf_height_smooth_sigma[domain]) / 10,
        mask_var=lambda ds: ds.ds["PSFC"].where(ds.ds["PSFC"] < 25000),
        title="250mb Wind [m s$^{-1}$] | Height [dam]",
        filename="wind250",
        stations=False,
        scales=2
    ),
    "WIND500": ProductDefault(
        name="WIND500",
        contourf_var=lambda ds: wrf_calc.wind_speed(
            wrf_calc.interp_ds(ds, "umet", "prs", 50000),
            wrf_calc.interp_ds(ds, "vmet", "prs", 50000)),
        cmap_label="[m s$^{-1}$]",
        barb_u=lambda ds: wrf_calc.interp_ds(ds, "umet", "prs", 50000),
        barb_v=lambda ds: wrf_calc.interp_ds(ds, "vmet", "prs", 50000),
        contour_var=lambda ds: gaussian_filter(
            wrf_calc.interp_ds(ds, "hgt", "prs", 50000),
            config.wrf_height_smooth_sigma[domain]) / 10,
        mask_var=lambda ds: ds.ds["PSFC"].where(ds.ds["PSFC"] < 50000),
        title="500mb Wind [m s$^{-1}$] | Height [dam]",
        filename="wind500",
        stations=False,
        scales=2
    ),
    "WIND700": ProductDefault(
        name="WIND700",
        contourf_var=lambda ds: wrf_calc.wind_speed(
            wrf_calc.interp_ds(ds, "umet", "prs", 70000),
            wrf_calc.interp_ds(ds, "vmet", "prs", 70000)),
        cmap_label="[m s$^{-1}$]",
        barb_u=lambda ds: wrf_calc.interp_ds(ds, "umet", "prs", 70000),
        barb_v=lambda ds: wrf_calc.interp_ds(ds, "vmet", "prs", 70000),
        contour_var=lambda ds: gaussian_filter(
            wrf_calc.interp_ds(ds, "hgt", "prs", 70000),
            config.wrf_height_smooth_sigma[domain]) / 10,
        mask_var=lambda ds: ds.ds["PSFC"].where(ds.ds["PSFC"] < 70000),
        title="700mb Wind [m s$^{-1}$] | Height [dam]",
        filename="wind700",
        stations=False,
        scales=2
    ),
    "WIND850": ProductDefault(
        name="WIND850",
        contourf_var=lambda ds: wrf_calc.wind_speed(
            wrf_calc.interp_ds(ds, "umet", "prs", 85000),
            wrf_calc.interp_ds(ds, "vmet", "prs", 85000)),
        cmap_label="[m s$^{-1}$]",
        barb_u=lambda ds: wrf_calc.interp_ds(ds, "umet", "prs", 85000),
        barb_v=lambda ds: wrf_calc.interp_ds(ds, "vmet", "prs", 85000),
        contour_var=lambda ds: gaussian_filter(
            wrf_calc.interp_ds(ds, "hgt", "prs", 85000),
            config.wrf_height_smooth_sigma[domain]) / 10,
        mask_var=lambda ds: ds.ds["PSFC"].where(ds.ds["PSFC"] < 85000),
        title="850mb Wind [m s$^{-1}$] | Height [dam]",
        filename="wind850",
        stations=False,
        scales=2
    ),
    "WIND925": ProductDefault(
        name="WIND925",
        contourf_var=lambda ds: wrf_calc.wind_speed(
            wrf_calc.interp_ds(ds, "umet", "prs", 92500),
            wrf_calc.interp_ds(ds, "vmet", "prs", 92500)),
        cmap_label="[m s$^{-1}$]",
        barb_u=lambda ds: wrf_calc.interp_ds(ds, "umet", "prs", 92500),
        barb_v=lambda ds: wrf_calc.interp_ds(ds, "vmet", "prs", 92500),
        contour_var=lambda ds: gaussian_filter(
            wrf_calc.interp_ds(ds, "hgt", "prs", 92500),
            config.wrf_height_smooth_sigma[domain]) / 10,
        mask_var=lambda ds: ds.ds["PSFC"].where(ds.ds["PSFC"] < 92500),
        title="925mb Wind [m s$^{-1}$] | Height [dam]",
        filename="wind925",
        stations=False,
        scales=2
    ),
}


class TimeSeriesDefaults:
    def __init__(self, name, var):
        self.name = name
        self.var = var


wrf_time_series_products = {
    "DPT2M": TimeSeriesDefaults(
        name="DPT2M", var=lambda ds: wrf_calc.dpt2m(ds)),
    "QPF1H": TimeSeriesDefaults(
        name="QPF1H", var=lambda ds: ds["PREC_ACC_C"] + ds["PREC_ACC_NC"]),
    "QPFT": TimeSeriesDefaults(
        name="QPFT", var=lambda ds: ds["RAINNC"] + ds["RAINC"]),
    "RH2M": TimeSeriesDefaults(name="RH2M", var=lambda ds: ds["rh2m"]),
    "SWDOWN": TimeSeriesDefaults(name="SWDOWN", var=lambda ds: ds["SWDOWN"]),
    "T2M": TimeSeriesDefaults(name="T2M", var=lambda ds: ds["T2"] - 273.15),
    "WIND10M": TimeSeriesDefaults(
        name="WIND10M",
        var=lambda ds: wrf_calc.wind_speed(
            ds["umet10"] * units("m/s"), ds["vmet10"] * units("m/s"))),
    "WIND10MDIR": TimeSeriesDefaults(
        name="WIND10MDIR",
        var=lambda ds: wrf_calc.wind_direction(
            ds["umet10"] * units("m/s"), ds["vmet10"] * units("m/s"))),
}
