import metpy.calc as mpcalc
from extra import HashDataset
from functools import cache
import wrf
from metpy.units import units
import scipy.interpolate as interpolate


def need_cache(func):
    func_cache = cache(lambda ds, *args: func(ds, hash=True, *args))

    def decorated(*args):
        if isinstance(args[0], HashDataset):
            return func_cache(*args)
        else:
            return func(*args)
    return decorated


def advection(scalar, u, v):
    """
    Calculates advection.

    Parameters
    ----------
    scalar : xarray.DataArray
        The quantity to be advected.
    u : xarray.DataArray
        The wind component in the x dimension.
    v : xarray.DataArray
        The wind component in the y dimension.

    Returns
    -------
    xarray.DataArray
        DataArray with advection.
    """
    if hash:
        return mpcalc.advection(scalar, u, v)

@need_cache
def apparent_temperature(ds, hash=False):
    """
    Calculates apparent temperature.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset with WRF data.

    Returns
    -------
    xarray.DataArray
        DataArray with apparent temperature.
    """
    if hash:
        return mpcalc.apparent_temperature(
            ds.ds["T2"], rh2m(ds), wind_speed(
                ds.ds["umet10"], ds.ds["vmet10"]),
            mask_undefined=False).metpy.dequantify()
    else:
        return mpcalc.apparent_temperature(
            ds["T2"], rh2m(ds), wind_speed(ds["umet10"], ds["vmet10"]),
            mask_undefined=False).metpy.dequantify()


@need_cache
def area_density(ds, hash=False):
    """
    Calculates mass per area of mass grid.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset with WRF data.

    Returns
    -------
    xarray.DataArray
        DataArray with grid mass per area.
    """
    if hash:
        return delta_z(ds) / ds.ds["ALT"]
    else:
        return delta_z(ds) / ds["ALT"]


def cape_2d(ds):
    """
    Calculates MUCAE, MCIN, LCL, and LFC.

    Parameters
    ----------
    ds : netCDF4.Dataset
        Dataset with WRF data.

    Returns
    -------
    xarray.DataArray
        DataArray with MUCAPE data.
    xarray.DataArray
        DataArray with MCIN data.
    xarray.DataArray
        DataArray with LCL data.
    xarray.DataArray
        DataArray with LFC data.
    """
    out = wrf.getvar(ds, "cape_2d")
    return out.sel({"mcape_mcin_lcl_lfc": "mcape"}), out.sel({"mcape_mcin_lcl_lfc": "mcin"}), out.sel({"mcape_mcin_lcl_lfc": "lcl"}), out.sel({"mcape_mcin_lcl_lfc": "lfc"}),


def cloud_cover(ds):
    """
    Calculates low, mid, and high cloud cover.

    Parameters
    ----------
    ds : netCDF4.Dataset
        Dataset with WRF data.

    Returns
    -------
    xarray.DataArray
        DataArray with high cloud cover data.
    xarray.DataArray
        DataArray with low cloud cover data.
    xarray.DataArray
        DataArray with mid cloud cover data.
    """
    cloud_cover = wrf.getvar(ds, "cloudfrac")
    return cloud_cover.sel({"low_mid_high": "high"}), cloud_cover.sel({"low_mid_high": "low"}), cloud_cover.sel({"low_mid_high": "mid"})


def ctt(ds):
    """
    Calculates cloud top temperature.

    Parameters
    ----------
    ds : netCDF4.Dataset
        Dataset with WRF data.

    Returns
    -------
    xarray.DataArray
        DataArray with cloud top temperature data.
    """
    return wrf.getvar(ds, "ctt")


def dbz(ds):
    """
    Calculates simulated reflectivity.

    Parameters
    ----------
    ds : netCDF4.Dataset
        Dataset with WRF data.

    Returns
    -------
    xarray.DataArray
        DataArray with simulated reflectivity data.
    """
    return wrf.getvar(ds, "dbz")


@need_cache
def delta_z(ds, hash=False):
    """
    Calculates the depth of mass grid.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset with WRF data.

    Returns
    -------
    xarray.DataArray
        DataArray with grid depth.
    """
    if hash:
        return ((ds.ds["PH"][1:] + ds.ds["PHB"][1:]) / 9.81 -
                (ds.ds["PH"][:-1] + ds.ds["PHB"][:-1]) / 9.81
                ).rename({"bottom_top_stag": "bottom_top"})
    else:
        return ((ds["PH"][1:] + ds["PHB"][1:]) / 9.81 -
                (ds["PH"][:-1] + ds["PHB"][:-1]) / 9.81
                ).rename({"bottom_top_stag": "bottom_top"})


@need_cache
def dewpoint(ds, hash=False):
    """
    Calculates dewpoint.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset with WRF data.

    Returns
    -------
    xarray.DataArray
        DataArray with dewpoint.
    """
    if hash:
        return mpcalc.dewpoint_from_specific_humidity(
            pressure(ds) * units("pascal"), temperature(ds) * units("kelvin"),
            mpcalc.specific_humidity_from_mixing_ratio(ds.ds["QVAPOR"])).metpy.dequantify()
    else:
        return mpcalc.dewpoint_from_specific_humidity(
            pressure(ds) * units("pascal"), temperature(ds) * units("kelvin"),
            mpcalc.specific_humidity_from_mixing_ratio(ds["QVAPOR"])).metpy.dequantify()


@need_cache
def dpt2m(ds, hash=False):
    """
    Calculates dew point at 2m.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset with WRF data.

    Returns
    -------
    xarray.DataArray
        DataArray with dew point.
    """
    if hash:
        return mpcalc.dewpoint_from_specific_humidity(
            ds.ds["PSFC"], ds.ds["T2"], ds.ds["Q2"]).metpy.dequantify()
    else:
        return mpcalc.dewpoint_from_specific_humidity(
            ds["PSFC"], ds["T2"], ds["Q2"]).metpy.dequantify()


def frontogenesis(theta, u, v):
    """
    Calculated frontogenesis.

    Parameters
    ----------
    theta : xarray.DataArray
        Potential temperature.
    u : xarray.DataArray
        The wind component in the x dimension.
    v : xarray.DataArray
        The wind component in the y dimension.

    Returns
    -------
    xarray.DataArray
        DataArray with advection.
    """
    return mpcalc.frontogenesis(theta, u, v)
    

@need_cache
def height(ds, hash=False):
    """
    Calculates height of mass grid.
    Values checked against wrf-python.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset with WRF data.

    Returns
    -------
    xarray.DataArray
        DataArray with heights.
    """
    if hash:
        full_mass = (ds.ds["PH"] + ds.ds["PHB"]) / 9.81
    else:
        full_mass = (ds["PH"] + ds["PHB"]) / 9.81
    return ((full_mass[1:] + full_mass[:-1]) / 2
            ).rename({"bottom_top_stag": "bottom_top"})


def interplevel(field3d, vert, desiredlev):
    """
    Return the three-dimensional field interpolated to a horizontal plane at the
    specified vertical level.

    Parameters
    ----------
    field3d : xarray.DataArray or numpy.ndarray
        A three-dimensional field to interpolate, with the rightmost dimensions
        of nz x ny x nx.
    vert : xarray.DataArray or numpy.ndarray 
        A three-dimensional array for the vertical coordinate, typically
        pressure or height. This array must have the same dimensionality as
        field3d.
    desiredlev : float, 1D sequence, or numpy.ndarray
        The desired vertical level(s). This can be a single value (e.g. 500), a
        sequence of values (e.g. [1000, 850, 700, 500, 250]), or a
        multidimensional array where the right two dimensions (ny x nx) must
        match field3d, and any leftmost dimensions match field3d.shape[:-3]
        (e.g. planetary boundary layer). Must be in the same units as the vert
        parameter.
    """
    return wrf.interplevel(field3d, vert, desiredlev)


def interpolate_1d(x, y, x_new, bounds_error=False, **kwargs):
    """
    Interpolate a 1-D function. x and y are arrays of values used to approximate
    some function f: y = f(x). This class returns a function whose call method
    uses interpolation to find the value of new points.

    Parameters
    ----------
    x : array-like
        A 1-D array of real values.
    y : array-like
        A N-D array of real values. The length of y along the interpolation axis
        must be equal to the length of x.
    x_new : array-like
        List of x values to interpolate to.
    bounds_error : bool (default=False)
        If True, a ValueError is raised any time interpolation is attempted on a
        value outside of the range of x (where extrapolation is necessary). If
        False, out of bounds values are assigned fill_value. By default, an
        error is raised unless fill_value="extrapolate".
    """
    f = interpolate.interp1d(x, y, bounds_error=bounds_error, **kwargs)
    return f(x_new)


@need_cache
def interp_ds(ds, desired_field, vertical_field, level, hash=False):
    """
    Interpolates 3 dimensional field to desired level.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset with WRF data.
    desired_field : str
        Variable to interpolate to.
    vertical_field : str
        Variable with vertical field.
    level : str
        Level to interpolate to.


    Returns
    -------
    xarray.DataArray
        DataArray with values at desired level
    """
    if hash:
        return interplevel(ds.ds[desired_field], ds.ds[vertical_field], level)
    else:
        return interplevel(ds[desired_field], ds[vertical_field], level)


@need_cache
def ivt(ds, hash=False):
    """
    Calculates integrated vapor transport.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset with WRF data.

    Returns
    -------
    xarray.DataArray
        DataArray with integrated vapor transport.
    xarray.DataArray
        DataArray with integrated vapor transport zonal vector.
    xarray.DataArray
        DataArray with integrated vapor transport meridional vector.
    """
    if hash:
        return (area_density(ds) * ds.ds["QVAPOR"] * wind_speed(ds.ds["umet"], ds.ds["vmet"])).sum(dim="bottom_top").metpy.dequantify(), (area_density(ds) * ds.ds["QVAPOR"] * ds.ds["umet"]).sum(dim="bottom_top"), (area_density(ds) * ds.ds["QVAPOR"] * ds.ds["vmet"]).sum(dim="bottom_top")
    else:
        return (area_density(ds) * ds["QVAPOR"] * wind_speed(ds["umet"], ds["vmet"])).sum(dim="bottom_top").metpy.dequantify(), (area_density(ds) * ds["QVAPOR"] * ds["umet"]).sum(dim="bottom_top"), (area_density(ds) * ds["QVAPOR"] * ds["vmet"]).sum(dim="bottom_top")


def omega(ds):
    """
    Calculates omega.

    Parameters
    ----------
    ds : netCDF4.Dataset
        Dataset with WRF data.

    Returns
    -------
    xarray.DataArray
        DataArray with omega.
    """
    return wrf.getvar(ds, "omega")


@need_cache
def potential_temperature(ds, hash=False):
    """
    Calculates potential temperature.
    Values checked against wrf-python.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset with WRF data.

    Returns
    -------
    xarray.DataArray
        DataArray with potential temperature.
    """
    if hash:
        return ds.ds["T"] + 300
    else:
        return ds["T"] + 300


@need_cache
def pressure(ds, hash=False):
    """
    Calculates pressure on mass grid.
    Values checked against wrf-python.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset with WRF data.

    Returns
    -------
    xarray.DataArray
        DataArray with perssure.
    """
    if hash:
        return ds.ds["P"] + ds.ds["PB"]
    else:
        return ds["P"] + ds["PB"]


@need_cache
def pwat(ds, hash=False):
    """
    Calculates precipitable water.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset with WRF data.

    Returns
    -------
    xarray.DataArray
        DataArray with precipitable water.
    """
    if hash:
        return (area_density(ds) * ds.ds["QVAPOR"]).sum(dim="bottom_top")
    else:
        return (area_density(ds) * ds["QVAPOR"]).sum(dim="bottom_top")


@need_cache
def relative_humidity(ds, hash=False):
    """
    Calculates relative humidity.
    Values checked against wrf-python.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset with WRF data.

    Returns
    -------
    xarray.DataArray
        DataArray with relative humidity.
    """
    if hash:
        return mpcalc.relative_humidity_from_mixing_ratio(
            pressure(ds) * units("pascal"), temperature(ds) * units("kelvin"),
            (ds.ds["QVAPOR"]) * units(""))
    else:
        return mpcalc.relative_humidity_from_mixing_ratio(
            pressure(ds) * units("pascal"), temperature(ds) * units("kelvin"),
            (ds["QVAPOR"]) * units(""))


@need_cache
def rh2m(ds, hash=False):
    """
    Calculates relative humidity at 2m.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset with WRF data.

    Returns
    -------
    xarray.DataArray
        DataArray with relative humidity.
    """
    if hash:
        return mpcalc.relative_humidity_from_specific_humidity(
            ds.ds["PSFC"] * units("pascal"), ds.ds["T2"] * units("kelvin"),
            ds.ds["Q2"]) * 100
    else:
        return mpcalc.relative_humidity_from_specific_humidity(
            ds["PSFC"] * units("pascal"), ds["T2"] * units("kelvin"),
            ds["Q2"]) * 100


def slp(ds):
    """
    Calculates sea level pressure.

    Parameters
    ----------
    ds : netCDF4.Dataset
        Dataset with WRF data.

    Returns
    -------
    xarray.DataArray
        DataArray with sea level pressure.
    """
    return wrf.getvar(ds, "slp")


@need_cache
def temperature(ds, hash=False):
    """
    Calculates temperature.
    Values checked against wrf-python.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset with WRF data.

    Returns
    -------
    xarray.DataArray
        DataArray with temperature.
    """
    if hash:
        return mpcalc.temperature_from_potential_temperature(
            pressure(ds) * units("pascal"),
            potential_temperature(ds) * units("kelvin")).metpy.dequantify()
    else:
        return mpcalc.temperature_from_potential_temperature(
            pressure(ds) * units("pascal"),
            potential_temperature(ds) * units("kelvin")).metpy.dequantify()


def uvmet(ds):
    """
    Calculates U and V component of wind relative to Earth coordinates.

    Parameters
    ----------
    ds : netCDF4.Dataset
        Dataset with WRF data.

    Returns
    -------
    xarray.DataArray
        DataArray with U component of wind.
    xarray.DataArray
        DataArray with V component of wind.
    """
    uvmet = wrf.getvar(ds, "uvmet")
    return uvmet.sel({"u_v": "u"}), uvmet.sel({"u_v": "v"})


def uvmet10(ds):
    """
    Calculates U and V component of 10 meter wind relative to Earth coordinates.

    Parameters
    ----------
    ds : netCDF4.Dataset
        Dataset with WRF data.

    Returns
    -------
    xarray.DataArray
        DataArray with U component of 10 meter wind.
    xarray.DataArray
        DataArray with V component of 10 meter wind.
    """
    uvmet = wrf.getvar(ds, "uvmet10")
    return uvmet.sel({"u_v": "u"}), uvmet.sel({"u_v": "v"})


def uvmet10_xr(ds):
    """
    Calculates U and V component of 10 meter wind relative to Earth coordinates
    with xarray.Dataset.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset with WRF data.

    Returns
    -------
    xarray.DataArray
        DataArray with U component of 10 meter wind.
    xarray.DataArray
        DataArray with V component of 10 meter wind.
    """
    u = ds["U10"] * ds["COSALPHA"] - ds["V10"] * ds["SINALPHA"]
    v = ds["V10"] * ds["COSALPHA"] + ds["U10"] * ds["SINALPHA"]
    return u, v


@need_cache
def vorticity(ds, hash=False):
    """
    Calculates vorticity.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset with WRF data.

    Returns
    -------
    xarray.DataArray
        DataArray with vorticity.
    """
    if hash:
        return mpcalc.vorticity(ds.ds["umet"], ds.ds["vmet"])
    else:
        return mpcalc.vorticity(ds["umet"], ds["vmet"])


def wind_direction(u, v):
    """
    Calculates wind direction from U and V components of wind.

    Parameters
    ----------
    u : xarray.DataArray
        Wind component in the X (East-West) direction.
    v : xarray.DataArray
        Wind component in the Y (North-South) direction.

    Returns
    -------
    xarray.DataArray
        DataArray with wind direction.
    """
    return mpcalc.wind_direction(u, v)


def wind_speed(u, v):
    """
    Calculates wind speed from U and V components of wind.

    Parameters
    ----------
    u : xarray.DataArray
        Wind component in the X (East-West) direction.
    v : xarray.DataArray
        Wind component in the Y (North-South) direction.

    Returns
    -------
    xarray.DataArray
        DataArray with wind speed.
    """
    return mpcalc.wind_speed(u, v)
