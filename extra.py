import config
import xarray as xr
import metpy
from netCDF4 import Dataset
import wrf_calc


class HashDataset:
    def __init__(self, ds):
        self.ds = ds

    def __hash__(self):
        return id(self.ds)

def read_data_mpas(init_time, i):
    """
    Reads MPAS output data.

    Parameters
    ----------
    init_time : datetime.datetime
        Initialization time of MPAS run.
    i : int
        Index of hour to read.
    """
    ds = xr.open_dataset(f"../mpas_out/{init_time:%Y%m%d%H}/latlon.nc")
    return ds.isel({"Time": i}).sel(
        {"latitude": (ds["latitude"] > config.mpas_trim["south"]) &
                     (ds["latitude"] < config.mpas_trim["north"]),
         "longitude": (ds["longitude"] > config.mpas_trim["west"]) &
                      (ds["longitude"] < config.mpas_trim["east"])})

def read_data_wrf(init_time, fcst_time, domain, ens, wrf_products):
    """
    Reads WRF data.

    Parameters
    ----------
    init_time : datetime.datetime
        Initialization time of WRF data.
    fcst_time : datetime.datetime
        Forecast time of WRF data.
    domain : int
        Domain of WRF data.
    ens : str
        Ensemble name.
    wrf_products : dict
        Dictionary of WRF products.

    Returns
    -------
    HashDataset
        Dataset of WRF data.
    """
    filename = f"../wrfout/{init_time:%Y%m%d%H}/{ens}/wrfout24_d{domain:02}_{fcst_time:%Y-%m-%d_%H_%M_%S}"
    ds_wrf = HashDataset(xr.open_dataset(filename).squeeze())
    nc_wrf = Dataset(filename)

    if True in [key in wrf_products for key in ["CAPEMU", "CIN", "LCL", "LFC"]]:
        ds_wrf.ds["capemu"], ds_wrf.ds["cin"], ds_wrf.ds["lcl"], ds_wrf.ds["lfc"] = wrf_calc.cape_2d(
            nc_wrf)
    if True in [key in wrf_products for key in [
            "CLOUDCOVERHIGH", "CLOUDCOVERLOW", "CLOUDCOVERMID"]]:
        ds_wrf.ds["cloud_cover_high"], ds_wrf.ds["cloud_cover_low"], ds_wrf.ds["cloud_cover_mid"] = wrf_calc.cloud_cover(
            nc_wrf)
    if True in [key in wrf_products for key in ["CTT"]]:
        ds_wrf.ds["ctt"] = wrf_calc.ctt(nc_wrf)
    if True in [key in wrf_products for key in ["DBZ", "DBZ1KM"]]:
        ds_wrf.ds["dbz"] = wrf_calc.dbz(nc_wrf)
    if True in [key in wrf_products for key in ["DPT2M"]]:
        ds_wrf.ds["dpt2m"] = wrf_calc.dpt2m(ds_wrf)
    if True in [key in wrf_products for key in [
        "DBZ1KM", "OMEGA500", "OMEGA700", "OMEGA850", "LR700500", "PWAT",
        "RH250", "RH500", "RH700", "RH850", "RH925",
        "T500", "T700", "T850", "T925", "VORT500", "VORT700", "VORT850",
            "WIND250", "WIND500", "WIND700", "WIND850", "WIND925"]]:
        ds_wrf.ds["hgt"] = wrf_calc.height(ds_wrf)
    if True in [key in wrf_products for key in [
            "OMEGA500", "OMEGA700", "OMEGA850"]]:
        ds_wrf.ds["omega"] = wrf_calc.omega(nc_wrf)
    if True in [key in wrf_products for key in [
        "CAPEMU", "LR700500", "OMEGA500", "OMEGA700", "OMEGA850", "PWAT",
        "RH250", "RH500", "RH700", "RH850", "RH925",
        "T500", "T700", "T850", "T925", "VORT500", "VORT700", "VORT850",
            "WIND250", "WIND500", "WIND700", "WIND850", "WIND925"]]:
        ds_wrf.ds["prs"] = wrf_calc.pressure(ds_wrf)
    if True in [key in wrf_products for key in ["PWAT"]]:
        ds_wrf.ds["pwat"] = wrf_calc.pwat(ds_wrf)
    if True in [key in wrf_products for key in [
            "RH250", "RH500", "RH700", "RH850", "RH925"]]:
        ds_wrf.ds["rh"] = wrf_calc.relative_humidity(ds_wrf)
    if True in [key in wrf_products for key in ["RH2M"]]:
        ds_wrf.ds["rh2m"] = wrf_calc.rh2m(ds_wrf)
    if True in [key in wrf_products for key in [
            "DBZ", "DBZ1KM", "DPT2M", "RH2M", "IVT", "T2M", "TA2M", "WIND10M"]]:
        ds_wrf.ds["slp"] = wrf_calc.slp(nc_wrf)
    if True in [key in wrf_products for key in [
            "LR700500", "T500", "T700", "T850", "T925"]]:
        ds_wrf.ds["tmp"] = wrf_calc.temperature(ds_wrf)
    if True in [key in wrf_products for key in [
        "CAPEMU", "IVT", "OMEGA500", "OMEGA700", "OMEGA850", "PWAT",
        "RH250", "RH500", "RH700", "RH850", "RH925",
        "T500", "T700", "T850", "T925", "VORT500", "VORT700", "VORT850",
            "WIND250", "WIND500", "WIND700", "WIND850", "WIND925"]]:
        ds_wrf.ds["umet"], ds_wrf.ds["vmet"] = wrf_calc.uvmet(nc_wrf)
    if True in [key in wrf_products for key in [
        "CAPEMU", "DBZ", "DBZ1KM", "DPT2M", "RH2M", "T2M", "TA2M", "WIND10M"]]:
        ds_wrf.ds["umet10"], ds_wrf.ds["vmet10"] = wrf_calc.uvmet10(nc_wrf)

    if True in [key in wrf_products for key in ["IVT"]]:
        ds_wrf.ds["ivt"], ds_wrf.ds["ivt_u"], ds_wrf.ds["ivt_v"] = wrf_calc.ivt(
            ds_wrf)
    if True in [key in wrf_products for key in ["TA2M"]]:
        ds_wrf.ds["ta2"] = wrf_calc.apparent_temperature(ds_wrf)
    if True in [key in wrf_products for key in [
            "VORT500", "VORT700", "VORT850"]]:
        ds_wrf.ds["vort"] = wrf_calc.vorticity(ds_wrf)

    return ds_wrf

def read_data_wrf_sounding(init_time, fcst_time, domain, ens):
    """
    Reads WRF data for soundings.

    Parameters
    ----------
    init_time : datetime.datetime
        Initialization time of WRF data.
    fcst_time : datetime.datetime
        Forecast time of WRF data.
    domain : int
        Domain of WRF data.
    ens : str
        Ensemble name.
    wrf_products : dict
        Dictionary of WRF products.

    Returns
    -------
    HashDataset
        Dataset of WRF data.
    """
    filename = f"../wrfout/{init_time:%Y%m%d%H}/{ens}/wrfout24_d{domain:02}_{fcst_time:%Y-%m-%d_%H_%M_%S}"
    ds_wrf = HashDataset(xr.open_dataset(filename).squeeze())
    nc_wrf = Dataset(filename)
    ds_wrf.ds["omega"] = wrf_calc.omega(nc_wrf)
    ds_wrf.ds["umet"], ds_wrf.ds["vmet"] = wrf_calc.uvmet(nc_wrf)
    ds_wrf.ds = ds_wrf.ds.drop_dims(["west_east_stag", "south_north_stag"])
    return ds_wrf

def add_cf_mpas(ds):
    """
    Adds projection data to MPAS dataset.

    Parameters
    ----------
    ds : xarray.dataset
        Dataset without projection data.

    Returns
    -------
    xarray.dataset
        Dataset with projection data.
    """
    return ds.metpy.parse_cf()

def add_cf_wrf(ds):
    """
    Adds projection data to WRF dataset.

    Parameters
    ----------
    ds : xarray.dataset
        Dataset without projection data.

    Returns
    -------
    ds : xarray.dataset
        Dataset with projection data.
    """
    if ds.MAP_PROJ == 1:
        ds = ds.metpy.assign_crs({
            "grid_mapping_name": "lambert_conformal_conic",
            "earth_radius": 6370000,
            "standard_parallel": [ds.attrs["TRUELAT1"], ds.attrs["TRUELAT2"]],
            "longitude_of_central_meridian": ds.attrs["STAND_LON"],
            "latitude_of_projection_origin": ds.attrs["MOAD_CEN_LAT"]})
    else:
        raise NotImplementedError("Projection not implemented.")
    return ds

def error_handler(e):
    """
    Error handler for multiprocessing.
    """
    print("-->{}<--".format(e.__cause__), flush=True)
