from datetime import datetime, timedelta
import pandas as pd
import numpy as np
import config
from regions_config import wrf_regions
from products import wrf_products
from extra import read_data_wrf, add_cf_wrf, error_handler
from mapper import SingleMap
from ccmaps import mask_cmap
import multiprocessing as mp
import matplotlib
matplotlib.rcParams["font.family"] = ["sans serif"]
import argparse


parser = argparse.ArgumentParser(description="Plots images of WRF output.")
parser.add_argument(
    "--date", type=int, required=True,
    help="Starting date of WRF run. (YYYYMMDDHH)")
parser.add_argument(
    "--hours", type=int, required=True, help="Total hours to plot data.")
parser.add_argument(
    "--ens", type=str, default=False, help="Ensemble member.")
args = parser.parse_args()

start_date = datetime.strptime(str(args.date), "%Y%m%d%H")
end_date = datetime.strptime(str(args.date), "%Y%m%d%H") + timedelta(hours=args.hours)

stations = pd.read_csv("stations.txt")
stations = stations[
    (stations["LAT"] > config.wrf_trim["south"]) &
    (stations["LAT"] < config.wrf_trim["north"]) &
    (stations["LON"] > config.wrf_trim["west"]) &
    (stations["LON"] < config.wrf_trim["east"])]

def plot_map(projection, extent, barb_skip,
             lon, lat,
             contourf_var, contourf_kwargs, cmap_label,
             mask_var,
             barb_u, barb_v, barb_kwargs,
             barb_u_2, barb_v_2, barb_kwargs_2,
             barb_u_3, barb_v_3, barb_kwargs_3,
             quiver_u, quiver_v, quiver_kwargs,
             contour_var, contour_kwargs,
             contour_var_2, contour_kwargs_2,
             title, init_time, fcst_time, 
             shapefile_kwargs,
             stations, station_kwargs,
             hilo_var, hilo_kwargs,
             region, filename, ens,
             ):
    map_img = SingleMap(projection=projection)
    map_img.set_extent(extent)

    map_img.draw_contourf(
        lon, lat, contourf_var, cmap_label=cmap_label, **contourf_kwargs)
    if mask_var is not False:
        map_img.draw_pcolormesh(
            lon, lat, mask_var * 0, cmap=mask_cmap(), plot_cbar=False)
    if barb_u is not False:
        map_img.draw_barbs(
            lon[::barb_skip, ::barb_skip].values,
            lat[::barb_skip, ::barb_skip].values,
            barb_u[::barb_skip, ::barb_skip].values,
            barb_v[::barb_skip, ::barb_skip].values,
            **barb_kwargs)
    if barb_u_2 is not False:
        map_img.draw_barbs(
            lon[::barb_skip, ::barb_skip].values,
            lat[::barb_skip, ::barb_skip].values,
            barb_u_2[::barb_skip, ::barb_skip].values,
            barb_v_2[::barb_skip, ::barb_skip].values,
            **barb_kwargs_2)
    if barb_u_3 is not False:
        map_img.draw_barbs(
            lon[::barb_skip, ::barb_skip].values,
            lat[::barb_skip, ::barb_skip].values,
            barb_u_3[::barb_skip, ::barb_skip].values,
            barb_v_3[::barb_skip, ::barb_skip].values,
            **barb_kwargs_3)
    if quiver_u is not False:
        map_img.draw_quiver(
            lon[::barb_skip, ::barb_skip].values,
            lat[::barb_skip, ::barb_skip].values,
            quiver_u[::barb_skip, ::barb_skip].values,
            quiver_v[::barb_skip, ::barb_skip].values,
            **quiver_kwargs)
    if contour_var is not False:
        map_img.draw_contour(lon, lat, contour_var, **contour_kwargs)
    if contour_var_2 is not False:
        map_img.draw_contour(lon, lat, contour_var_2, **contour_kwargs_2)
    if stations is not False:
        map_img.draw_station_values_wrf(
            contourf_var, stations, **station_kwargs)
    if hilo_var is not False:
        map_img.draw_hilo_wrf(hilo_var, **hilo_kwargs)

    map_img.draw_title(
        f"{title}\nInit: {init_time:%Y-%m-%d %HZ} | FcstHr: [{int((fcst_time - init_time).total_seconds() / 3600)}] | Valid: {fcst_time:%Y-%m-%d %H:%MZ}",
        loc="left")
    map_img.draw_shapefiles(**shapefile_kwargs)

    map_img.save_image(
        f"../images_wrf/{init_time:%Y%m%d%H}/{ens}/{filename.lower()}/{region.lower()}",
        f"{filename.lower()}_{region.lower()}_{fcst_time:%Y%m%d%H%M}")


def plot_hour(init_time, fcst_time, domain, ens):
    print("Reading data for", fcst_time, "domain", domain, flush=True)
    ds_wrf = read_data_wrf(init_time, fcst_time, domain, ens, wrf_products(domain))
    ds_wrf.ds = add_cf_wrf(ds_wrf.ds)
    for product in wrf_products(domain).values():
        print(fcst_time, product.filename, flush=True)
        for region in wrf_regions[domain].values():
            if region.scale <= product.scales:
                plot_map(
                    projection=region.projection,
                    extent=region.extent,
                    barb_skip=region.barb_skip,
                    lon=ds_wrf.ds["XLONG"],
                    lat=ds_wrf.ds["XLAT"],
                    contourf_var=product.contourf_var(ds_wrf),
                    contourf_kwargs=product.contourf_kwargs,
                    mask_var=product.mask_var(ds_wrf) if product.mask_var else False,
                    cmap_label=product.cmap_label,
                    barb_u=product.barb_u(ds_wrf) if product.barb_u else False,
                    barb_v=product.barb_v(ds_wrf) if product.barb_v else False,
                    barb_kwargs=product.barb_kwargs,
                    barb_u_2=product.barb_u_2(ds_wrf) if product.barb_u_2 else False,
                    barb_v_2=product.barb_v_2(ds_wrf) if product.barb_v_2 else False,
                    barb_kwargs_2=product.barb_kwargs_2,
                    barb_u_3=product.barb_u_3(ds_wrf) if product.barb_u_3 else False,
                    barb_v_3=product.barb_v_3(ds_wrf) if product.barb_v_3 else False,
                    barb_kwargs_3=product.barb_kwargs_3,
                    quiver_u=product.quiver_u(ds_wrf) if product.quiver_u else False,
                    quiver_v=product.quiver_v(ds_wrf) if product.quiver_v else False,
                    quiver_kwargs=product.quiver_kwargs,
                    contour_var=product.contour_var(ds_wrf) if product.contour_var else False,
                    contour_kwargs=product.contour_kwargs,
                    contour_var_2=product.contour_var_2(ds_wrf) if product.contour_var_2 else False,
                    contour_kwargs_2=product.contour_kwargs_2,
                    title=product.title,
                    init_time=init_time,
                    fcst_time=fcst_time,
                    shapefile_kwargs=product.shapefile_kwargs,
                    stations=stations if product.stations else False,
                    station_kwargs={
                        "priority": region.station_priority,
                        **product.station_kwargs},
                    hilo_var=product.hilo_var(ds_wrf) if product.hilo_var else False,
                    hilo_kwargs={
                        "size": region.hilo_size,
                        **product.hilo_kwargs},
                    region=region.name,
                    filename=product.filename,
                    ens=ens,
                )


def main():
    with mp.Pool(processes=24, maxtasksperchild=1) as pool:
        for domain in np.arange(
                config.wrf_max_domain, config.wrf_min_domain - 1, -1):
            for time in pd.date_range(
                    start_date, end_date, freq=config.wrf_freq[domain]):
                pool.apply_async(
                    plot_hour,
                    args=(start_date, time, domain, args.ens, ),
                    error_callback=error_handler)
        pool.close()
        pool.join()


if __name__ == "__main__":
    mp.set_start_method("fork")
    main()
