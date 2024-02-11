from metpy.plots import SkewT, Hodograph
from metpy.units import units
from pathlib import Path
import sharppy.sharptab.profile as profile
import sharppy.sharptab.utils as utils
import sharppy.sharptab.thermo as thermo
import sharppy.sharptab.params as params
import sharppy.sharptab.interp as interp
import sharppy.sharptab.winds as winds
import config
import wrf_calc
from extra import add_cf_wrf, error_handler, read_data_wrf_sounding
from regions_config import wrf_sounding_stations
import pandas as pd
import numpy as np
from datetime import datetime, timedelta
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
import matplotlib.patheffects as path_effects
import multiprocessing as mp
import argparse


parser = argparse.ArgumentParser(description="Plot skewt of WRF output.")
parser.add_argument(
    "--date", type=int, required=True,
    help="Starting date of WRF run. (YYYYMMDDHH)")
parser.add_argument(
    "--hours", type=int, required=True, help="Total hours to plot data.")
parser.add_argument(
    "--ens", type=str, default="", help="Ensemble member.")
args = parser.parse_args()

start_date = datetime.strptime(str(args.date), "%Y%m%d%H")
end_date = start_date + timedelta(hours=args.hours)

h_sep = 0.29
top_sep = 0.95
skewt_r = 0.44
barb_r = 0.49
tadv_r = 0.54
skewt_shift = 0.02
right = 0.98

class Skewt():
    def __init__(self, figsize=[12, 8]):
        self.fig = plt.figure(figsize=figsize)
        self.skewt = SkewT(
            self.fig, aspect=100,
            rect=(0 + skewt_shift, h_sep, skewt_r, top_sep - h_sep))

        self.skewt.plot_dry_adiabats(linewidth=0.5)
        self.skewt.plot_mixing_lines(linewidth=0.5)
        self.skewt.plot_moist_adiabats(linewidth=0.5)
        self.skewt.ax.axvline(-20, color="blue", linewidth=0.5, linestyle="--")
        self.skewt.ax.axvline(0, color="blue", linewidth=0.5, linestyle="--")
        for i in range(-100, 50, 20):
            self.skewt.ax.axvspan(i, i + 10, color="black", alpha=0.05)

        self.omega_ax = plt.axes(
            (0.08, h_sep, 0.05, top_sep - h_sep), sharey=self.skewt.ax,
            clip_on=False)
        self.omega_ax.set_xlim(2, -6)
        self.omega_ax.set_xticks([2, -2, -4])
        self.omega_ax.tick_params(
            axis="x", direction="in", pad=-10, colors="magenta", labelsize=6)
        self.omega_ax.axvline(0, color="magenta", linestyle="--", linewidth=0.5)
        self.omega_ax.spines["left"].set_visible(False)
        self.omega_ax.spines["right"].set_visible(False)
        self.omega_ax.yaxis.set_visible(False)
        self.omega_ax.patch.set_alpha(0)
        #self.omega_ax.axis("off")

        self.barb_ax = plt.axes(
            (skewt_r, h_sep, barb_r - skewt_r, top_sep - h_sep),
            sharey=self.skewt.ax, clip_on=False)
        self.barb_ax.axis("off")

        self.tadv_ax = plt.axes(
            (barb_r, h_sep, tadv_r - barb_r, top_sep - h_sep),
            sharey=self.skewt.ax)
        self.tadv_ax.axis("off")
        self.tadv_ax.axvline(0, color="grey", linestyle="--", linewidth=0.5)

        self.hodo_ax = plt.axes(
            (tadv_r, h_sep, right - tadv_r, top_sep - h_sep))
        self.hodo = Hodograph(self.hodo_ax)
        self.hodo.add_grid(increment=10, linestyle="--", linewidth=0.5)
        for i in range(10, 161, 10):
            self.hodo.ax.text(
                i, 0, f"\n{i}", ha="center", va="center", color="grey",
                fontsize=6, clip_box=self.hodo_ax, clip_on=True)
            self.hodo.ax.text(
                -i, 0, f"\n{i}", ha="center", va="center", color="grey",
                fontsize=6, clip_box=self.hodo_ax, clip_on=True)
            self.hodo.ax.text(
                0, i, f" {i}", ha="left", va="center", color="grey",
                fontsize=6, clip_box=self.hodo_ax, clip_on=True)
            self.hodo.ax.text(
                0, -i, f" {i}", ha="left", va="center", color="grey",
                fontsize=6, clip_box=self.hodo_ax, clip_on=True)
        self.hodo.ax.set_xticks([])
        self.hodo.ax.set_yticks([])
        self.hodo_ax.set_adjustable("datalim")

        thetae_left = 0.8
        thetae_bottom = 0.02
        self.thetae_ax = plt.axes(
            (thetae_left, thetae_bottom, right - thetae_left,
                h_sep - thetae_bottom))
        self.thetae_ax.tick_params(axis="x", direction="in", pad=-15)
        self.thetae_ax.tick_params(axis="y", direction="in", pad=-25)
        self.thetae_ax.set_ylim(999, 501)
        self.thetae_ax.text(
            0.02, 0.98, "Θ-e", ha="left", va="top", fontsize=10,
            transform=self.thetae_ax.transAxes)
    
    def plot_barbs(self, p, u, v):
        """
        Plots wind barbs on side of skewt.

        Parameters
        ----------
        p : list
            Pressure values in millibars.
        u : list
            U component of wind in knots.
        v : list
            V component of wind in knots.
        """
        p_barb = [
            i for i in np.arange(1000, 99, -50) if (i > p[-1]) & (i < p[0])]
        u_barb = wrf_calc.interpolate_1d(np.log(p), u, np.log(p_barb))
        v_barb = wrf_calc.interpolate_1d(np.log(p), v, np.log(p_barb))
        self.barb_ax.barbs(
            [0] * len(p_barb), p_barb, u_barb, v_barb, length=6, linewidth=1,
            clip_on=False)

    def plot_dgz(self, b_prs, t_prs):
        """
        Plots DGZ layer on skewt.

        Parameters
        ----------
        b_prs : float
            Bottom pressure of inflow layer in millibars.
        t_prs : float
            Top pressure of inflow layer in millibars.
        """
        self.omega_ax.axhline(
            b_prs, color="cyan", linewidth=3, linestyle="--")
        self.omega_ax.axhline(
            t_prs, color="cyan", linewidth=3, linestyle="--")

    def plot_heights(self, p, h):
        """
        Plots height markers on skewt.

        Parameters
        ----------
        p : list
            Pressure values in millibars.
        h : list
            Height values in meters.
        """
        hgts = [i for i in [1, 3, 6, 9, 12, 15] if (i < max(h) / 1000)]
        h_prs = wrf_calc.interpolate_1d(h / 1000, p, hgts)
        trans = transforms.blended_transform_factory(
            self.skewt.ax.transAxes, self.skewt.ax.transData)
        self.skewt.ax.text(
            0.01, p.max(), f"SFC", color="red", ha="left", va="center",
            fontsize=8, transform=trans)
        for prs, hgt in zip(h_prs, hgts):
            self.skewt.ax.text(
                0.01, prs, f"{hgt} km", color="red", ha="left", va="center",
                fontsize=8, transform=trans)

    def plot_hodo(self, h, u, v):
        """
        Plot hodograph line.

        Parameters
        ----------
        h : list
            Height values in meters.
        u : list
            U component of wind in knots.
        v : list
            V component of wind in knots.
        """
        mask = h < 12000
        h = h[mask]
        u = u[mask]
        v = v[mask]
        u_min = min(u)
        u_max = max(u)
        v_min = min(v)
        v_max = max(v)
        u_range = u_max - u_min
        v_range = v_max - v_min
        u_min_dis = u_min - u_range * 0.2
        u_max_dis = u_max + u_range * 0.2
        v_min_dis = v_min - v_range * 0.2
        v_max_dis = v_max + v_range * 0.2
        if u_min_dis > 0:
            u_min_dis = -5
        if v_min_dis > 0:
            v_min_dis = -5
        if u_max_dis < 0:
            u_max_dis = 5
        if v_max_dis < 0:
            v_max_dis = 5
        self.hodo.ax.set_xlim(u_min_dis, u_max_dis)
        self.hodo.ax.set_ylim(v_min_dis, v_max_dis)

        self.hodo.plot_colormapped(
            u, v, h / 1000, intervals=np.array([0, 3, 6, 9, 99]) * units(""),
            colors=["red", "green", "orange", "blue"], linewidth=3,
            capstyle="round", zorder=7)

        h_label = [
            i for i in [1, 2, 3, 4, 5, 6, 7, 8, 9] if (i > h[0] / 1000) & (i < h[-1] / 1000)]
        u_label = wrf_calc.interpolate_1d(h / 1000, u, h_label)
        v_label = wrf_calc.interpolate_1d(h / 1000, v, h_label)
        for h, x, y in zip(h_label, u_label, v_label):
            txt = self.hodo.ax.text(
                x, y, f"{h}", ha="center", va="center", zorder=9)
            txt.set_path_effects(
                [path_effects.Stroke(linewidth=3, foreground="white"),
                    path_effects.Normal()])

    def plot_hodo_inflow(self, bottom_u, bottom_v, top_u, top_v, movement_u,
                        movement_v):
        """
        Plots inflow layer on hodograph.

        Parameters
        ----------
        bottom_u : float
            U component of wind at bottom of inflow layer.
        bottom_v : float
            V component of wind at bottom of inflow layer.
        top_u : float
            U component of wind at top of inflow layer.
        top_v : float
            V component of wind at top of inflow layer.
        movement_u : float
            U component of storm motion.
        movement_v : float
            V component of storm motion.
        """
        self.hodo.ax.plot(
            [bottom_u, movement_u, top_u], [bottom_v, movement_v, top_v],
            color="cyan", zorder=6)

    def plot_inflow(self, b_prs, t_prs, b_hgt, t_hgt, srh):
        """
        Plots inflow layer on skewt.

        Parameters
        ----------
        b_prs : float
            Bottom pressure of inflow layer in millibars.
        t_prs : float
            Top pressure of inflow layer in millibars.
        b_hgt : float
            Bottom height of inflow layer in meters.
        t_hgt : float
            Top height of inflow layer in meters.
        """
        x_pos = 0.35
        trans = transforms.blended_transform_factory(
            self.skewt.ax.transAxes, self.skewt.ax.transData)
        self.skewt.ax.plot(
            [x_pos - 0.025, x_pos + 0.025], [t_prs, t_prs], color="purple",
            transform=trans)
        self.skewt.ax.plot(
            [x_pos, x_pos], [b_prs, t_prs], color="purple", transform=trans)
        self.skewt.ax.plot(
            [x_pos - 0.025, x_pos + 0.025], [b_prs, b_prs], color="purple",
            transform=trans)
        self.skewt.ax.text(
            x_pos - 0.05, t_prs, f"{utils.INT2STR(t_hgt)} m   ", fontsize=8,
            ha="right", va="bottom", color="purple", transform=trans)
        self.skewt.ax.text(
            x_pos - 0.05, b_prs, f"{utils.INT2STR(b_hgt)} m   ",
            fontsize=8, ha="right", va="top", color="purple", transform=trans)
        self.skewt.ax.text(
            x_pos, t_prs, f"{utils.INT2STR(srh)} m$^{2}$ s$^{{{-2}}}$",
            fontsize=8, ha="center", va="bottom", color="purple",
            transform=trans)

    def plot_line_skewt(self, p, t, **kwargs):
        """
        Plots line on skewt.

        Parameters
        ----------
        p : list
            Pressure values in millibars.
        t : list
            Temperature values in celsius.
        """
        self.skewt.plot(p, t, **kwargs)

    def plot_max_lr(self, b_prs, t_prs, lr):
        """
        Plots maximum lapse rate range on skewt.

        Parameters
        ----------
        b_prs : float
            Bottom pressure of maximum lapse rate in millibars.
        t_prs : float
            Top pressure of maximum lapse rate in millibars.
        lr : float
            Lapse rate of layer in K/km.
        """
        x_pos = 0.775
        trans = transforms.blended_transform_factory(
            self.skewt.ax.transAxes, self.skewt.ax.transData)
        self.skewt.ax.plot(
            [x_pos - 0.025, x_pos + 0.025], [t_prs, t_prs], color="peru",
            transform=trans)
        self.skewt.ax.plot(
            [x_pos, x_pos], [b_prs, t_prs], color="peru", transform=trans)
        self.skewt.ax.plot(
            [x_pos - 0.025, x_pos + 0.025], [b_prs, b_prs], color="peru",
            transform=trans)
        self.skewt.ax.text(
            x_pos, t_prs, f"{utils.FLOAT2STR(lr, 1)}°C km$^{{{-1}}}$",
            fontsize=8, ha="center", va="bottom", color="peru", transform=trans)

    def plot_omega(self, p, omega):
        """
        Plots omega on skewt.

        Parameters
        ----------
        p : list
            Pressure values in millibars.
        omega : list
            Omega values in Pa / s
        """
        width = lambda prs, w: 10**(np.log10(prs) + w / 2) - 10**(np.log10(prs) - w / 2)
        mask_pos = (omega > 0) & (p> 100)
        mask_neg = (omega < 0) & (p> 100)
        self.omega_ax.barh(
            p[mask_pos], omega[mask_pos], color="blue",
            height=width(p[mask_pos], 0.01), clip_on=False, zorder=1)
        self.omega_ax.barh(
            p[mask_neg], omega[mask_neg], color="salmon",
            height=width(p[mask_neg], 0.01), clip_on=False, zorder=1)

    def plot_sounding_levels(self, lcl, lfc, el):
        """
        Plots LCL, LFC, and EL levels on skewt.

        Parameters
        ----------
        lcl : float
            Pressure of LCL in millibars.
        lfc : float
            Pressure of LFC in millibars.
        el : float
            Pressure of EL in millibars.
        """
        x_pos = 0.875
        trans = transforms.blended_transform_factory(
            self.skewt.ax.transAxes, self.skewt.ax.transData)
        self.skewt.ax.plot(
            [x_pos - 0.025, x_pos + 0.025], [lcl, lcl], color="lime",
            linewidth=2, transform=trans)
        self.skewt.ax.text(
            x_pos, lcl, "LCL", color="lime", ha="center", va="top",
            transform=trans)
        self.skewt.ax.plot(
            [x_pos - 0.025, x_pos + 0.025], [lfc, lfc], color="orange",
            linewidth=2, transform=trans)
        self.skewt.ax.text(
            x_pos, lfc, "LFC", color="orange", ha="center", va="bottom",
            transform=trans)
        self.skewt.ax.plot(
            [x_pos - 0.025, x_pos + 0.025], [el, el], color="purple",
            linewidth=2, transform=trans)
        self.skewt.ax.text(
            x_pos, el, "EL", color="purple", ha="center", va="bottom",
            transform=trans)

    def plot_storm_movers(self, right_u, right_v, left_u, left_v, mean_u,
                          mean_v):
        """
        Plots right, left, and mean storm motion on hodograph.

        Parameters
        ----------
        right_u : float
            U compnent of motion for right moving storm in knots.
        right_v : float
            V compnent of motion for right moving storm in knots.
        left_u : float
            U compnent of motion for left moving storm in knots.
        left_v : float
            V compnent of motion for left moving storm in knots.
        mean_u : float
            U compnent of mean motion in knots.
        mean_v : float
            V compnent of mean motion in knots.
        """
        self.hodo.ax.plot(
            [right_u, left_u], [right_v,left_v], linestyle="", marker="o",
            fillstyle="none", color="black", zorder=8)
        self.hodo.ax.plot(
            [mean_u], [mean_v], linestyle="none", marker="s", fillstyle="none",
            color="brown", zorder=8)
        self.hodo.ax.text(
            right_u, right_v, f"  RM", fontsize=8, ha="left", va="center",
            zorder=8)
        self.hodo.ax.text(
            left_u, left_v, f"  LM", fontsize=8, ha="left", va="center",
            zorder=8)
        self.hodo.ax.text(
            mean_u, mean_v, f"  MN", fontsize=8, ha="left", va="center",
            color="brown", zorder=8)
        trans = self.hodo.ax.transAxes
        right_vec = utils.comp2vec(right_u, right_v)
        left_vec = utils.comp2vec(left_u, left_v)
        mean_vec = utils.comp2vec(mean_u, mean_v)
        self.hodo.ax.text(
            0.01, 0.99,
            f"RM {utils.INT2STR(right_vec[0])}°/{utils.INT2STR(right_vec[1])}kts\nLM {utils.INT2STR(left_vec[0])}°/{utils.INT2STR(left_vec[1])}kts\nMN {utils.INT2STR(mean_vec[0])}°/{utils.INT2STR(mean_vec[1])}kts",
            fontsize=10, ha="left", va="top", zorder=10, transform=trans)

    def plot_surface_values(self, p, t, td):
        """
        Plots surface values at bottom of skewt.

        Parameters
        ----------
        p : float
            Pressure at surface in millibars.
        t : float
            Temperature as surface in celsius.
        td : float
            Dew point at surface in celsius.
        """
        self.skewt.ax.text(
            td, p, f"{utils.INT2STR(thermo.ctof(td))}°F", ha="right", va="top",
            color="green", weight="bold")
        self.skewt.ax.text(
            t, p, f"{utils.INT2STR(thermo.ctof(t))}°F", ha="left", va="top",
            color="red", weight="bold")

    def plot_temperature_advection(self, t_adv, bounds):
        """
        Plots temperature advection bar graph.

        Parameters
        ----------
        t_adv : list
            Temperature advection in C/km.
        bounds : list
            Boundaries of advection layers in millibars.
        """
        for t, bound in zip(t_adv, bounds):
            if t < 0:
                self.tadv_ax.barh(
                    (bound[1] + bound[0]) / 2, t, bound[1] - bound[0],
                    color="blue", alpha=0.25)
                self.tadv_ax.text(
                    0, (bound[1] + bound[0]) / 2, f"{round(t, 1)} ", ha="right",
                    va="center", fontsize=6)
            elif t > 0:
                self.tadv_ax.barh(
                    (bound[1] + bound[0]) / 2, t, bound[1] - bound[0],
                    color="red", alpha=0.25)
                self.tadv_ax.text(
                    0, (bound[1] + bound[0]) / 2, f" {round(t, 1)}", ha="left",
                    va="center", fontsize=6)
        x_mag = abs(max(self.tadv_ax.get_xlim()))
        self.tadv_ax.set_xlim(xmin=-x_mag, xmax=x_mag)

    def plot_thetae(self, p, thetae):
        """
        Plots theta-e profile.

        Parameters
        ----------
        p : list
            Pressure values in millibars.
        thetae : list
            Theta-e values in kelvin.
        """
        mask = p > 490
        self.thetae_ax.plot(thetae[mask], p[mask], color="blue")

    def plot_title(self, station_name, lon, lat, fcst_time, init_time):
        """
        Plots title and time into on skewt.

        Parameters
        ----------
        station_name : str
            Name of station.
        lon : float
            Longitude of station.
        lat : float
            Latitude of station.
        fcst_time : datetime.datetime
            Forecast time of skewt.
        init_time : datetime.datetime
            Initalization time of model.
        """
        self.skewt.ax.set_title(
            f"{station_name} Lon:{round(lon, 2)} Lat:{round(lat, 2)}",
            loc="left")
        self.hodo_ax.set_title(
            f"Init: {init_time:%Y-%m-%d %H:%MZ} | Valid: {fcst_time:%Y-%m-%d %H:%MZ}",
            loc="right")

    def save_image(self, path, fname, dpi=100, facecolor="white",
                   transparent=False, **kwargs):
        """
        Saves skewt as PNG image.

        Parameters
        ----------
        path : str
            Path of pirecotry relative to script.
        fname : str
            Filename of image.
        dpi : float or 'figure' (default=100)
            The resolution in dots per inch. If 'figure', use the figure's dpi
            value.
        facecolor : color or 'auto' (default="white")
            The facecolor of the figure. If 'auto', use the current figure
            facecolor.
        transparent : bool (default=False)
            Set background transparency.
        """
        Path(path).mkdir(parents=True, exist_ok=True)
        if transparent:
            self.fig.savefig(
                f"{path}/{fname}.png", dpi=dpi, transparent=True, **kwargs)
        else:
            self.fig.savefig(
                f"{path}/{fname}.png", dpi=dpi, facecolor=facecolor, **kwargs)
        plt.close()

def plot_sounding(data, lon, lat, station_name, fcst_time, init_time, ens):
    point = data["T2"].metpy.cartopy_crs.transform_point(
        lon, lat, ccrs.PlateCarree())
    data = data.metpy.assign_y_x(tolerance=10 * units("meter"))
    point_data = data.interp(
        {"west_east": point[0], "south_north": point[1]}, method="linear")
    prof = profile.create_profile(
        profile="convective", pres=wrf_calc.pressure(point_data) / 100,
        hght=wrf_calc.height(point_data),
        tmpc=thermo.ktoc(wrf_calc.temperature(point_data)),
        dwpc=wrf_calc.dewpoint(point_data),
        u=utils.MS2KTS(point_data["umet"]),
        v=utils.MS2KTS(point_data["vmet"]),
        omeg=point_data["omega"])
    t_adv, t_adv_bounds = params.inferred_temp_adv(prof, lat=lat)

    skewt = Skewt()

    skewt.plot_line_skewt(
        prof.pres, prof.wetbulb, color="cyan", linewidth=1)
    skewt.plot_line_skewt(
        prof.dpcl_ptrace, prof.dpcl_ttrace, color="purple", linewidth=1,
        linestyle="--")
    skewt.plot_line_skewt(
        prof.pres, prof.vtmp, color="red", linewidth=1,
        linestyle="--")
    skewt.plot_line_skewt(
        prof.pres, prof.tmpc, color="red", linewidth=2)
    skewt.plot_line_skewt(
        prof.pres, prof.dwpc, color="green", linewidth=2)

    if prof.mlpcl.bplus > 0:
        skewt.plot_line_skewt(
            prof.mlpcl.ptrace, prof.mlpcl.ttrace, color="orange", linewidth=1,
            linestyle="--")

    skewt.plot_heights(prof.pres, prof.hght - prof.hght.min())
    skewt.plot_surface_values(prof.pres.max(), prof.tmpc[0], prof.dwpc[0])
    if prof.mlpcl.bplus > 0:
        skewt.plot_sounding_levels(
            prof.mlpcl.lclpres, prof.mlpcl.lfcpres, prof.mlpcl.elpres)
    skewt.plot_max_lr(
        prof.max_lapse_rate_2_6[1], prof.max_lapse_rate_2_6[2],
        prof.max_lapse_rate_2_6[0])
    skewt.plot_inflow(
        prof.ebottom, prof.etop, prof.ebotm, prof.etopm, prof.esrh[0])

    skewt.plot_omega(prof.pres, prof.omeg)
    skewt.plot_dgz(prof.dgz_pbot, prof.dgz_ptop)

    skewt.plot_barbs(prof.pres, prof.u, prof.v)

    skewt.plot_temperature_advection(t_adv, t_adv_bounds)

    skewt.plot_hodo(prof.hght - prof.hght.min(), prof.u, prof.v)
    skewt.plot_storm_movers(
        prof.bunkers[0], prof.bunkers[1], prof.bunkers[2], prof.bunkers[3],
        *utils.vec2comp(prof.mean_lcl_el[0], prof.mean_lcl_el[1]))
    ebot_wind = interp.components(prof, prof.ebottom)
    etop_wind = interp.components(prof, prof.etop)
    skewt.plot_hodo_inflow(
        ebot_wind[0], ebot_wind[1], etop_wind[0], etop_wind[1], prof.bunkers[0],
        prof.bunkers[1])

    fontsize = 10
    thermo_x = [0.03, 0.07, 0.115, 0.155, 0.185, 0.22, 0.275, 0.325]
    thermo_y = np.linspace(0.23, 0.11, 5)
    thermo_vars = [
        "", "CAPE", "3CAPE", "CIN", "LI", "LCL(m)", "LFC(m)", "EL(m)"]
    for x, var in zip(thermo_x, thermo_vars):
        skewt.fig.text(
            x, thermo_y[0], var, weight="bold", fontsize=fontsize, ha="center",
            va="center")
    parcels = ["", "SFC", "ML", "MU", "FCST"]
    for y, pcl in zip(thermo_y, parcels):
        skewt.fig.text(
            thermo_x[0], y, pcl, weight="bold", fontsize=fontsize, ha="center",
            va="center")

    for i, parcel in enumerate([prof.sfcpcl, prof.mlpcl, prof.mupcl, prof.fcstpcl]):
        skewt.fig.text(
            thermo_x[1], thermo_y[i+1], utils.INT2STR(parcel.bplus),
            fontsize=fontsize, ha="center", va="center")
        skewt.fig.text(
            thermo_x[2], thermo_y[i+1], utils.INT2STR(parcel.b3km),
            fontsize=fontsize, ha="center", va="center")
        skewt.fig.text(
            thermo_x[3], thermo_y[i+1], utils.INT2STR(parcel.bminus),
            fontsize=fontsize, ha="center", va="center")
        skewt.fig.text(
            thermo_x[4], thermo_y[i+1], utils.INT2STR(parcel.li5),
            fontsize=fontsize, ha="center", va="center")
        skewt.fig.text(
            thermo_x[5], thermo_y[i+1], utils.INT2STR(parcel.lclhght),
            fontsize=fontsize, ha="center", va="center")
        skewt.fig.text(
            thermo_x[6], thermo_y[i+1], utils.INT2STR(parcel.lfchght),
            fontsize=fontsize, ha="center", va="center")
        skewt.fig.text(
            thermo_x[7], thermo_y[i+1], utils.INT2STR(parcel.elhght),
            fontsize=fontsize, ha="center", va="center")

    wind_x = [0.4, 0.46, 0.515, 0.58, 0.64]
    wind_y = thermo_y
    wind_vars = ["", "SRH", "Shear(kt)", "MW(kt)", "SRW(kt)"]
    for x, var in zip(wind_x, wind_vars):
        skewt.fig.text(
            x, wind_y[0], var, weight="bold", fontsize=fontsize, ha="center",
            va="center")
    layers = ["", "SFC-1km", "SFC-3km", "Eff Inflow", "SFC-6km"]
    for y, layer in zip(wind_y, layers):
        skewt.fig.text(
            wind_x[0], y, layer, weight="bold", fontsize=fontsize, ha="center",
            va="center")
    mean_eff = winds.mean_wind(prof, prof.ebottom, prof.etop)
    mean_eff = utils.comp2vec(prof.mean_eff[0], prof.mean_eff[1])
    srw_eff = utils.comp2vec(prof.srw_eff[0], prof.srw_eff[1])
    wind_values = [
        [utils.INT2STR(prof.srh1km[0]),
         utils.INT2STR(utils.mag(prof.sfc_1km_shear[0], prof.sfc_1km_shear[1])),
         f"{utils.INT2STR(prof.mean_1km[0])}/{utils.INT2STR(prof.mean_1km[1])}",
         f"{utils.INT2STR(prof.srw_1km[0])}/{utils.INT2STR(prof.srw_1km[1])}"],
        [utils.INT2STR(prof.srh3km[0]),
         utils.INT2STR(utils.mag(prof.sfc_3km_shear[0], prof.sfc_3km_shear[1])),
         f"{utils.INT2STR(prof.mean_3km[0])}/{utils.INT2STR(prof.mean_3km[1])}",
         f"{utils.INT2STR(prof.srw_3km[0])}/{utils.INT2STR(prof.srw_3km[1])}"],
        [utils.INT2STR(prof.esrh[0]),
         utils.INT2STR(utils.mag(prof.eff_shear[0], prof.eff_shear[1])),
         f"{utils.INT2STR(mean_eff[0])}/{utils.INT2STR(mean_eff[1])}",
         f"{utils.INT2STR(srw_eff[0])}/{utils.INT2STR(srw_eff[1])}"],
        ["",
         utils.INT2STR(utils.mag(prof.sfc_6km_shear[0], prof.sfc_6km_shear[1])),
         f"{utils.INT2STR(prof.mean_6km[0])}/{utils.INT2STR(prof.mean_6km[1])}",
         f"{utils.INT2STR(prof.srw_6km[0])}/{utils.INT2STR(prof.srw_6km[1])}"],
    ]
    for ix, x in enumerate(wind_x[1:]):
        for iy, y in enumerate(wind_y[1:]):
            skewt.fig.text(
                x, y, wind_values[iy][ix], fontsize=fontsize, ha="center",
                va="center")

    lambda_x = [0.725, 0.78]
    lambda_y = thermo_y
    skewt.fig.text(
        lambda_x[1], lambda_y[0], "Γ",weight="bold",
        fontsize=fontsize, ha="center", va="center")
    layers = ["", "SFC-3km", "3-6km", "850-500mb", "700-500mb"]
    lambda_values = [
        "", prof.lapserate_3km, prof.lapserate_3_6km, prof.lapserate_850_500,
        prof.lapserate_700_500]
    for y, layer, value in zip(wind_y, layers, lambda_values):
        skewt.fig.text(
            lambda_x[0], y, layer, weight="bold", fontsize=fontsize,
            ha="center", va="center")
        skewt.fig.text(
            lambda_x[1], y, utils.FLOAT2STR(value, 1), fontsize=fontsize,
            ha="center", va="center")

    lower_y = [0.07, 0.04]
    lower_x = [0.03, 0.08, 0.14, 0.19, 0.25, 0.3, 0.41, 0.5, 0.56, 0.61]
    skewt.fig.text(
        lower_x[0], lower_y[0], "DCAPE:", weight="bold", fontsize=fontsize,
        ha="center", va="center")
    skewt.fig.text(
        lower_x[0], lower_y[1], "DownT:", weight="bold", fontsize=fontsize,
        ha="center", va="center")
    skewt.fig.text(
        lower_x[1], lower_y[0], utils.INT2STR(prof.dcape), fontsize=fontsize,
        ha="center", va="center")
    skewt.fig.text(
        lower_x[1], lower_y[1],
        f"{utils.INT2STR(thermo.ctof(prof.dpcl_ttrace.max()))}°F",
        fontsize=fontsize, ha="center", va="center")

    skewt.fig.text(
        lower_x[2], lower_y[0], "ConvT:", weight="bold", fontsize=fontsize,
        ha="center", va="center")
    skewt.fig.text(
        lower_x[2], lower_y[1], "MaxT:", weight="bold", fontsize=fontsize,
        ha="center", va="center")
    skewt.fig.text(
        lower_x[3], lower_y[0], f"{utils.INT2STR(prof.convT)}°F",
        fontsize=fontsize, ha="center", va="center")
    skewt.fig.text(
        lower_x[3], lower_y[1], f"{utils.INT2STR(prof.maxT)}°F",
        fontsize=fontsize, ha="center", va="center")

    skewt.fig.text(
        lower_x[4], lower_y[0], "PWAT:", weight="bold", fontsize=fontsize,
        ha="center", va="center")
    skewt.fig.text(
        lower_x[4], lower_y[1], "Precip:", weight="bold", fontsize=fontsize,
        ha="center", va="center")
    skewt.fig.text(
        lower_x[5], lower_y[0], f"{utils.FLOAT2STR(prof.pwat, 2)}in",
        fontsize=fontsize, ha="center", va="center")
    skewt.fig.text(
        lower_x[5], lower_y[1], prof.precip_type, fontsize=fontsize,
        ha="center", va="center")

    skewt.fig.text(
        lower_x[6], lower_y[0], "150-350mb AGL RH:", weight="bold",
        fontsize=fontsize, ha="center", va="center")
    skewt.fig.text(
        lower_x[6], lower_y[1], "SFC-150mb AGL RH:", weight="bold",
        fontsize=fontsize, ha="center", va="center")
    skewt.fig.text(
        lower_x[7], lower_y[0], f"{utils.INT2STR(prof.mid_rh)}%",
        fontsize=fontsize, ha="center", va="center")
    skewt.fig.text(
        lower_x[7], lower_y[1], f"{utils.INT2STR(prof.low_rh)}%",
        fontsize=fontsize, ha="center", va="center")

    skewt.fig.text(
        lower_x[8], lower_y[0], "FZL:", weight="bold",
        fontsize=fontsize, ha="center", va="center")
    skewt.fig.text(
        lower_x[8], lower_y[1], "WBZ:", weight="bold",
        fontsize=fontsize, ha="center", va="center")
    skewt.fig.text(
        lower_x[9], lower_y[0],
        f"{utils.INT2STR(interp.hght(prof, params.temp_lvl(prof, 0)))}m",
        fontsize=fontsize, ha="center", va="center")
    skewt.fig.text(
        lower_x[9], lower_y[1],
        f"{utils.INT2STR(interp.hght(prof, params.temp_lvl(prof, 0, wetbulb=True)))}m",
        fontsize=fontsize, ha="center", va="center")

    haz_x = 0.71
    skewt.fig.text(
        haz_x, lower_y[0], "Hazzard:", weight="bold", fontsize=fontsize,
        ha="center", va="center")
    skewt.fig.text(
        haz_x, lower_y[1], prof.watch_type, weight="bold", fontsize=fontsize,
        ha="center", va="center")

    skewt.plot_thetae(prof.pres, prof.thetae)

    skewt.plot_title(station_name, lon, lat, fcst_time, init_time)

    skewt.save_image(
        f"../images_wrf/{init_time:%Y%m%d%H}/{ens}/skewt/{station_name}/",
        f"skewt_{station_name}_{fcst_time:%Y%m%d%H%M}")

def plot_hour(init_time, fcst_time, domain, ens):
    print("Reading data for", fcst_time, "domain", domain, flush=True)
    ds_wrf = read_data_wrf_sounding(init_time, fcst_time, domain, ens)
    ds_wrf.ds = add_cf_wrf(ds_wrf.ds)
    for station in wrf_sounding_stations[domain].values():
        plot_sounding(
            data=ds_wrf.ds,
            lon=station.lon,
            lat=station.lat,
            station_name=station.name,
            fcst_time=fcst_time,
            init_time=start_date,
            ens=ens)

def main():
    with mp.Pool(processes=1, maxtasksperchild=1) as pool:
        for domain in np.arange(
                config.wrf_max_domain, config.wrf_min_domain - 1, -1):
            if wrf_sounding_stations[domain] != {}:
                for time in pd.date_range(
                        start_date, end_date, freq=config.wrf_freq[domain]):
                    print(time, flush=True)
                    pool.apply_async(
                        plot_hour,
                        args=(start_date, time, domain, args.ens, ),
                        error_callback=error_handler)
        pool.close()
        pool.join()


if __name__ == "__main__":
    mp.set_start_method("fork")
    main()
