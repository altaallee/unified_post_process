import cartopy.crs as ccrs
import matplotlib.pyplot as plt
plt.switch_backend("Agg")
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patheffects as path_effects
from pathlib import Path
import cartopy.feature as cfeat
import metpy
from metpy.units import units
from scipy import ndimage
import numpy as np


class SingleMap:
    def __init__(self, figsize=[12, 8], projection=ccrs.Mercator()):
        """
        Creates an image with a single map.

        Parameters
        ----------
        figsize : (float, float) (default=[12, 8])
            Width, height in inches.
        projection : cartopy.crs (default=cartopy.crs.Mercator())
            The projection type of the subplot (Axes).
        """
        self.fig = plt.figure(figsize=figsize)
        self.ax = self.fig.add_subplot(111, projection=projection)
        self.ax.set_adjustable("datalim")
        plt.subplots_adjust(bottom=0.07, left=0.01, top=0.92, right=0.99)
        self.has_colorbar = False

    def draw_barbs(self, x, y, u, v, **kwargs):
        """
        Draws wind barbs on map.

        Parameters
        ----------
        x, y : 1D or 2D array-like
            The x and y coordinates of the barb locations. See pivot for how the
            barbs are drawn to the x, y positions. If x and y are 1D but u, v
            are 2D, x, y are expanded to 2D using x, y = np.meshgrid(x, y). In
            this case len(x) and len(y) must match the column and row dimensions
            of u and v.
        u, v : 1D or 2D array-like
            The x and y components of the barb shaft.
        """
        self.ax.barbs(x, y, u, v, transform=ccrs.PlateCarree(), **kwargs)

    def draw_contour(self, x, y, z, label=True, format="%1.0f", **kwargs):
        """"
        Draws contours on map.

        Parameters
        ----------
        x, y : array-like
            The coordinates of the values in z. x and y must both be 2D with the
            same shape as z (e.g. created via numpy.meshgrid), or they must both
            be 1-D such that len(x) == N is the number of columns in z and
            len(Y) == M is the number of rows in z. x and y must both be ordered
            monotonically.
        z : (M, N) array-like
            The height values over which the contour is drawn.
        label : boolean (default=True)
            Label contour lines.
        format : str (default=f"%1.0f")
            Format of contour labels.
        """
        cs = self.ax.contour(
            x, y, z, transform=ccrs.PlateCarree(), transform_first=True,
            **kwargs)
        if label:
            self.ax.clabel(cs, cs.levels, fmt=format)

    def draw_contourf(self, x, y, z, plot_cbar=True, labelsize=12,
                      cmap_label="", ticks=None, **kwargs):
        """
        Draws filled contours on map.

        Parameters
        ----------
        x, y : array-like
            The coordinates of the values in z. x and y must both be 2D with the
            same shape as z (e.g. created via numpy.meshgrid), or they must both
            be 1-D such that len(x) == N is the number of columns in z and
            len(Y) == M is the number of rows in z. x and y must both be ordered
            monotonically.
        z : (M, N) array-like
            The height values over which the contour is drawn.
        plot_cbar : boolean (default=True)
            Choose to plot colorbar.
        labelsize : int (default=12)
            Size of colorbar label.
        cmap_label : str (default=None)
            The label on the colorbar's long axis.
        ticks : None or list of ticks or Locator
            If None, ticks are determined automatically from the input.
        """
        cf = self.ax.contourf(
            x, y, z, transform=ccrs.PlateCarree(), transform_first=True,
            **kwargs)
        if (not self.has_colorbar) and plot_cbar:
            self.has_colorbar = True
            divider = make_axes_locatable(self.ax)
            cax = divider.append_axes(
                "bottom", size="2%", pad=0, axes_class=plt.Axes)
            cb = self.fig.colorbar(
                cf, cax=cax, orientation="horizontal", extendfrac=0,
                ticks=ticks)
            cb.ax.tick_params(labelsize=labelsize)
            cb.set_label(label=cmap_label, size=labelsize)

    def draw_hilo_mpas(self, data, gaussian_sigma, size):
        """
        Draws L and H at low and high points for MPAS.

        Parameters
        ----------
        data : xarray.DataArray
            Data for L and H markers.
        gaussian_sigma : float
            Standard deviation for Gaussian kernel.
        size : int
            Minimum distance between each L and H marker.
        """
        data_smooth = ndimage.gaussian_filter(data, gaussian_sigma)
        data_min = ndimage.minimum_filter(data_smooth, size, mode="nearest")
        miny, minx = np.where(data_smooth == data_min)
        for x, y in zip(minx, miny):
            data_point = data.isel({"longitude": x, "latitude": y})
            self.draw_text(
                data_point["longitude"], data_point["latitude"], "L",
                color="red", ha="center", va="center", fontsize=16,
                weight="bold", clip_on=True, clip_box=self.ax.bbox)
            self.draw_text(
                data_point["longitude"], data_point["latitude"],
                f"\n{int(np.round(data_point.values))}", color="red",
                ha="center", va="top", weight="bold", clip_on=True,
                clip_box=self.ax.bbox)

        data_max = ndimage.maximum_filter(data_smooth, size, mode="nearest")
        maxy, maxx = np.where(data_smooth == data_max)
        for x, y in zip(maxx, maxy):
            data_point = data.isel({"longitude": x, "latitude": y})
            self.draw_text(
                data_point["longitude"], data_point["latitude"], "H",
                color="blue", ha="center", va="center", fontsize=16,
                weight="bold", clip_on=True, clip_box=self.ax.bbox)
            self.draw_text(
                data_point["longitude"], data_point["latitude"],
                f"\n{int(np.round(data_point.values))}", color="blue",
                ha="center", va="top", weight="bold", clip_on=True,
                clip_box=self.ax.bbox)

    def draw_hilo_wrf(self, data, gaussian_sigma, size):
        """
        Draws L and H at low and high points for WRF.

        Parameters
        ----------
        data : xarray.DataArray
            Data for L and H markers.
        gaussian_sigma : float
            Standard deviation for Gaussian kernel.
        size : int
            Minimum distance between each L and H marker.
        """
        data_smooth = ndimage.gaussian_filter(data, gaussian_sigma)
        data_min = ndimage.minimum_filter(data_smooth, size, mode="nearest")
        miny, minx = np.where(data_smooth == data_min)
        for x, y in zip(minx, miny):
            data_point = data.isel({"west_east": x, "south_north": y})
            self.draw_text(
                data_point["XLONG"], data_point["XLAT"], "L",
                color="red", ha="center", va="center", fontsize=16,
                weight="bold", clip_on=True, clip_box=self.ax.bbox)
            self.draw_text(
                data_point["XLONG"], data_point["XLAT"],
                f"\n{int(np.round(data_point.values))}", color="red",
                ha="center", va="top", weight="bold", clip_on=True,
                clip_box=self.ax.bbox)

        data_max = ndimage.maximum_filter(data_smooth, size, mode="nearest")
        maxy, maxx = np.where(data_smooth == data_max)
        for x, y in zip(maxx, maxy):
            data_point = data.isel({"west_east": x, "south_north": y})
            self.draw_text(
                data_point["XLONG"], data_point["XLAT"], "H",
                color="blue", ha="center", va="center", fontsize=16,
                weight="bold", clip_on=True, clip_box=self.ax.bbox)
            self.draw_text(
                data_point["XLONG"], data_point["XLAT"],
                f"\n{int(np.round(data_point.values))}", color="blue",
                ha="center", va="top", weight="bold", clip_on=True,
                clip_box=self.ax.bbox)

    def draw_pcolormesh(self, x, y, c, plot_cbar=True, **kwargs):
        """
        Create a pseudocolor plot with a non-regular rectangular grid. X and Y
        can be used to specify the corners of the quadrilaterals.

        Parameters:
        -----------
        x, y : array-like, optional
            The coordinates of the corners of quadrilaterals of a pcolormesh:
                (X[i+1, j], Y[i+1, j])       (X[i+1, j+1], Y[i+1, j+1])
                                      +-----+
                                      |     |
                                      +-----+
                    (X[i, j], Y[i, j])       (X[i, j+1], Y[i, j+1])
            Note that the column index corresponds to the x-coordinate, and the
            row index corresponds to y. For details, see the Notes section
            below. If shading='flat' the dimensions of X and Y should be one
            greater than those of C, and the quadrilateral is colored due to the
            value at C[i, j]. If X, Y and C have equal dimensions, a warning
            will be raised and the last row and column of C will be ignored. If
            shading='nearest' or 'gouraud', the dimensions of X and Y should be
            the same as those of C (if not, a ValueError will be raised). For
            'nearest' the color C[i, j] is centered on (X[i, j], Y[i, j]). For
            'gouraud', a smooth interpolation is caried out between the
            quadrilateral corners. If X and/or Y are 1-D arrays or column
            vectors they will be expanded as needed into the appropriate 2D
            arrays, making a rectangular grid.
        c : 2D array-like
            The color-mapped values.
        plot_cbar : boolean (default=True)
            Choose to plot colorbar.
        """
        cm = self.ax.pcolormesh(x, y, c, transform=ccrs.PlateCarree(), **kwargs)
        if (not self.has_colorbar) and plot_cbar:
            self.has_colorbar = True
            divider = make_axes_locatable(self.ax)
            cax = divider.append_axes(
                "right", size="3%", pad=0, axes_class=plt.Axes)
            cb = self.fig.colorbar(
                cm, cax=cax, orientation="vertical", **kwargs)
            cb.ax.tick_params(labelsize=labelsize)
            cb.set_label(label=label, size=labelsize)

    def draw_quiver(self, x, y, u, v, **kwargs):
        """"
        Draws quivers on map.

        Parameters
        ----------
        x, y: 1D or 2D array-like
            The x and y coordinates of the arrow locations. If x and y are 1D
            but u, v are 2D, x, y are expanded to 2D using
            x, y = np.meshgrid(x, y). In this case len(x) and len(y) must match
            the column and row dimensions of u and v.
        U, V: 1D or 2D array-like
            The x and y direction components of the arrow vectors. They must
            have the same number of elements, matching the number of arrow
            locations. u and v may be masked. Only locations unmasked in u, v,
            and C will be drawn.
        """
        self.ax.quiver(x, y, u, v, transform=ccrs.PlateCarree(), **kwargs)

    def draw_shapefiles(self, hires=False, **kwargs):
        """
        Draws shapefiles onto map.

        Parameters
        ----------
        hires : cartopy.feature.ShapelyFeature (default=False)
            Features to plot other than cartopy.
        """
        if hires:
            self.ax.add_feature(hires, facecolor="none", **kwargs)
        else:
            self.ax.add_feature(cfeat.COASTLINE, **kwargs)
            self.ax.add_feature(cfeat.BORDERS, **kwargs)

    def draw_station_values_mpas(self, data, stations, priority, decimal=0,
                                 station_condition=False):
        """
        Draws values of stations onto map for MPAS.

        Parameters
        ----------
        data : xarray.DataArray
            Data to plot values from.
        stations : pandas.DataFrame
            Dataframe with longitude, latitude, and priority values.
        priority : int
            Value of maximum priority to plot.
        decimal : int (default=0)
            Number of decimal places to print.
        station_condition : function (default=False)
            Function returning bool to plot station value or not.
        """
        stations = stations.loc[stations["PRIORITY"] <= priority]
        out = data.metpy.cartopy_crs.transform_points(
            ccrs.PlateCarree(), stations["LON"], stations["LAT"])
        for proj_lon, proj_lat in zip(out[:, 0], out[:, 1]):
            point_data = data.interp(
                {"longitude": proj_lon, "latitude": proj_lat},
                method="linear")
            if not point_data.isnull() and (
                    (station_condition == False) or
                    station_condition(point_data.values)):
                txt = self.draw_text(
                    point_data["longitude"], point_data["latitude"],
                    f"{np.round(point_data.values, decimal):.{decimal}f}",
                    fontweight="bold", fontsize=12, ha="center", va="center",
                    clip_on=True, clip_box=self.ax.bbox)
                txt.set_path_effects(
                    [path_effects.Stroke(linewidth=2, foreground="white"),
                        path_effects.Normal()])

    def draw_station_values_wrf(self, data, stations, priority, decimal=0,
                                station_condition=False):
        """
        Draws values of stations onto map.

        Parameters
        ----------
        data : xarray.DataArray
            Data to plot values from.
        stations : pandas.DataFrame
            Dataframe with longitude, latitude, and priority values.
        priority : int
            Value of maximum priority to plot.
        decimal : int (default=0)
            Number of decimal places to print.
        station_condition : function (default=False)
            Function returning bool to plot station value or not.
        """
        stations = stations.loc[stations["PRIORITY"] <= priority]
        out = data.metpy.cartopy_crs.transform_points(
            ccrs.PlateCarree(), stations["LON"], stations["LAT"])
        data = data.metpy.assign_y_x(tolerance=10 * units("meter"))
        for proj_lon, proj_lat in zip(out[:, 0], out[:, 1]):
            point_data = data.interp(
                {"west_east": proj_lon, "south_north": proj_lat},
                method="linear")
            if not point_data.isnull() and (
                    (station_condition == False) or
                    station_condition(point_data.values)):
                txt = self.draw_text(
                    point_data["XLONG"], point_data["XLAT"],
                    f"{np.round(point_data.values, decimal):.{decimal}f}",
                    fontweight="bold", fontsize=12, ha="center", va="center",
                    clip_on=True, clip_box=self.ax.bbox)
                txt.set_path_effects(
                    [path_effects.Stroke(linewidth=2, foreground="white"),
                        path_effects.Normal()])

    def draw_text(self, x, y, s, **kwargs):
        """
        Add text to the Axes.

        Parameters
        ----------
        x, y : float
            The position to place the text. By default, this is in data
            coordinates.
        s : str
            The text.

        Returns
        -------
        Text
            The created Text instance.
        """
        return self.ax.text(x, y, s, transform=ccrs.PlateCarree(), **kwargs)

    def draw_title(self, label, loc="left", fontsize=16, **kwargs):
        """
        Draws title for map.

        Parameters
        ----------
        label : str
            Text to use for the title.
        loc : {'center', 'left', 'right'} (default="left")
            Which title to set.
        fontsize : int (default=16)
            Size of title.
        """
        self.ax.set_title(label, loc=loc, fontsize=fontsize, **kwargs)

    def save_image(self, path, fname, dpi=100, facecolor="white",
                   transparent=False, **kwargs):
        """
        Saves map as PNG image.

        Parameters
        ----------
        path : str
            Path of directory relative to script.
        fname : str
            Filename of image.
        dpi : float or 'figure' (default=100)
            The resolution in dots per inch. If 'figure', use the figure's dpi
            value.
        facecolor : color or 'auto' (default="white")
            The facecolor of the figure. If 'auto', use the current figure
            facecolor.
        transparent : bool (default=False)
            Set backgroud transparency.
        """
        Path(path).mkdir(parents=True, exist_ok=True)
        if transparent:
            self.fig.savefig(
                f"{path}/{fname}.png", dpi=dpi, transparent=True, **kwargs)
        else:
            self.fig.savefig(
                f"{path}/{fname}.png", dpi=dpi, facecolor=facecolor, **kwargs)
        plt.close()

    def set_extent(self, extent):
        """
        Sets extent of map.

        Parameters
        ----------
        extent : 4-tuple of float
            The position and size of the image as tuple (left, right, bottom,
            top) in data coordinates.
        """
        self.ax.set_extent(extent, crs=ccrs.PlateCarree())
