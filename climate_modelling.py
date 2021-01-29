# climate_modelling

"""
Created: 9 January 2021
Last modified: 28 January 2021

This module contains functions that are likely to be used repeatedly in computer practicals in the course Climate Modelling at VU Amsterdam, Jan 2021.

@author: AJSchiller
"""


import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import pandas as pd
from cartopy.util import add_cyclic_point
from matplotlib.colors import TwoSlopeNorm


############################################################

def list_vars(data):
    """Prints each variable in the dataset and its long-form name.

    Keyword arguments:
    data -- a NetCDF dataset.
    """

    for var_name in data.variables.keys():
        variable = data.variables[var_name]
        if 'long_name' in variable.ncattrs():
            long_name = getattr(variable, 'long_name')
        else:
            long_name = "N/A"
        print('{} = {}'.format(var_name, long_name))


############################################################

def get_var(data, variable):
    """Extracts the variable of interest from a NetCDF dataset by latitude and longitude.

    Keyword arguments:
    data -- a NetCDF dataset.
    variable -- name of the variable to be extracted (other than latitude and longitude).

    Assumes the variable of interest has either 2 or 3 dimensions (2 of them being 'lat' and 'lon'). If 3 dimensions are present, the function will average the data along the 3rd dimension (i.e. not 'lat' or 'lon')
    Returns variable of interest, latitude, and longitude values stored in masked arrays.

    Adapted from Didier M. Roche a.k.a. dmr.
    """

    var = data.variables[variable]
    var_plot =  np.ma.masked_equal(var, var.missing_value)

    lons = data.variables['lon'][...]
    lats = data.variables['lat'][...]

    return var_plot, lons, lats


############################################################

def plot_nearside(var_plot, lons, lats, title, legend_label, cmap='viridis', contours=True, central_latitude=50.72, central_longitude=-3.53, satellite_height=10000000.0, out_file=False):
    """Produces a map of the output data from the `get_var` function. Uses the `NearsidePerspective` projection.

    Keyword arguments:
    var_plot -- a masked array of the data. This must be 2-dimensions (lat, lon).
    lons -- a masked array of longitude values.
    lats -- a masked array of latitude values.
    title -- figure title for the map.
    legen_label -- label for the map's colorbar.
    cmap -- colormap to be used (default 'viridis', suggested alternative plt.cm.coolwarm).
    contours -- boolean, determines whether contour lines will be included in the map (default True).
    central_latitude -- the latitude at which the output map will be centre (default 50.72).
    central_longitude -- the longitude at which the output map will be centre (default -3.53).
    satellite_height -- the height from which the map is viewed (default 10000000.0).
    out_file -- optional filepath to save the output image (default False).

    Adapted from Didier M. Roche a.k.a. dmr.
    """

    crs = ccrs.PlateCarree()

    min_bounds = np.min(var_plot)
    max_bounds = np.max(var_plot)
    nbs_bounds = 15
    fix_bounds = np.linspace(min_bounds, max_bounds, nbs_bounds)

    var_plot, lons = add_cyclic_point(var_plot, coord=lons)

    fig = plt.figure(figsize=(10,10))

    ax = plt.axes(projection=ccrs.NearsidePerspective(central_latitude=central_latitude, central_longitude=central_longitude, satellite_height=satellite_height))

    ax.set_title(title)
    ax.set_global()
    fig.set_facecolor("grey")

    mesh = ax.pcolormesh(lons, lats, var_plot, cmap=cmap,
                         transform=crs, edgecolor="white",linewidths=0.05,
                         vmin=min_bounds, vmax=max_bounds)

    if contours:
        CS = ax.contour(lons, lats, var_plot, fix_bounds, transform=crs, colors='k', linewidths=1.0, vmin=min_bounds, vmax=max_bounds)
        ax.clabel(CS, fontsize=9, inline=1)

    plt.colorbar(mesh, orientation='horizontal', shrink=0.75, label=legend_label)
    ax.gridlines()
    ax.coastlines()
    plt.show()

    if out_file:
        fig.savefig(out_file)


############################################################

def map_data(var_map, lons, lats, title, legend_label, cmap='viridis', contours=True, out_file=False):
    """Generates a map of output data from `get_var` function.

    Keyword arguments:
    var_plot -- a masked array of the data. This must be 2-dimensions (lat, lon).
    lons -- a masked array of longitude values.
    lats -- a masked array of latitude values.
    title -- figure title for the map.
    legen_label -- label for the map's colorbar.
    cmap -- colormap to be used (default 'viridis', suggested alternative plt.cm.coolwarm).
    contours -- boolean, determines whether contour lines will be included in the map (default True).
    out_file -- optional filepath to save the output image (default False).

    Adapted from Didier M. Roche a.k.a. dmr.
    """
    crs = ccrs.PlateCarree()

    min_bounds = np.min(var_map)
    max_bounds = np.max(var_map)
    nbs_bounds = 15
    fix_bounds = np.linspace(min_bounds, max_bounds, nbs_bounds)

    var_map, lons = add_cyclic_point(var_map, coord=lons)

    fig, ax = plt.subplots(figsize=(10,5), subplot_kw=dict(projection=crs))
    ax.set_title(title)
    ax.set_global()

    mesh = ax.pcolormesh(lons, lats, var_map, norm=TwoSlopeNorm(0), cmap=cmap,
                         transform=crs,
                         vmin=min_bounds, vmax=max_bounds)

    if contours:
        CS = ax.contour(lons, lats, var_map, fix_bounds, transform=crs, colors='k', linewidths=1.0, vmin=min_bounds, vmax=max_bounds)
        ax.clabel(CS, fontsize=9, inline=1)

    plt.colorbar(mesh, orientation='vertical', shrink=0.61, label=legend_label)
    ax.set_xticks(np.linspace(-180, 180, num=7), crs=crs)
    ax.set_yticks(np.linspace(-90, 90, num=7), crs=crs)
    ax.gridlines()
    ax.coastlines()
    fig.show()

    if out_file:
        fig.savefig(out_file)


############################################################

def read_book(filename):
    """Reads data from filename into a pandas dataframe.
    The input data is assumed to be daily.
    Returns a pandas series of yearly averages.

    Keyword arguments:
    filename -- path to input file.
    """

    df = pd.read_csv(filename, sep='\s+', header=None, names=['year', 'day', 'var'])

    results = df.groupby('year')['var'].mean()

    return results


############################################################

def plot_timeline(data, title, ylabel, xlabel='Time (model years)', window_length=10, window='hann', legend='on', out_file=False):
    """Graphs a scatterplot of the data (assumed to be annual) with a 10-year running average.

    Keyword arguments:
    data -- annual data (suggested: output of the `read_book` function).
    window -- string, the type of window from the list: https://docs.scipy.org/doc/scipy/reference/signal.windows.html (default 'hann').
    ylabel -- label for the vertical (y) axis.
    xlabel -- label for the horizontal (x) axis (default 'Time (model years)').
    window_length -- the dimension of the smoothing window (default 10).
    legend -- includes a legend with the default labels "Yearly data" and "x-year running average" (default 'on'). To disable the legend, set to 'off'.
    out_file -- optional, filepath to which to save the output graph (default False).
    """

    df = pd.DataFrame(data)
    var = df.columns[0]
    df['rolling'] = df[var].rolling(10, win_type='hann', center=True).mean()

    fig, ax = plt.subplots()
    ax.set_title(title)
    ax.plot(df.index, df[var], '.', color='black', markersize=1, label="Yearly data")
    ax.plot(df.index[10:], df['rolling'][10:], '-', color='indianred', linewidth=2, label="{}-year running average".format(window_length))

    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)

    if legend == 'on':
        ax.legend()

    if out_file:
        fig.savefig(out_file)

    fig.show()


############################################################
