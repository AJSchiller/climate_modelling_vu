# climate_modelling

"""
Created: 9 January 2020

This module contains functions that are likely to be used repeatedly in computer practicals in the course Climate Modelling at VU Amsterdam, Jan 2021.

@author: AJSchiller
"""


import numpy as np
from numpy import ma
import netCDF4
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import glob


def read_temp(filename):
    """Reads a NetCDF file and extracts temperature, latitude, and longitude data.

    Keyword arguments:
    filename -- a string containing the path to the file to be read.

    The temperature variable in the file should have the name "ts".
    Returns temperature, latitude, and longitude values stored in masked arrays.
    """

    f = netCDF4.Dataset(filename)
    var = f.variables["ts"]
    var_plot =  ma.masked_equal(var,var.missing_value)-273.15

    lons = f.variables['lon'][...]
    lats = f.variables['lat'][...]

    return var_plot, lons, lats


def temp_map(var_plot, lons, lats, central_latitude=50.72, central_longitude=-3.53, satellite_height=10000000.0, out_file=False):
    """Produces a map of the output data from the `read_temp` function.

    Keyword arguments:
    var_plot -- a masked array of temperature data, as returned by the function `read_temp`.
    lons -- a masked array of longitude values, as returned by the function `read_temp`.
    lats -- a masked array of latitude values, as returned by the function `read_temp`.
    central_latitude -- the latitude at which the output map will be centre (default 50.72).
    central_longitude -- the longitude at which the output map will be centre (default -3.53).
    satellite_height -- the height from which the map is viewed (default 10000000.0).
    out_file -- optional filepath to save the output image (default False).
    """

    min_bounds = np.min(var_plot)
    max_bounds = np.max(var_plot)
    nbs_bounds = 15
    fix_bounds = np.linspace(min_bounds,max_bounds,nbs_bounds)

    the_chosen_map = plt.cm.coolwarm

    fig = plt.figure(figsize=(10,10))

    ax = plt.axes(projection=ccrs.NearsidePerspective(central_latitude=central_latitude,
                                                    central_longitude=central_longitude,
                                                    satellite_height=satellite_height))

    ax.set_global()
    fig.set_facecolor("grey")

    mesh = ax.pcolormesh(lons, lats, var_plot, cmap=the_chosen_map,
                         transform=ccrs.PlateCarree(), edgecolor="white",linewidths=0.05,
                         vmin=min_bounds, vmax=max_bounds)

    CS = ax.contour(lons, lats, var_plot, fix_bounds, transform=ccrs.PlateCarree(), colors='k', linewidths=1.0, vmin=min_bounds, vmax=max_bounds)
    ax.clabel(CS, fontsize=9, inline=1)

    plt.colorbar(mesh, orientation='horizontal', shrink=0.75)
    ax.gridlines()
    ax.coastlines()
    plt.show()

    if out_file:
        fig.savefig(out_file)


def smooth(x, window_len, window):
    """Smooths the input data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal (with the window size) at both ends, so that transient parts are minimized at the begining and end parts of the output signal.

    Keyword arguments:
    x -- the input signal (array).
    window_len -- the dimension of the smoothing window; should be an odd integer.
    window -- the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
                flat window will produce a moving average smoothing.

    Adapted from Didier M. Roche a.k.a. dmr.
    """

    if x.ndim != 1:
        raise(ValueError, "smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise(ValueError, "Input vector needs to be bigger than window size.")


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise(ValueError, "Window must be one of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")


    s=numpy.r_[2*x[0]-x[window_len:1:-1],x,2*x[-1]-x[-1:-window_len:-1]]
    if window == 'flat': #moving average
        w=numpy.ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')

    y=numpy.convolve(w/w.sum(),s,mode='same')
    y2=y[window_len-1:-window_len+1]
    return y2


def lecture_data(filename):
    """Reads data from input into an array.

    Keyword arguments:
    filename -- path to input file.
    """

    input_f = filename
    f = open(input_f)
    f.seek(0)
    nblines = f.readlines()
    compteur_line = 0
    compteur_line2 = 0
    dada = 0
    nmax = len(nblines) - 1
    temp = np.zeros(nmax)
    nmax2 = int(nmax/360)
    temp2 = np.zeros(nmax2)
    for i in range(0, len(nblines)-1):
        line = nblines[i]
        temp[compteur_line] = float(line.strip().split()[2])
        compteur_line += 1
        if dada == 0:
            toto = float(line.strip().split()[2])
        if dada < 359:
            toto += float(line.strip().split()[2])
        if dada == 359:
            toto += float(line.strip().split()[2])
            toto = toto/360.
            temp2[compteur_line2] = toto
            compteur_line2 += 1
        dada += 1
        if dada == 360:
            dada = 0
    results = [temp2]
    return results
