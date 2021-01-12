# climate_modelling

"""
Created: 9 January 2021
Last modified: 12 January 2021

This module contains functions that are likely to be used repeatedly in computer practicals in the course Climate Modelling at VU Amsterdam, Jan 2021.

@author: AJSchiller, heavily based on code by Didier M. Roche a.k.a. dmr.
"""


import numpy as np
from numpy import ma
import netCDF4 as nc
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import glob


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
    """

    var = data.variables[variable]
    var_plot =  ma.masked_equal(var, var.missing_value)

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
    """

    min_bounds = np.min(var_plot)
    max_bounds = np.max(var_plot)
    nbs_bounds = 15
    fix_bounds = np.linspace(min_bounds, max_bounds, nbs_bounds)

    the_chosen_map = plt.cm.coolwarm

    fig = plt.figure(figsize=(10,10))

    ax = plt.axes(projection=ccrs.NearsidePerspective(central_latitude=central_latitude, central_longitude=central_longitude, satellite_height=satellite_height))

    ax.set_title(title)
    ax.set_global()
    fig.set_facecolor("grey")

    mesh = ax.pcolormesh(lons, lats, var_plot, cmap=the_chosen_map,
                         transform=ccrs.PlateCarree(), edgecolor="white",linewidths=0.05,
                         vmin=min_bounds, vmax=max_bounds)

    if contours:
        CS = ax.contour(lons, lats, var_plot, fix_bounds, transform=ccrs.PlateCarree(), colors='k', linewidths=1.0, vmin=min_bounds, vmax=max_bounds)
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
    """
    crs = ccrs.PlateCarree()

    min_bounds = np.min(var_map)
    max_bounds = np.max(var_map)
    nbs_bounds = 15
    fix_bounds = np.linspace(min_bounds, max_bounds, nbs_bounds)

    fig, ax = plt.subplots(figsize=(10,5), subplot_kw=dict(projection=crs))
    ax.set_title(title)
    ax.set_global()

    mesh = ax.pcolormesh(lons, lats, var_map, cmap=cmap,
                         transform=crs,
                         vmin=min_bounds, vmax=max_bounds)

    if contours:
        CS = ax.contour(lons, lats, var_map, fix_bounds, transform=crs, colors='k', linewidths=1.0, vmin=min_bounds, vmax=max_bounds)
        ax.clabel(CS, fontsize=9, inline=1)

    plt.colorbar(mesh, orientation='vertical', shrink=0.61, label=legend_label)
    ax.gridlines()
    ax.coastlines()
    fig.show()

    if out_file:
        fig.savefig(out_file)


############################################################

def smooth(x, window_len, window):
    """Smooths the input data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal (with the window size) at both ends, so that transient parts are minimized at the begining and end parts of the output signal.

    Keyword arguments:
    x -- the input signal (array).
    window_len -- the dimension of the smoothing window; should be an odd integer.
    window -- the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
                flat window will produce a moving average smoothing. See numpy documentation for explanation of these window types.
    """

    if x.ndim != 1:
        raise(ValueError, "smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise(ValueError, "Input vector needs to be bigger than window size.")


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise(ValueError, "Window must be one of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")


    s=np.r_[2*x[0]-x[window_len:1:-1],x,2*x[-1]-x[-1:-window_len:-1]]
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='same')
    y2=y[window_len-1:-window_len+1]
    return y2


############################################################

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


############################################################

def plot_lecture_data(lecture_data, out_file=False):
    """Graphs a scatterplot of the output of `lecture_data` with a 10-year running average.

    Keyword arguments:
    lecture_data -- output of the `lecture_data` function.
    out_file -- optional, filepath to which to save the output graph.
    """

    sizef=3
    orange='#FF7F00'
    yellow='#EEC900'
    sizefont=13

    fig=plt.figure(1, figsize=(12,7))

    plt.subplot(1, 1, 1)

    results=lecture_data
    temp=results[0]

    temp_all=temp
    var = np.mean(temp)
    nb_max = len(temp)
    my_window_len = 10
    time_all=np.arange(0,0+len(temp_all))
    temp_all_smooth=smooth(temp_all, window_len=my_window_len, window='hanning')
    temp_all_smooth=temp_all_smooth[my_window_len-1:-my_window_len+1]
    time_all_smooth=time_all[my_window_len-1:-my_window_len+1]
    print('size',len(temp_all_smooth),len(time_all_smooth))

    compteur=0
    line2=plt.plot(time_all,temp,'.',color='black', markersize=1,label="Yearly data")
    line2=plt.plot(time_all_smooth,temp_all_smooth,'-',color='indianred', linewidth=2,label=""+str(my_window_len)+" years running average")


    plt.ylabel('Global mean temperature (°C)',{'fontsize': sizefont})
    plt.xlabel('Time (simulation years)',{'fontsize': sizefont})
    plt.legend(fontsize='large')
    plt.subplots_adjust(bottom=0.1, right=0.90, top=0.92, left=0.15, hspace=0.2)

    if out_file:
        fig.savefig(out_file)

    plt.show()


############################################################
