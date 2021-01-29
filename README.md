# Climate Modelling at VU Amsterdam

[![](https://img.shields.io/badge/python-3.7-blue.svg)](https://www.python.org/downloads/)

## Description

This module is intended for use in the graduate course Climate Modelling at VU Amsterdam. It was created based on the practical work assigned during the January 2021 fully-online presentation of the course. Mileage in future presentations may vary, and relevant contributions are of course welcome.

## Getting started

### Dependencies

This module requires that the following packages also be installed:

* [numpy](https://numpy.org/install/)
* [pandas](https://pandas.pydata.org/docs/getting_started/install.html)
* [matplotlib](https://matplotlib.org/users/installing.html)
* [netCDF4](https://pypi.org/project/netCDF4/)
* [cartopy](https://scitools.org.uk/cartopy/docs/latest/installing.html#installing)

If you are using Anaconda, you will already have some of these packages installed. Use Anaconda to install any missing packages.<br>
Otherwise, all packages are also available from PyPi. If you are unsure how to install a specific package, follow the respective link from the list above.

### Installation

To install this module, you can either click the green 'Code' button on the right of the [GitHub page](https://github.com/AJSchiller/climate_modelling_vu) and then click 'Download ZIP':
<p align="right">
<img style="float: right;" src="https://docs.github.com/assets/images/help/repository/code-button.png" alt="The 'Code' button (source: GitHub)" width="250"/>
</p>

<br><br>Or run the following command from your Terminal (linux/MacOS) or Command Prompt (Windows):

`git clone https://github.com/AJSchiller/climate_modelling_vu.git`

Before you can import this package into your python code, you will have to do one of the following:

1. Place the file `climate_modelling.py` in the same folder as the rest of your python packages;
2. Add the folder containing `climate_modelling.py` to your system's `PATH` variable (search the internet for instructions for your specific operating system); or
3. Place the file `climate_modelling.py` in the same folder as the python script to which you want to import it (this is undesireable, as you would need to do this for every script in which you want to import the module).

You should now be able to import the `climate_modelling` module into your python script using:

``` python
import climate_modelling as cm
```

If, for some reason, none of these options works, you can add the following to the start of your script (again, you would need to do this in every script that uses the module):

``` python
import sys
sys.path.append(`path/to/the/folder/containing/the/module`)
import climate_modelling as cm
```

## How to use this module

This module is intended to help students quickly visualise selected model outputs for the practicals in the course Climate Modelling, allowing you to go through large amounts of information relatively quickly. However, for the final report on the practicals, you will most likely want to use your own code to produce the necessary figures with appropriate annotation - you can of course adapt code contained in this module, if you wish.

I will add an example script as soon as I have clearance to include the associated data files, and more detailed documentation will eventually follow. For now though, read on for a brief description of what this module can do, and refer to the comments in the `.py` file itself for further details.

The module currently contains functions for handling two types of data:

### 1. Spatial data contained in a netCDF file

Start with the function `cm.list_vars`, which will simply print a list of the variables contained in the dataset, together with their full names, so that you can get a quick overview of the data you have.

Once you've found the name of a variable you want to explore further, try the following to extract the relevant data:

``` python
var, lons, lats = cm.get_var(dataset, 'variable_name')
```

You now have arrays containing the variable values, as well as the corresponding longitudes and latitudes (you will need these later to map the data).<br>
Note that there is no need to declare (create) new longitude and latitude variables every time you extract another variable from a dataset. If you already have these variables for the current dataset, you can get just the variable of interest using something like:

```python
var = cm.get_var(dataset, 'variable_name')[0]
```

Note that these variables usually have at least two dimensions, often more. Before you can produce a map of the data, you will have to reduce them to two dimensions - latitude and longitude only. Depending on what information you want, you might select an individual year of variable data or get an average over several years, to name two examples.<br>
Once you have two-dimensional (lat-lon) data, you can quickly draw a map of it using either `cm.map_data` (draws a rectangular map of the whole world) or `cm.plot_nearside` (draws a globe, centred on a user-specified coordinate).

### 2. Time series data

The main model used in the course is iLOVECLIM [[1]](#1). One of its outputs is a time series of daily data in a file called something like `book*****` (with no file extension). The function `cm.read_book` reads this file and returns a time series of the annual mean of the data. You can use `cm.plot_timeline` to graph the annual mean together with a 10-year running mean, to see how the variable changes over the model run.


## References

<a name="1"></a>
[1] Goosse, H., Brovkin, V., Fichefet, T., Haarsma, R., Huybrechts, P., Jongma, J., Mouchet, A., Selten, F., Barriat, P.-Y., Campin, J.-M., Deleersnijder, E., Driesschaert, E., Goelzer, H., Janssens, I., Loutre, M.-F., Morales Maqueda, M.A., Opsteegh, T., Mathieu, P.-P., Munhoven, G., Pettersson, E.J., Renssen, H., Roche, D.M., Schaeffer, M., Tartinville, B., Timmermann, A., Weber, S.L., 2010. Description of the Earth system model of intermediate complexity LOVECLIM version 1.2. Geoscientific Model Development 3, 603–633. [DOI: 10.5194/gmd-3-603-2010](https://doi.org/10.5194/gmd-3-603-2010)
