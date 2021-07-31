#!/usr/bin/env python3.9.5
# -*- Coding: UTF-8 -*-

from helpers.utilities_nwp import Utilities_nwp  # Plot Numeric Weather Models

# INPUT VARIABLES - extent [min. lon, min. lat, max. lon, max. lat]
# Initial forecast time (run model/initial conditions)
run_date = '2021072700'
# Number of hours since run to forecast
init_time = 0
end_time = 0 #72
interval = 12
fcst_hours = list(range(init_time, end_time+1, interval))
# Select model: GFS or WRF
model = 'WRF'
path = f'input/{model}_{run_date}'
if model == 'GFS':
    extent = [-93.0, -60.00, -25.00, 18.00]
    for hour in fcst_hours:
        files = Utilities_nwp.download_gfs(path, run_date, extent, '50', hour)
elif model == 'WRF':
    #extent = [-83.75, -55.75, -10.0, 14.3]
    extent = [-60, -35, -43.0, -20]
    for hour in fcst_hours:
        files = Utilities_nwp.download_wrf(path, run_date, '5', hour)
# Variable name
var_name = '2 metre temperature'

hour = 0
for file_in in files:
    print(f'Getting variables to map...: {file_in}')
    # Get info from file - variable names
    #file_txt = f'output/variables_{model}.txt'
    #Utilities_nwp.get_info(file_in, file_txt)
    # Read general informations from data file
    grib, lats, lons, valid = Utilities_nwp.grib_info(file_in, var_name, extent)
    # Get 2 metre temperature
    temp = Utilities_nwp.get_temperature(grib, extent)

    # Smooth
    #grb = Utilities_nwp.zoom(temp)

    # Plot map
    print(f'Plotting map...')
    # Define properties to plot map - 'None' to auto
    properties = {
        'colormap': 'jet',
        'label': '2 m Temperature (°C)',
        'model': model,
        'vmin': -15,
        'vmax': 32
    }
    file_out = f'output/{model}_{run_date}_{str(hour).zfill(2)}_var.png'
    Utilities_nwp.plot_map(temp, valid, extent, file_out, properties)
    hour = hour + interval
    #exit()
