#!/usr/bin/env python3.9.5
# -*- Coding: UTF-8 -*-

from helpers.utilities_nwp import Utilities_nwp  # Plot Numeric Weather Models

# INPUT VARIABLES - extent [min. lon, min. lat, max. lon, max. lat]
# Initial forecast time (run model/initial conditions)
run_date = '2021072700'
# Number of hours since run to forecast
init_time = 0
end_time = 72
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
    extent = [-83.75, -55.75, -10.0, 14.3]
    for hour in fcst_hours:
        files = Utilities_nwp.download_wrf(path, run_date, '5', hour)

# Plot maps (4 level's maps) for each forecast
hour = 0
for file_in in files:
    print(f'Getting variables to map...: {file_in}')
    # Read general informations from data file - could be any var_name
    grib, lats, lons, valid = Utilities_nwp.grib_info(file_in, 'Temperature', extent)
    # Get 2 metre temperature
    temp = Utilities_nwp.get_temperature(grib, extent)
    # Get surface pressure to isobars
    prmls = Utilities_nwp.get_pressure(grib, extent)
    # Calculate 1000-500 hPa thickness
    thick = Utilities_nwp.get_thickness(grib, extent)
    # Get 850 level temperature
    temp850 = Utilities_nwp.get_temp850(grib, extent)
    # Get precipitable water
    pw = Utilities_nwp.get_pw(grib, extent)
    # Get vorticity - does not exist in WRF data
    #vort = Utilities_nwp.get_vort(grib, extent)
    # Get humidity
    humidity = Utilities_nwp.get_humidity(grib, extent)
    # Get/Calculate wind from level
    levels_hPa = [1000, 850, 500, 250]
    ucomp = {}
    vcomp = {}
    ws = {}
    for level in levels_hPa[1:]:
        u, v, w = Utilities_nwp.get_wind(grib, extent, level)
        ucomp[level] = u
        vcomp[level] = v
        ws[level] = w

    # Start plot 2x2
    fig, axs = Utilities_nwp.start_subplots(2,2)
    l = 0
    for i, j in zip([0,0,1,1], [0,1,0,1]):
        # Define subplot's level
        level = levels_hPa[l]
        print(f'- Making subplot: level {level} hPa')
        # Draw base map
        axs[i,j] = Utilities_nwp.base_map(axs[i, j], extent)
        # Subplot at [i,j] position
        if l == 0:
            title = f'{model}: 2m Temperature + 1000-500 hPa Thickness (m) + PMSL (hPa)'
            axs = Utilities_nwp.plot_sfc(axs, i, j, lats, lons, temp, prmls, thick, title, valid)
        elif l == 1:
            title = f'{model}: Streamlines/Barbs and Temperature (850 hPa) + Precipitable Water (kg/m²)'
            axs = Utilities_nwp.plot_850(axs, i, j, lats, lons, ucomp[850], vcomp[850], temp850, pw, title, valid)
        elif l == 2:
            title = f'{model}: Streamlines (500 hPa) and 2 metre specific humidity'
            axs = Utilities_nwp.plot_500(axs, i, j, lats, lons, ucomp[500], vcomp[500], humidity, title, valid)
        elif l == 3:
            title = f'{model}: Streamlines (250 hPa)'
            axs = Utilities_nwp.plot_250(axs, i, j, lats, lons, ucomp[250], vcomp[250], ws[250], title, valid)
        # Update levels counter
        l = l + 1

    # Save plot from forecast
    file_out = f'output/{model}_{run_date}_{str(hour).zfill(2)}.png'
    hour = hour + interval
    Utilities_nwp.save_plot(fig, file_out)
    #exit()
