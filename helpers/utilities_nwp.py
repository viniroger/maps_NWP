#!/usr/bin/env python3.9.5
# -*- Coding: UTF-8 -*-

import glob
from datetime import datetime, timedelta   # Basic Dates and time types
import os                                  # Operating system interfaces
import requests                            # HTTP library for Python
import numpy as np                         # Scientific computing with Python
import pygrib                              # Reading GRIB files
import matplotlib.pyplot as plt            # Plotting library
from matplotlib.axes import Axes           # Axes for plot
from matplotlib.colors import ListedColormap # Customized colorbar
import cartopy, cartopy.crs as ccrs        # Plot maps
import cartopy.io.shapereader as shpreader # Import shapefiles

class Utilities_nwp():

    @staticmethod
    def download(path, file_name, url):
        '''
        Download file, if doesn't exist
        '''
        # Check/create directory
        os.makedirs(path, exist_ok=True)
        # Check if file exists
        path_file = f'{path}/{file_name}'
        #print('URL:', url)
        print('File name:', path_file)
        if os.path.exists(path_file):
            print('File exists')
        else:
            print('Downloading file...')
            # Sends a GET request to the specified url
            myfile = requests.get(url)
            # Download the file
            x = f'{path}//{file_name}'
            open(x, 'wb').write(myfile.content)

    @classmethod
    def download_gfs(self, path, run_date, extent, resolution, hour):
        '''
        Download GFS files from NOMADS
        '''
        # Create the URL (add 360 for longitudes)
        hour_run = run_date[-2:]
        min_lon = str(int(extent[0]+360)) # leftlon
        min_lat = str(int(extent[1])) # bottomlat
        max_lon = str(int(extent[2]+360)) # rightlon
        max_lat = str(int(extent[3])) # toplat
        date = run_date[:8]
        if (resolution == '25'):
            url = 'https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p'\
            +resolution+'.pl?file=gfs.t'+hour_run+'z.pgrb2.0p'+resolution\
            +'.f'+str(hour).zfill(3)\
            +'&all_lev=on&all_var=on&subregion=&leftlon='+min_lon\
            +'&rightlon='+max_lon+'&toplat='+max_lat+'&bottomlat='\
            +min_lat+'&dir=%2Fgfs.'+date+'%2F00%2Fatmos'
            file_name = 'gfs.t'+hour_run+'z.pgrb2.0p'+resolution+'.f'+str(hour).zfill(3)
        elif (resolution == '50'):
            url = 'https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p'\
            +resolution+'.pl?file=gfs.t'+hour_run+'z.pgrb2full.0p'+resolution\
            +'.f'+str(hour).zfill(3)\
            +'&all_lev=on&all_var=on&subregion=&leftlon='+min_lon\
            +'&rightlon='+max_lon+'&toplat='+max_lat+'&bottomlat='\
            +min_lat+'&dir=%2Fgfs.'+date+'%2F00%2Fatmos'
            file_name = 'gfs.t'+hour_run+'z.pgrb2.0p'+resolution+'.f'+str(hour).zfill(3)
        elif (resolution == '1'):
            url = 'https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_'\
            +resolution+'p00.pl?file=gfs.t'+hour_run+'z.pgrb2.'+resolution\
            +'p00.f'+str(hour).zfill(3)\
            +'&all_lev=on&all_var=on&subregion=&leftlon='+min_lon+'&rightlon='\
            +max_lon+'&toplat='+max_lat+'&bottomlat='+min_lat+'&dir=%2Fgfs.'\
            +date+'%2F00%2Fatmos'
            file_name = 'gfs.t'+hour_run+'z.pgrb2.'+resolution+'p00.f'+str(hour).zfill(3)
        # Download file, if doesn't exist
        self.download(path, file_name, url)
        # List all files with pattern at directory
        pattern = f'{path}/*'
        files = sorted([f for f in glob.glob(pattern, recursive=True)])
        return files

    @classmethod
    def download_wrf(self, path, run_date, resolution, hour):
        '''
        Download WRF files from INPE
        '''
        # Check/create directory
        os.makedirs(path, exist_ok=True)
        # Create the URL
        res = resolution.zfill(2)
        datetime_fmt = datetime.strptime(run_date, '%Y%m%d%H')
        run_date_path = datetime_fmt.strftime('%Y/%m/%d/%H')
        fsct_date = datetime_fmt + timedelta(hours=hour)
        fsct_date = fsct_date.strftime('%Y%m%d%H')
        file_name = f'WRF_cpt_{res}KM_{run_date}_{fsct_date}.grib2'
        url = f'http://ftp.cptec.inpe.br/modelos/tempo/WRF/'\
        f'ams_{res}km/brutos/{run_date_path}/{file_name}'
        # Download file, if doesn't exist
        self.download(path, file_name, url)
        # List all files with pattern at directory
        pattern = f'{path}/*'
        files = sorted([f for f in glob.glob(pattern, recursive=True)])
        return files

    @staticmethod
    def get_info(file_in, file_out):
        '''
        Get variables names from file_in and save to file_out
        '''
        # Open the GRIB file
        grib = pygrib.open(file_in)
        # Print all variables in the terminal and save them in a text file
        f = open(file_out, 'w')
        for variables in grib:
            # Print the variables in the terminal
            print(variables)
            # Put the variables in the text file
            print(variables, file=f)
        f.close()

    @staticmethod
    def smooth(data):
        '''
        Smooth
        '''
        print("\nArray dimensions before smoothing:")
        print(data.shape)
        import scipy.ndimage
        data = scipy.ndimage.zoom(data, 3)
        lats = scipy.ndimage.zoom(lats, 3)
        lons = scipy.ndimage.zoom(lons, 3)
        print("Array dimensions after smoothing:")
        print(data.shape)
        return data

    @staticmethod
    def plot_map(data, valid, extent, file_out, properties):
        '''
        Plot map
        '''
        # Choose the plot size (width x height, in inches)
        plt.figure(figsize=(10,10))
        # Use the Cilindrical Equidistant projection in cartopy
        ax = plt.axes(projection=ccrs.PlateCarree())
        # Define the image extent [min. lon, max. lon, min. lat, max. lat]
        img_extent = [extent[0], extent[2], extent[1], extent[3]]

        # Plot the image
        plt.imshow(data, origin='lower', extent=img_extent,\
         vmin=properties['vmin'], vmax=properties['vmax'],\
         cmap=properties['colormap'])
        # Add a shapefile
        shapefile = list(shpreader.Reader('helpers/BR_UF_2019.shp').geometries())
        ax.add_geometries(shapefile, ccrs.PlateCarree(), edgecolor='black',\
         facecolor='none', linewidth=0.3)

        # Add coastlines, borders and gridlines
        ax.coastlines(resolution='10m', color='black', linewidth=0.8)
        ax.add_feature(cartopy.feature.BORDERS, edgecolor='black', linewidth=0.5)
        gl = ax.gridlines(crs=ccrs.PlateCarree(), color='gray', alpha=1.0,\
         linestyle='--', linewidth=0.25, xlocs=np.arange(-180, 180, 5),\
          ylocs=np.arange(-90, 90, 5), draw_labels=True)
        gl.top_labels = False
        gl.right_labels = False

        # Add a colorbar
        plt.colorbar(label=properties['label'], extend='both',\
         orientation='horizontal', pad=0.05, fraction=0.05)
        # Add a title
        title = f'{properties["model"]}: {properties["label"]}'
        plt.title(title , fontweight='bold', fontsize=10, loc='left')
        plt.title('Valid: ' + valid, fontsize=10, loc='right')

        # Save the image
        os.makedirs('output', exist_ok=True)
        plt.savefig(file_out, bbox_inches='tight')
        # Show the image
        #plt.show()

    @staticmethod
    def plot_maxmin_points(ax, lon, lat, data, extrema, nsize, symbol, color='k',
                       plotValue=True, transform=None):
        """
        This function will find and plot relative maximum and minimum for a 2D grid. The function
        can be used to plot an H for maximum values (e.g., High pressure) and an L for minimum
        values (e.g., low pressue). It is best to used filetered data to obtain  a synoptic scale
        max/min value. The symbol text can be set to a string value and optionally the color of the
        symbol and any plotted value can be set with the parameter color
        lon = plotting longitude values (2D)
        lat = plotting latitude values (2D)
        data = 2D data that you wish to plot the max/min symbol placement
        extrema = Either a value of max for Maximum Values or min for Minimum Values
        nsize = Size of the grid box to filter the max and min values to plot a reasonable number
        symbol = String to be placed at location of max/min value
        color = String matplotlib colorname to plot the symbol (and numerica value, if plotted)
        plot_value = Boolean (True/False) of whether to plot the numeric value of max/min point
        The max/min symbol will be plotted on the current axes within the bounding frame
        (e.g., clip_on=True)
        """
        from scipy.ndimage.filters import maximum_filter, minimum_filter

        if (extrema == 'max'):
            data_ext = maximum_filter(data, nsize, mode='nearest')
        elif (extrema == 'min'):
            data_ext = minimum_filter(data, nsize, mode='nearest')
        else:
            raise ValueError('Value for hilo must be either max or min')

        mxy, mxx = np.where(data_ext == data)

        for i in range(len(mxy)):
             txt1 = ax.annotate(symbol, xy=(lon[mxy[i], mxx[i]], lat[mxy[i], mxx[i]]), xycoords=ccrs.PlateCarree()._as_mpl_transform(ax), color=color, size=24,
                    clip_on=True, annotation_clip=True, horizontalalignment='center', verticalalignment='center',
                    transform=ccrs.PlateCarree())

             txt2 = ax.annotate('\n' + str(int(data[mxy[i], mxx[i]])), xy=(lon[mxy[i], mxx[i]], lat[mxy[i], mxx[i]]), xycoords=ccrs.PlateCarree()._as_mpl_transform(ax),
                    color=color, size=12, clip_on=True, annotation_clip=True, fontweight='bold', horizontalalignment='center', verticalalignment='top',
                    transform=ccrs.PlateCarree())

    @staticmethod
    def grib_info(file_in, var_name, extent):
        '''
        Read grib file and get information
        '''
        # Open the GRIB file
        grib = pygrib.open(file_in)
        # Select the variable
        grb = grib.select(name=var_name)[0]
        # Get information from the file
        init  = str(grb.analDate)      # Init date / time
        run   = str(grb.hour).zfill(2) # Run
        ftime = str(grb.forecastTime)  # Forecast hour
        valid = str(grb.validDate)     # Valid date / time
        #print('Init: ' + init + ' UTC')
        #print('Run: ' + run + 'Z')
        #print('Forecast: +' + ftime)
        #print('Valid: ' + valid + ' UTC')
        # Get latitudes and longitudes
        grb, lats, lons = grb.data(lat1=extent[1], lat2=extent[3],\
         lon1=extent[0]+360, lon2=extent[2]+360)
        return grib, lats, lons, valid

    @staticmethod
    def get_temperature(grib, extent):
        '''
        Get 2 metre temperature and convert from Kelvin to Celsius
        '''
        # Select the variable
        temp = grib.select(name='2 metre temperature')[0]
        # Read the data for a specific region
        temp, lats, lons = temp.data(lat1=extent[1],lat2=extent[3],lon1=extent[0]+360,lon2=extent[2]+360)
        # Convert from K to °C
        temp = temp - 273.15
        return temp

    @staticmethod
    def get_temp850(grib, extent):
        '''
        Get 850 level temperature and convert from Kelvin to Celsius
        '''
        # Select the variable
        temp = grib.select(name='Temperature', typeOfLevel = 'isobaricInhPa', level = 850)[0]
        # Read the data for a specific region
        temp = temp.data(lat1=extent[1],lat2=extent[3],lon1=extent[0]+360,lon2=extent[2]+360)[0]
        # Convert from K to °C
        temp = temp - 273.15
        return temp

    @staticmethod
    def get_pressure(grib, extent):
        '''
        Get Pressure reduced to MSL in Pa and convert to hPa
        '''
        # Select the variable
        prmls = grib.select(name='Pressure reduced to MSL')[0]
        # Read the data for a specific region
        prmls, lats, lons = prmls.data(lat1=extent[1],lat2=extent[3],lon1=extent[0]+360,lon2=extent[2]+360)
        # Convert to hPa
        prmls = prmls / 100
        return prmls

    @staticmethod
    def get_thickness(grib, extent):
        '''
        Get Geopotential Height at 1000 and 500 hPa and
        Calculate 1000-500 hPa thickness
        '''
        # Select the variable
        hght_1000 = grib.select(name='Geopotential Height', typeOfLevel = 'isobaricInhPa', level = 1000)[0]
        # Read the data for a specific region
        hght_1000 = hght_1000.data(lat1=extent[1],lat2=extent[3],lon1=extent[0]+360,lon2=extent[2]+360)[0]

        # Select the variable
        hght_500 = grib.select(name='Geopotential Height', typeOfLevel = 'isobaricInhPa', level = 500)[0]
        # Read the data for a specific region
        hght_500 = hght_500.data(lat1=extent[1],lat2=extent[3],lon1=extent[0]+360,lon2=extent[2]+360)[0]

        # Calculate and smooth 1000-500 hPa thickness
        thickness_1000_500 = hght_500 - hght_1000
        return thickness_1000_500

    @staticmethod
    def get_pw(grib, extent):
        '''
        Get precipitable water
        '''
        # Select the variable
        pw = grib.select(name='Precipitable water')[0]
        # Read the data for a specific region
        pw = pw.data(lat1=extent[1],lat2=extent[3],lon1=extent[0]+360,lon2=extent[2]+360)[0]
        return pw

    @staticmethod
    def get_vort(grib, extent):
        '''
        Get vorticity
        '''
        # Select the variable
        vort_500 = grib.select(name='Absolute vorticity', typeOfLevel = 'isobaricInhPa', level = 500)[0]
        # Read the data for a specific region
        vort_500 = vort_500.data(lat1=extent[1],lat2=extent[3],lon1=extent[0]+360,lon2=extent[2]+360)[0]
        # Scale: 10^4
        vort_500 = vort_500 * 10000
        return vort_500

    @staticmethod
    def get_humidity(grib, extent):
        '''
        Get Specific humidity or 2 metre specific humidity
        '''
        # Select the variable
        #humidity = grib.select(name='Specific humidity', typeOfLevel = 'isobaricInhPa', level = 500)[0]
        humidity = grib.select(name='2 metre specific humidity')[0]
        # Read the data for a specific region
        humidity = humidity.data(lat1=extent[1],lat2=extent[3],lon1=extent[0]+360,lon2=extent[2]+360)[0]
        # Scale: 10^2
        humidity = humidity * 100
        return humidity

    @staticmethod
    def get_wind(grib, extent, level):
        '''
        Get wind (u,v) from grib file and calculate w
        '''
        # Select the variable - u
        u = grib.select(name='U component of wind',\
         typeOfLevel = 'isobaricInhPa', level = level)[0]
        # Read the data for a specific region
        u = u.data(lat1=extent[1], lat2=extent[3],\
         lon1=extent[0]+360, lon2=extent[2]+360)[0]

        # Select the variable - v
        v = grib.select(name='V component of wind',\
         typeOfLevel = 'isobaricInhPa', level = level)[0]
        # Read the data for a specific region
        v = v.data(lat1=extent[1], lat2=extent[3],\
         lon1=extent[0]+360, lon2=extent[2]+360)[0]

        # Calculate the wind speed
        ws = np.sqrt(u**2 + v**2)

        return u, v, ws

    @staticmethod
    def custom_colormap():
        '''
        Customized colormap
        '''
        # Define de contour interval
        data_min = 0
        data_max = 60
        interval = 5
        levels = np.arange(data_min,data_max,interval)

        # Create a custom color palette
        colors = ["#e7f2f4", "#ceeaee", "#b6e2e8", "#abdcff", "#a4d685", "#9cd04e",
                  "#abcf2a", "#c9d21b", "#e8d50c", "#ffd100", "#ffba00", "#ffa200"]
        cmap = ListedColormap(colors)
        cmap.set_over('#ff8c00')
        cmap.set_under('#fffafa')
        return cmap, levels

    @staticmethod
    def start_subplots(nrows, ncols):
        '''
        Start plot nrows x ncols
        '''
        fig, axs = plt.subplots(nrows, ncols, figsize=(20,20),\
         sharex = False, sharey = False, gridspec_kw ={'left':0, 'bottom':0,\
         'right':1, 'top':1, 'hspace':0, 'wspace':0.1},\
         subplot_kw=dict(projection=ccrs.PlateCarree()))
        return fig, axs

    @staticmethod
    def base_map(ax, extent):
        '''
        Draw map based at extent and using:
        - shapefile with brazilian states borders
        - coastlines, borders and gridlines
        '''
        # Define the image extent
        ax.set_extent([extent[0], extent[2], extent[1], extent[3]], ccrs.PlateCarree())

        # Add a shapefile
        shapefile = list(shpreader.Reader('helpers/BR_UF_2019.shp').geometries())
        ax.add_geometries(shapefile, ccrs.PlateCarree(), edgecolor='black',\
         facecolor='none', linewidth=0.3, linestyles='dashed')

        # Add coastlines, borders and gridlines
        ax.coastlines(resolution='10m', color='black', linewidth=0.8)
        ax.add_feature(cartopy.feature.BORDERS, edgecolor='black', linewidth=0.5)
        gl = ax.gridlines(crs=ccrs.PlateCarree(), color='gray', alpha=1.0,\
         linestyle='--', linewidth=0.25, xlocs=np.arange(-180, 180, 5),\
          ylocs=np.arange(-90, 90, 5), draw_labels=True)
        gl.top_labels = False
        gl.right_labels = False

        return ax

    @classmethod
    def plot_sfc(self, axs, i, j, lats, lons, temp, prmls, thick, title, valid):
        '''
        Plot Surface/1000 hPa:
        - 2 metre temperature
        - 1000-500 hPa thickness
        - Isobars with Higs/Lows
        '''
        # TEMPERATURE
        # Define de contour interval
        data_min = -15
        data_max = 30
        interval = 5
        levels = np.arange(data_min,data_max,interval)
        # Plot the contours
        img1 = axs[i,j].contourf(lons, lats, temp,\
         cmap='jet', levels=levels, extend='both', alpha = 0.6)
        plt.colorbar(img1, label='2 m Temperature (°C)',\
         orientation='vertical', pad=0.02, fraction=0.05, ax=axs[i,j], shrink=0.85)

        # THICKNESS
        # Define de contour interval
        data_min = 4900
        data_max = 5900
        interval = 20
        levels = np.arange(data_min,data_max,interval)
        # Plot the contours
        img2 = axs[i,j].contour(lons, lats, thick, cmap='seismic',\
         linestyles='dashed', linewidths=1.0, levels=levels)
        axs[i,j].clabel(img2, inline=1, inline_spacing=0, fontsize='10',\
         fmt = '%1.0f')#, colors= 'black')
        # Get the index of elements with value "5400"
        mid_value = int(np.where(levels == 5400)[0])
        img2.collections[mid_value].set_linewidth(4)
        img2.collections[mid_value].set_color('blue')

        # PRESSURE
        # Define de contour interval
        data_min = 500
        data_max = 1050
        interval = 2
        levels = np.arange(data_min,data_max,interval)
        # Plot the contours
        img3 = axs[i,j].contour(lons, lats, prmls, colors='black',\
         linewidths=1, levels=levels)
        axs[i,j].clabel(img3, inline=1, inline_spacing=0, fontsize='10',\
         fmt = '%1.0f', colors= 'black')
         # Define nsize to plot_maxmin_points considering higher WRF definition
        if 'WRF' in title:
            nsizeH = 300
            nsizeL = 250
        else:
            nsizeH = 50
            nsizeL = 25
        # Use definition to plot H/L symbols
        self.plot_maxmin_points(axs[i,j], lons, lats, prmls, 'max', nsizeH,\
         symbol='H', color='b', transform=ccrs.PlateCarree())
        self.plot_maxmin_points(axs[i,j], lons, lats, prmls, 'min', nsizeL,\
         symbol='L', color='r', transform=ccrs.PlateCarree())

        # Add a title
        axs[i,j].set_title(title, fontweight='bold', fontsize=10, loc='left')
        axs[i,j].set_title('Valid: ' + valid + ' UTC', fontsize=10, loc='right')

        return axs

    @staticmethod
    def plot_850(axs, i, j, lats, lons, u, v, temp, pw, title, valid):
        '''
        Plot 850 hPa:
        - Streamlines and Barbs
        - Temperature lines
        - Precipitable water
        '''
        # STREAMLINES
        img1 = Axes.streamplot(axs[i,j], lons, lats, u, v, density=[4, 4],\
         linewidth=1, color='gray', transform=ccrs.PlateCarree())

        # BARBS
        # Create a flag to determine which barbs are flipped
        flip_flag = np.zeros((u.shape[0],u.shape[1]))
        # All flags below the equator will be flipped
        flip_flag[lats < 0] = 1
        # Plot the barbs - consider that WRF is higher resolution
        if 'WRF' in title:
            n = 64
        else:
            n = 8
        img3 = axs[i,j].barbs(lons[::n,::n], lats[::n,::n], u[::n,::n],\
         v[::n,::n], length = 7.0,\
         sizes = dict(emptybarb=0.0, spacing=0.2, height=0.5),\
         linewidth=1.5, pivot='middle', barbcolor='orange',\
         flip_barb = flip_flag[::n,::n])

        # PRECIPITABLE WATER
        # Define de contour interval
        data_min = 10
        data_max = 60
        interval = 10
        levels = np.arange(data_min,data_max,interval)
        # Plot the contours
        img2 = axs[i,j].contourf(lons, lats, pw,\
         cmap='Blues', levels = levels, extend='both', alpha = 0.5)
        plt.colorbar(img2, label='Precipitable water (kg/m²)',\
         orientation='vertical', pad=0.02, fraction=0.05, ax=axs[i,j], shrink=0.85)

        # TEMPERATURE
        # Define de contour interval
        data_min = -20
        data_max = 20
        interval = 2
        levels = np.arange(data_min,data_max,interval)
        # Plot the contours
        img3 = axs[i,j].contour(lons, lats, temp, cmap='seismic',\
         linewidths=1, levels=levels)
        axs[i,j].clabel(img3, inline=1, inline_spacing=0, fontsize='10',\
         fmt = '%1.0f')
        # Get the index of elements with value "0"
        mid_value = int(np.where(levels == 0)[0])
        img3.collections[mid_value].set_linewidth(2)
        img3.collections[mid_value].set_color('blue')

        # Add a title
        axs[i,j].set_title(title, fontweight='bold', fontsize=10, loc='left')
        axs[i,j].set_title('Valid: ' + valid + ' UTC', fontsize=10, loc='right')

        return axs

    @classmethod
    def plot_500(self, axs, i, j, lats, lons, u, v, humidity, title, valid):
        '''
        Plot 500 hPa:
        - Streamlines
        - Vorticity or Humidity
        '''
        # STREAMLINES
        img1 = Axes.streamplot(axs[i,j], lons, lats, u, v, density=[4, 4],\
         linewidth=1, color='gray', transform=ccrs.PlateCarree())

        # VORTICITY
        #img2 = axs[i,j].contourf(lons, lats, vort,\
        # cmap='bwr', extend='both', alpha = 0.9)
        #plt.colorbar(img2, label='Vorticity (10^-4 s-1)',\
        # orientation='vertical', pad=0.02, fraction=0.05, ax=axs[i,j], shrink=0.85)

        # HUMIDITY
        img2 = axs[i,j].contourf(lons, lats, humidity,\
         cmap='rainbow_r', extend='both', alpha = 0.9)
        plt.colorbar(img2, label='2 metre specific humidity (10^-2 kg/kg)',\
         orientation='vertical', pad=0.02, fraction=0.05, ax=axs[i,j], shrink=0.85)

        # Add a title
        axs[i,j].set_title(title, fontweight='bold', fontsize=10, loc='left')
        axs[i,j].set_title('Valid: ' + valid + ' UTC', fontsize=10, loc='right')

        return axs

    @classmethod
    def plot_250(self, axs, i, j, lats, lons, u, v, w, title, valid):
        '''
        Plot 250 hPa:
        - Streamlines
        - Wind fields
        '''
        cmap, levels = self.custom_colormap()
        # Plot the contours
        img1 = axs[i,j].contourf(lons, lats, w, cmap=cmap, levels=levels,\
         extend='both')
        img2 = axs[i,j].contour(lons, lats, w, colors='white', linewidths=0.3,\
         levels=levels)
        axs[i,j].clabel(img2, inline=1, inline_spacing=0, fontsize='10',\
         fmt = '%1.0f', colors= 'black')
        # Add a colorbar
        plt.colorbar(img1, label='Isotachs (kt)', orientation='vertical',\
         pad=0.02, fraction=0.05, ax=axs[i,j], shrink=0.85)

        # STREAMLINES
        img3 = Axes.streamplot(axs[i,j], lons, lats, u, v, density=[4, 4],\
         linewidth=1, color='gray', transform=ccrs.PlateCarree())

        # Add a title
        axs[i,j].set_title(title, fontweight='bold', fontsize=10, loc='left')
        axs[i,j].set_title('Valid: ' + valid + ' UTC', fontsize=10, loc='right')

        return axs

    @staticmethod
    def save_plot(fig, file_out):
        '''
        Save plot
        '''
        # Check/create directory
        path = ''.join(file_out.split('/')[:-1])
        os.makedirs(path, exist_ok=True)
        plt.savefig(file_out, bbox_inches='tight')
