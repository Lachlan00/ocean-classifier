# ocean classifier - training data builder
"""
WARNING: This code is currentl;y very messey and hacked together
using method that are generally consdiered bad practice (e.g. 
setting global variables insdie functions). Be very cautiois if
extending out this code to other purposes.
"""

from netCDF4 import Dataset # reads netCDF file
from os import listdir
from os.path import isfile, join
import matplotlib.pyplot as plt
import pylab
from mpl_toolkits.basemap import Basemap # basemap tools
from datetime import datetime, timedelta # for working with datetimes
from random import randint
import numpy as np
import pandas as pd

# __FUNCTIONS__

# for getting time data
def grab_sst_time(time_idx):
    """
    gets datetime object for sst map projection
    """
    dtcon_days = time[time_idx]
    dtcon_start = datetime(1990,1,1) # This is the "days since" part
    dtcon_delta = timedelta(dtcon_days/24/60/60) # Create a time delta object from the number of days
    dtcon_offset = dtcon_start + dtcon_delta # Add the specified number of days to 1990
    frame_time = dtcon_offset
    return frame_time

def plot_temp(temp, time_idx, fig_no):
    """
    Make maps of temperature and salinity
    """

    # make frame__idx an integer to avoid slicing errors
    frame_idx = int(time_idx)

    # get 'frame_time'
    frame_time = grab_sst_time(frame_idx)

    # map setup
    fig = plt.figure()
    fig.subplots_adjust(left=0., right=1., bottom=0., top=0.9)
    # Setup the map
    m = Basemap(projection='merc', llcrnrlat=-38.050653, urcrnrlat=-34.453367,\
            llcrnrlon=147.996456, urcrnrlon=152.457344, lat_ts=20, resolution='h') # full range
    # draw stuff
    m.drawcoastlines() # comment out when using shapefile
    m.fillcontinents(color='black')
    # plot salt
    cs = m.pcolor(lons,lats,np.squeeze(temp), latlon = True ,vmin=temp_min, vmax=temp_max, cmap='plasma')
    # plot colourbar
    plt.colorbar()
    # datetime title
    plt.title('Regional - Temperature (Celcius)\n' + frame_time.strftime("%Y-%m-%d %H:%M:%S") + ' | ' + str(fname) + '_idx: ' + str(frame_idx))
    # stop axis from being cropped
    plt.tight_layout()

    # allow clicks to return lon/lat pairs
    lon, lat = 150.2269, -36.25201 # Location of Lighthouse
    xpt,ypt = m(lon,lat)
    lonpt, latpt = m(xpt,ypt,inverse=True)
    point, = m.plot(xpt,ypt,'+', markersize=12)

    annotation = plt.annotate('%5.1fW,%3.1fS' % (lon, lat), xy=(xpt,ypt),
             xytext=(20,35), textcoords="offset points", 
             bbox={"facecolor":"w", "alpha":0.3}, 
             arrowprops=dict(arrowstyle="->", connectionstyle="arc3"))

    # this is a hack method for now that should never be used in good code
    # sets all fucntion variables as global variables so 'm' can be called in 
    # the onlclick event. This is bad code practice. 
    globals().update(locals())

    return fig

def plot_salt(salt, time_idx, fig_no):
    """
    Make maps of temperature and salinity
    """

    # make frame__idx an integer to avoid slicing errors
    frame_idx = int(time_idx)

    # get 'frame_time'
    frame_time = grab_sst_time(frame_idx)

    # map setup
    fig = plt.figure()
    fig.subplots_adjust(left=0., right=1., bottom=0., top=0.9)
    # Setup the map
    m = Basemap(projection='merc', llcrnrlat=-38.050653, urcrnrlat=-34.453367,\
            llcrnrlon=147.996456, urcrnrlon=152.457344, lat_ts=20, resolution='h') # full range
    # draw stuff
    m.drawcoastlines() # comment out when using shapefile
    m.fillcontinents(color='black')
    # plot salt
    cs = m.pcolor(lons,lats,np.squeeze(salt), latlon = True ,vmin=salt_min, vmax=salt_max, cmap='viridis')
    # plot colourbar
    plt.colorbar()
    # datetime title
    plt.title('Regional - Salinity (PSU)\n' + frame_time.strftime("%Y-%m-%d %H:%M:%S") + ' | ' + str(fname) + '_idx: ' + str(frame_idx))
    # stop axis from being cropped
    plt.tight_layout()

    # add updating points
    lon, lat = 150.2269, -36.25201 # Location of Lighthouse

    # make salt map update with temp map
    xpt,ypt = m(lon,lat)
    lonpt, latpt = m(xpt,ypt,inverse=True)
    global point2
    point2, = m.plot(xpt,ypt,'+', markersize=12)

    global annotation2
    annotation2 = plt.annotate('%5.1fW,%3.1fS' % (lon, lat), xy=(xpt,ypt),
             xytext=(20,35), textcoords="offset points", 
             bbox={"facecolor":"w", "alpha":0.3}, 
             arrowprops=dict(arrowstyle="->", connectionstyle="arc3"))

    return fig


def onclick(event):
	ix, iy = event.xdata, event.ydata
	xpti, ypti = m(ix, iy,inverse=True)
	string = '(%5.1fW,%3.1fN)' % (xpti, ypti)
	print(string)
	annotation.xy = (ix, iy)
	point.set_data([ix], [iy])
	annotation.set_text(string)
	fig.canvas.draw_idle()
	annotation2.xy = (ix, iy)
	point2.set_data([ix], [iy])
	fig2.canvas.draw_idle()

	# convert to lat/lon
	# datlon, datlat = m(xpti,ypti,inverse=True)

	global train_data
	train_data.append((xpti, ypti, mass, time_value))



# __Setup__

# set colour scale variables
temp_min = 14
temp_max = 24
salt_min = 35.3
salt_max = 35.7

# list to hold data collected
train_data = []

# __Random Data_

# get list of files in data directory
directory = "/Users/lachlanphillips/PhD_Large_Data/ROMS/Montague_subset"
file_ls = [f for f in listdir(directory) if isfile(join(directory, f))]
file_ls = list(filter(lambda x:'naroom_avg' in x, file_ls))


# set randomness seed
plot_num = 1
np.random.seed(1010)
rnd_file = np.random.randint(len(file_ls), size=plot_num)
rnd_times = np.random.randint(29, size=plot_num)

# __Make plots__

for i in range(0, plot_num):
	# grab file
	file_no = rnd_file[i]
	file_path = directory + "/" + file_ls[file_no]
	fname = str(file_ls[i])[11:16]
	# grab time
	time_idx = rnd_times[i]
	fh = Dataset(file_path, mode='r')
    # extract data
	lats = fh.variables['lat_rho'][:]
	lons = fh.variables['lon_rho'][:]
	time = fh.variables['ocean_time'][:]
	temp = fh.variables['temp'][time_idx,29,:,:] 
	salt = fh.variables['salt'][time_idx,29,:,:] 

	# time output
	time_value = grab_sst_time(time_idx)

	# make interactive plot	
	fig = plot_temp(temp, time_idx, i)
	fig2 = plot_salt(salt, time_idx, i)
	cid = fig.canvas.mpl_connect('button_press_event', onclick)

	# Select EAC
	mass = 'EAC'
	print('Select EAC water')
	input("Press Enter to continue...")
	plt.show()

	# Rebuild plot
	fig = plot_temp(temp, time_idx, i)
	fig2 = plot_salt(salt, time_idx, i)
	cid = fig.canvas.mpl_connect('button_press_event', onclick)

	# Select BS 
	mass = 'BS'
	print('Select bass strait water')
	input("Press Enter to continue...")
	plt.show()

# convert to dataframe
print('__TEST__')
df = pd.DataFrame(train_data)
df.columns = ['lon', 'lat', 'class', 'datetime']
print(df)















