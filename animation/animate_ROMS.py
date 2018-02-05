"""
Produce animation of ROMS data 
"""
print('----------------------------')
print('NetCDF animation builder v-2.0')
print('----------------------------')
print('----Program Start----')

print('loading packages...')
#Packages
from netCDF4 import Dataset # reads netCDF file
import numpy as np # manipulates arrays
import pandas as pd # for dataframes and reading track csv data
import matplotlib.pyplot as plt # for plotting map
from mpl_toolkits.basemap import Basemap # basemap tools
from datetime import datetime, timedelta #for working with datetimes
import moviepy.editor as mpy # creates animation
from moviepy.video.io.bindings import mplfig_to_npimage # converts map to numpy array
from matplotlib.backends.backend_agg import FigureCanvasAgg # draws canvas so that map can be converted

# define animation buiding functions
##########__FUNCTIONS__##########

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

def make_frame_temp(frame_idx):
    """
    Make animation of just temperature
    """
    # set start frame
    start_frame = 0

    # make frame__idx an integer to avoid slicing errors
    frame_idx = int(frame_idx)

    # import salt at timestamp
    temp = fh.variables['temp'][start_frame + frame_idx,29,:,:] 

    # get 'frame_time'
    frame_time = grab_sst_time(start_frame + frame_idx)
   
    # map setup
    fig = plt.figure()
    fig.subplots_adjust(left=0., right=1., bottom=0., top=0.9)
    # Setup the map
    m = Basemap(projection='merc', llcrnrlat=-38.050653, urcrnrlat=-34.453367,\
            llcrnrlon=147.996456, urcrnrlon=152.457344, lat_ts=20, resolution='h')
    # draw stuff
    m.drawcoastlines()
    m.fillcontinents(color='black')
    # plot salt
    cs = m.pcolor(lons,lats,np.squeeze(temp), latlon = True ,vmin=temp_min, vmax=temp_max, cmap='plasma')
    # plot colourbar
    plt.colorbar()
    # datetime title
    plt.title('Regional - Temperature (Celcius)\n' + frame_time.strftime("%Y-%m-%d %H:%M:%S") + ' | ' + str(fname) + '_idx: ' + str(frame_idx))
    # stop axis from being cropped
    plt.tight_layout() 
    
    #convert to array
    canvas = FigureCanvasAgg(fig)
    canvas.draw()
    frame = np.fromstring(canvas.tostring_rgb(), dtype='uint8')
    frame = frame.reshape(fig.canvas.get_width_height()[::-1] + (3,))
    
    return frame

def make_frame_salt(frame_idx):
    """
    Make animation of just salinity
    """
    # set start frame
    start_frame = 0

    # make frame__idx an integer to avoid slicing errors
    frame_idx = int(frame_idx)

    # import salt at timestamp
    salt = fh.variables['salt'][start_frame + frame_idx,29,:,:] 

    # get 'frame_time'
    frame_time = grab_sst_time(start_frame + frame_idx)
   
    # map setup
    fig = plt.figure()
    fig.subplots_adjust(left=0., right=1., bottom=0., top=0.9)
    # Setup the map
    m = Basemap(projection='merc', llcrnrlat=-38.050653, urcrnrlat=-34.453367,\
            llcrnrlon=147.996456, urcrnrlon=152.457344, lat_ts=20, resolution='h')
    # draw stuff
    m.drawcoastlines()
    m.fillcontinents(color='black')
    # plot salt
    cs = m.pcolor(lons,lats,np.squeeze(salt), latlon = True ,vmin=salt_min, vmax=salt_max, cmap='viridis')
    # plot colourbar
    plt.colorbar()
    # datetime title
    plt.title('Regional - Salinity (PSU)\n' + frame_time.strftime("%Y-%m-%d %H:%M:%S") + ' | ' + str(fname) + '_idx: ' + str(frame_idx))
    # stop axis from being cropped
    plt.tight_layout() 
    
    #convert to array
    canvas = FigureCanvasAgg(fig)
    canvas.draw()
    frame = np.fromstring(canvas.tostring_rgb(), dtype='uint8')
    frame = frame.reshape(fig.canvas.get_width_height()[::-1] + (3,))
    
    return frame

# import file
print('loading data...')
nc_file = '/Users/lachlanphillips/PhD_Large_Data/ROMS/Montague_subset/naroom_avg_02601.nc'
fname = '02601'
fh = Dataset(nc_file, mode='r')

#extract lats, lons and time
lats = fh.variables['lat_rho'][:] 
lons = fh.variables['lon_rho'][:]
time = fh.variables['ocean_time'][:]

# set colour scale variables
temp_min = 14
temp_max = 24 
salt_min = 35.1
salt_max = 35.7

# make animations
output_name = 'test_animation_salt.gif'
animation = mpy.VideoClip(make_frame_salt, duration=30)
animation.write_gif(output_name, fps=1)

output_name2 = 'test_animation_temp.gif'
animation = mpy.VideoClip(make_frame_temp, duration=30)
animation.write_gif(output_name2, fps=1)
#animation.write_videofile('animation.mp4', fps=1)

# close nc file
print('closing file...')
fh.close()
print('-----Program End-----')

print('')
print('´´´´´´´´´´´´´´´´´´´´´´¶¶¶¶¶¶¶¶¶……..')
print('´´´´´´´´´´´´´´´´´´´´¶¶´´´´´´´´´´¶¶……')
print('´´´´´´¶¶¶¶¶´´´´´´´¶¶´´´´´´´´´´´´´´¶¶……….')
print('´´´´´¶´´´´´¶´´´´¶¶´´´´´¶¶´´´´¶¶´´´´´¶¶…………..')
print('´´´´´¶´´´´´¶´´´¶¶´´´´´´¶¶´´´´¶¶´´´´´´´¶¶…..')
print('´´´´´¶´´´´¶´´¶¶´´´´´´´´¶¶´´´´¶¶´´´´´´´´¶¶…..')
print('´´´´´´¶´´´¶´´´¶´´´´´´´´´´´´´´´´´´´´´´´´´¶¶….')
print('´´´´¶¶¶¶¶¶¶¶¶¶¶¶´´´´´´´´´´´´´´´´´´´´´´´´¶¶….')
print('´´´¶´´´´´´´´´´´´¶´¶¶´´´´´´´´´´´´´¶¶´´´´´¶¶….')
print('´´¶¶´´´´´´´´´´´´¶´´¶¶´´´´´´´´´´´´¶¶´´´´´¶¶….')
print('´¶¶´´´¶¶¶¶¶¶¶¶¶¶¶´´´´¶¶´´´´´´´´¶¶´´´´´´´¶¶…')
print('´¶´´´´´´´´´´´´´´´¶´´´´´¶¶¶¶¶¶¶´´´´´´´´´¶¶….')
print('´¶¶´´´´´´´´´´´´´´¶´´´´´´´´´´´´´´´´´´´´¶¶…..')
print('´´¶´´´¶¶¶¶¶¶¶¶¶¶¶¶´´´´´´´´´´´´´´´´´´´¶¶….')
print('´´¶¶´´´´´´´´´´´¶´´¶¶´´´´´´´´´´´´´´´´¶¶….')
print('´´´¶¶¶¶¶¶¶¶¶¶¶¶´´´´´¶¶´´´´´´´´´´´´¶¶…..')
print('´´´´´´´´´´´´´´´´´´´´´´´¶¶¶¶¶¶¶¶¶¶¶…….)')
print('')
print('Tip: To increase fame speed, use the following commands to use imagemagick from terminal:')
print('convert -delay 25x100 ' + output_name + ' ' + output_name)
print('convert -delay 25x100 ' + output_name2 + ' ' + output_name2)
print('NetCDF animation builder v-2.0 does not currently support subplots.')
print('For now use "https://ezgif.com/combine" to join the gifs.')

# Additional function to make both plots on one GIF. Needs to be debugged. 
def make_frame_dual(frame_idx):
    """
    Make animation of just salinity
    NOTE: Until I work out the solution for plotting colorbars using `ax` 
    variables this website is a great alternative solution to combine the gifs: 
    https://ezgif.com/combine
    """
    # set start frame
    start_frame = 0

    # make frame__idx an integer to avoid slicing errors
    frame_idx = int(frame_idx)

    # import salt and temp at timestamp
    salt = fh.variables['salt'][start_frame + frame_idx,29,:,:] 
    temp = fh.variables['temp'][start_frame + frame_idx,29,:,:]

    # get 'frame_time'
    frame_time = grab_sst_time(start_frame + frame_idx)

    # set up figure
    fig = plt.figure()
    fig.subplots_adjust(left=0., right=1., bottom=0., top=0.9)

    # Temperature figure
    plt.subplot(1, 2, 1)
    # fig1.subplots_adjust(left=0., right=1., bottom=0., top=0.9)
    # Setup the map
    m = Basemap(projection='merc', llcrnrlat=-38.050653, urcrnrlat=-34.453367,\
            llcrnrlon=147.996456, urcrnrlon=152.457344, lat_ts=20, resolution='h')
    # draw stuff
    m.drawcoastlines()
    m.fillcontinents(color='black')
    # plot salt
    cs = m.pcolor(lons,lats,np.squeeze(temp), latlon = True ,vmin=temp_min, vmax=temp_max, cmap='plasma')
    # plot colourbar
    plt.colorbar()
    # datetime title
    plt.title('Regional - Temperature (Celcius)\n' + frame_time.strftime("%Y-%m-%d %H:%M:%S") + ' | ' + str(fname) + '_idx: ' + str(frame_idx))
   
    # Salinity figure
    plt.subplot(1, 2, 2)
    # fig2.subplots_adjust(left=0., right=1., bottom=0., top=0.9)
    # Setup the map
    m = Basemap(projection='merc', llcrnrlat=-38.050653, urcrnrlat=-34.453367,\
            llcrnrlon=147.996456, urcrnrlon=152.457344, lat_ts=20, resolution='h')
    # draw stuff
    m.drawcoastlines()
    m.fillcontinents(color='black')
    # plot salt
    cs = m.pcolor(lons,lats,np.squeeze(salt), latlon = True ,vmin=salt_min, vmax=salt_max, cmap='viridis')
    # plot colourbar
    plt.colorbar()
    # datetime title
    plt.title('Regional - Salinity (PSU)\n' + frame_time.strftime("%Y-%m-%d %H:%M:%S") + ' | ' + str(fname) + '_idx: ' + str(frame_idx))
    
    # make layout nice
    plt.tight_layout()
    
    # convert to array
    canvas = FigureCanvasAgg(fig)
    canvas.draw()
    frame = np.fromstring(canvas.tostring_rgb(), dtype='uint8')
    frame = frame.reshape(fig.canvas.get_width_height()[::-1] + (3,))
    
    return frame
