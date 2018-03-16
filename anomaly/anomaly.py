"""
Script to calculate mean anomalies for winter/summer for each year

STEPS:
1 - Extract all data into 1 dataframe
2 - Calculate means for each season/year pairs
3 - Calulate deviation from mean
4 - Make plots
"""

# animate_class_fast_new.py

from sklearn.linear_model import LogisticRegression
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap # basemap tools
import numpy as np
import pandas as pd
from datetime import datetime, timedelta # for working with datetimes
from netCDF4 import Dataset # reads netCDF file
from os import listdir
from os.path import isfile, join
import argparse
from scipy import stats

# __SETUP__

parser = argparse.ArgumentParser(description=__doc__)
# add a positional writable file argument
parser.add_argument('input_csv_file', type=argparse.FileType('r'))
# parser.add_argument('output_csv_file', type=argparse.FileType('w')) 
args = parser.parse_args()

# how many files
start = 0
end = 277

# for 
prob_ls = [] 
time_ls = []

# get probs
def get_prob(temp, salt, lr_model):
    # ravel to 1D array
    temp1d = temp.ravel()
    salt1d = salt.ravel()
    # make data frame and replace NaNs
    data = {'var1': temp1d, 'var2': salt1d}
    data = pd.DataFrame(data=data)
    data = data.fillna(-9999)
    # calculate probabilities
    probs = lr_model.predict_proba(data[['var1','var2']])
    prob_TSW, prob_EAC = zip(*probs)
    # convert tuples to list
    prob_EAC = list(prob_EAC)
    # sub back in nans
    prob_EAC = [x if x != 0.0 else np.nan for x in prob_EAC]
    # make 1D array
    prob_EAC = np.asarray(prob_EAC)
    # make 2D array
    prob_EAC = np.reshape(prob_EAC, (-1, 165))

    return prob_EAC

# add to list witout append
def add(lst, obj, index): return lst[:index] + [obj] + lst[index:]

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

def make_plot(data,lons,lats,title):
    fig = plt.figure()
    fig.subplots_adjust(left=0., right=1., bottom=0., top=0.9)
    # draw stuff
    m.drawcoastlines()
    m.fillcontinents(color='black')
    # plot color
    m.pcolor(lons,lats,np.squeeze(data), latlon = True ,vmin=-0.2, vmax=0.2, cmap='bwr')
    # datetime title
    plt.title(title)

    return fig

# __Make_Model__

# Read in classification training data
csv_data = pd.read_csv(args.input_csv_file, parse_dates = ['datetime'], 
                        infer_datetime_format = True) #Read as DateTime obsject

# make training dataset model
# create data points for training algorithm
# 1 = classA (EAC), 0 = classB (TS)
var1 = list(csv_data['temp'])
var2 = list(csv_data['salt'])
water_class = list(csv_data['class'])
# make data frame
train_data = {'var1': var1, 'var2': var2, 'class': water_class}
train_data = pd.DataFrame(data=train_data)
# replace current data strings with binary integers
train_data['class'] = train_data['class'].replace(to_replace='EAC', value=1)
train_data['class'] = train_data['class'].replace(to_replace='BS', value=0)
# fit logistic regression to the training data
lr_model = LogisticRegression()
lr_model = lr_model.fit(train_data[['var1','var2']], np.ravel(train_data[['class']]))

# __Grab_Data__

# get list of files in data directory
in_directory = "/Users/lachlanphillips/PhD_Large_Data/ROMS/Montague_subset"
file_ls = [f for f in listdir(in_directory) if isfile(join(in_directory, f))]
file_ls = list(filter(lambda x:'naroom_avg' in x, file_ls))
file_ls = sorted(file_ls)

# set output directory

# load one file to get lats and lons and data for initial plot
nc_file = in_directory + '/' + file_ls[1]
fh = Dataset(nc_file, mode='r')
lats = fh.variables['lat_rho'][:] 
lons = fh.variables['lon_rho'][:]

# get all the probability arrays
# itterate through all files
idx = 0
for i in range(start, end):
    # import file
    nc_file = in_directory + '/' + file_ls[i]
    fh = Dataset(nc_file, mode='r')
    print('Grabbing data from: '+file_ls[i]+' | '+str(i+1).zfill(3)+' of '+ str(len(file_ls)))
    fname = str(file_ls[i])[11:16]

    # extract time
    time = fh.variables['ocean_time'][:]

    # iterate through all time steps
    for j in range(0, len(time)):
        frame_time = grab_sst_time(int(j))
        # get data
        temp = fh.variables['temp'][j,29,:,:] 
        salt = fh.variables['salt'][j,29,:,:]
        # get probs 
        prob_EAC = get_prob(temp, salt, lr_model)
        # add data to lists
        prob_ls = add(prob_ls, prob_EAC, idx)
        time_ls = add(time_ls, frame_time, idx)
        idx += 1

    # close file
    fh.close()

# make time_ls into dataframe and calculate months
time_df = {'datetime':time_ls, 'month':[x.month for x in time_ls]}
time_df = pd.DataFrame(data=time_df)
summer = [1,2,3,10,11,12]
#time_df['season'] = ['summer' if x['month'] in summer else 'winter' for x in time_df]
time_df['season'] = time_df.month.apply(lambda x: 'summer' if x in summer else 'winter')
# index lists
summer_ls = time_df.index[time_df.season == 'summer']
winter_ls = time_df.index[time_df.season == 'winter']

# __Calculate_Anomalies__

# stack 2d arrays to one 3d array
data = np.dstack(prob_ls)
# calculate mean array
total_mean = np.mean(data, axis=2)
summer_mean = np.mean(data[:,:,summer_ls], axis=2)
winter_mean = np.mean(data[:,:,winter_ls], axis=2)

summer_anom = summer_mean - total_mean
winter_anom = winter_mean - total_mean

# make plots
# Setup map
##############################################################################
m = Basemap(projection='merc', llcrnrlat=-38.050653, urcrnrlat=-34.453367,\
        llcrnrlon=147.996456, urcrnrlon=152.457344, lat_ts=20, resolution='h')
##############################################################################

fig1 = make_plot(summer_anom, lons, lats, 'summer')
fig2 = make_plot(winter_anom, lons, lats, 'winter')
fig2 = make_plot(total_mean, lons, lats, 'mean')

plt.show()


"""
The plots don't give any useful information... need to rethink this... 
"""
















