#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 29 19:05:44 2017

"""

# Import libraries
from datetime import datetime
import netCDF4 as nc
import os
from calculate_fetch import FetchModel
import utm
from joblib import Parallel, delayed
import numpy as np
import pickle
from get_coord import get_coord
from scipy.optimize import minimize
import build_dataset

# Set directory locations
base_dir = '/Volumes/Stella/'

raw_sat_directory = base_dir + 'landsat_order/landsat8_convert/'
sat_filter_directory = base_dir + 'sat_on_mesh/'

sat = 'l8'
if sat == 'l8':
    l8_keys = [ 'LAT', 'LON',  'RHOW_443', 'RHOW_483', 'RHOW_561', 'RHOW_655', 'RHOW_865'] #Landsat 8
    l8_test_key = 'RHOW_655'



# A function that writes a combined satellite and tidal dataset
def write_sat_data(satnc, tide_data, date):
    """
    Combine satellite and tide data then write to file
    
    :param satnc: Acolite corrected satellite remote sensing netcdf file
    :param tide_data: dict of tidal data at each point in domain
    :param date: date of data
    """
    sat_data = {}
    final_data = {}
    for k in l8_keys:
        final_data[k] = []
    for k in tide_data.keys():
        final_data[k] = []
        
    for s in l8_keys:
        sat_data[s] = np.asarray(satnc[s][:,:])

    if np.min(sat_data['RHOW_655']>= 9.96921e35):
        return None
    
    maxlat = np.max(sat_data['LAT'])
    minlat = np.min(sat_data['LAT'])
    maxlon = np.max(sat_data['LON'])
    minlon = np.min(sat_data['LON'])
    m = tide_data['lat'].shape
    n = sat_data['LAT'].shape
    A = None
    lat = tide_data['lat']
    lon = tide_data['lon']
    
    
    for i in range(0,m[0]):

        if lat[i] > maxlat or lat[i] < minlat or lon[i] > maxlon or lon[i] < minlon:
            continue
        x,y, A = get_coord(sat_data,lat[i],lon[i],n[0],n[1],A)
        if x >= n[0] or y >= n[1]:
            continue
        if sat_data['RHOW_655'][x,y] >= 9.96921e35 or np.isnan(sat_data[l8_test_key][x,y]):
            continue
        
        linedata = {k:sat_data[k][x,y] for k in l8_keys}
        for k in tide_data.keys():
                linedata[k] = tide_data[k][i]
    
        for k in linedata.keys():
            final_data[k].append(linedata[k])

    for k in final_data.keys():
        final_data[k] = np.asarray(final_data[k])
    
    this_filename = sat_filter_directory+sat+'_timelapse_sat_on_grid_tide'+date+'.dk'
    with open(this_filename,'wb') as f:
        pickle.dump(final_data,f)
        print(this_filename)
    return final_data

   


def sat_worker(filename, tide_data):
    model = nc.Dataset(sat_directory+filename)
    if not tide_data == None:
        sat_tide_data = write_sat_data(model, tide_data, model.DATE+'_'+model.HOUR+model.MINUTE+model.SECOND[0:2])
        if not sat_tide_data == None:
            print('Sat: ' + sat_directory+filename)
            print(len(tide_data['lat']), len(sat_tide_data['lat']))
    model.close()
    return    
    
if __name__ == "__main__":
    # Load tide model data (requires current velocities, water level, and depth)
    blank_date = datetime.strptime('20160301 18 30 00','%Y%m%d %H %M %S')        
    tide_data = build_dataset.find_tide_data_by_time(blank_date)
    sat_inputs = [k for k in os.listdir(sat_directory) if k.endswith('.nc')]
    sat_worker(sat_inputs[1], tide_data)
#    results = Parallel(n_jobs=3)(delayed(sat_worker)(filename) for filename in sat_inputs)
#   

##                     
    full_data = {}
    full_data['date'] = np.array([])              
    for filename in os.listdir(sat_filter_directory):
        if not filename[0:2] == sat:
            continue
        combined_sat_data = {}
        print(sat_filter_directory+filename)
        with open(sat_filter_directory+filename,'rb') as f:
            combined_sat_data = pickle.load(f)
            for k in combined_sat_data.keys():
                if not k in full_data.keys():
                    full_data[k] = np.array([])
                full_data[k] = np.append(full_data[k], combined_sat_data[k])
            full_data['date'] = np.append(full_data['date'],[datetime.strptime(filename[29:44],'%Y%m%d_%H%M%S')]*combined_sat_data[k])
            f.close()
##            
#    with open('/Users/jadelson/Dropbox/phdResearch/AllOptical/timelapse_work/full_dataset.dk','wb') as f:
#        pickle.dump(full_data, f)
##        