#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 29 19:05:44 2017

"""

# Import libraries
from datetime import datetime
import netCDF4 as nc
import os
#from joblib import Parallel, delayed
import multiprocessing
import numpy as np
import pickle
from get_coord import get_coord
import build_dataset
from bay_remote_sensing_init import *
import itertools

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
    for k in sat_keys:
        final_data[k] = []
    for k in tide_data.keys():
        final_data[k] = []
        
    for s in sat_keys:
        sat_data[s] = np.asarray(satnc[s][:,:])

    if np.min(sat_data[sat_test_key]>= 9.96921e35):
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
#            print('fail1')
            continue
        x,y, A = get_coord(sat_data,lat[i],lon[i],n[0],n[1],A)
        if x >= n[0] or y >= n[1]:
#            print('fail2')
            continue
        if sat_data[sat_test_key][x,y] >= 9.96921e35 or np.isnan(sat_data[sat_test_key][x,y]):
#            print('fail3')
            continue
        
        linedata = {k:sat_data[k][x,y] for k in sat_keys}
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
    
    """
    Extract satelite data onto the tidal model grid
    
    :param filename: name of the Acolite corrected netcdf satellite file
    :param tide_data: dict of tidal data at each point in domain
    """
    print('processing ' + filename + ' on CPU: ', str(multiprocessing.current_process()))
    
    model = nc.Dataset(raw_sat_directory+filename)
    if not tide_data == None:
        combined_sat_data = write_sat_data(model, tide_data, model.DATE+'_'+model.HOUR+model.MINUTE+model.SECOND[0:2])
        if not combined_sat_data == None:
            print('Sat: ' + raw_sat_directory+filename)
            print(len(tide_data['lat']), len(combined_sat_data['lat']))
    model.close()
    return
 


if __name__ == "__main__":
    print('Satellite data filteration running...')
    # Load tide model data (requires current velocities, water level, and depth)
    blank_date = datetime.strptime('20160301 18 30 00','%Y%m%d %H %M %S')        
    tide_data = build_dataset.find_tide_data_by_time(blank_date, base_dir + 'raw_rusty_tides/')
    sat_inputs = [k for k in os.listdir(raw_sat_directory) if k.endswith('.nc')]
#    blah = sat_worker(sat_inputs[3], tide_data)

    jobs = []
    with multiprocessing.Pool() as p:
        p.starmap(sat_worker, zip(sat_inputs, itertools.repeat(tide_data)))
        p.close()
        p.join()

    print('Satellite Data filteration Complete.')

##  
    print('Satellite data combination running...')             
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
            full_data['date'] = np.append(full_data['date'],[datetime.strptime(filename[29:44],'%Y%m%d_%H%M%S')]*len(combined_sat_data[k]))
            f.close()
##  
    print('Satellite data combination complete.')             
    with open(working_dir+sat+'_no_stress_dataset.dk','wb') as f:
        pickle.dump(full_data, f)
    print('Data dump complete.')
##        