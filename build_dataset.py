#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 29 19:05:44 2017

@author: jadelson
"""

# Import libraries
from find_rusty_data_by_time import *
from datetime import datetime, timedelta
import netCDF4 as nc
import os
from sat_grabber import *
from simple_sfbay_windterp import simple_sfbay_windterp
from SFBayRemoteSensing.conversions.calculate_fetch import FetchModel
import utm
from calculate_stress3 import calculate_stress
from joblib import Parallel, delayed
import numpy as np
import pickle

# Set directory locations
sat_directory = '/Volumes/Stella/landsat_order/landsat8_convert/'
sat_rusty_directory = '/Volumes/Stella/sat_and_rusty/'
wind_sat_rusty_directory = '/Volumes/Stella/wind_sat_and_rusty/'
stress_sat_directory = '/Volumes/Stella/stress_sat/'

# Load tide model data (requires current velocities, water level, and depth)
for filename in os.listdir(sat_directory):
    if filename.endswith(".nc"): 
        model = nc.Dataset(sat_directory+filename)
        sat_date = datetime.strptime(model.DATE + ' '+ model.TIME[0:8],'%Y%m%d %H %M %S')        
        rusty_data = find_rusty_data_by_time(sat_date)
        if not rusty_data == None:
            print(filename, model.DATE, model.TIME)
#            sat_rusty_data = write_sat_data(model, rusty_data, model.DATE+model.HOUR+model.MINUTE+':+model.SECOND[0:2])         
        model.close()
#
#fetch_model = FetchModel()
#    wind_station_data = fetch_model.downloadwinds(2014,2017)
#    fetch_model.savedata(wind_station_data,'/Volumes/Stella/weather_data/weather_data_2014-2017.dk')
#wind_station_data = fetch_model.loadwindfromfile('/Volumes/Stella/weather_data/weather_data_2014-2017.dk')
#    
#def stress_worker(filename):
#    if not filename.endswith(".dk"): 
#        return 'FAIL' + filename
#    date_string = filename[19:27]
#    file_date = datetime.strptime(date_string,'%Y%m%d')
#    satf = [s for s in os.listdir(sat_directory) if file_date.strftime('%Y%j') in s]
#    satf = satf[0]
#    model = nc.Dataset(sat_directory+satf)
#    sat_datetime = datetime.strptime(date_string + ' ' + model.TIME[0:8],'%Y%m%d %H %M %S')
#    save_date = sat_datetime.strftime('%Y%m%d_%H%M%S')
#    epoch = datetime.utcfromtimestamp(0)
#    time = (sat_datetime - epoch).total_seconds()/60
#    with open(sat_rusty_directory+filename, 'rb') as f:
#        sat_rusty_data = pickle.load(f)
#    winds = fetch_model.data
#    
#    
#    x = []
#    y = []
#    u = []
#    v = []
#    for w,k in zip(winds.values(),winds.keys()):
#        if not (k == '994034-99999' or k == '724940-23234'or k == '994036-99999 ' 
#                or k == '998477-99999' or k == '724955-93227'):#'720646-99999' or k == '994016-99999':
#            continue
#       
#        temp = w['Uwindf'](time)
#        u.append(temp)
#        temp = w['Vwindf'](time)
#        v.append(temp)
#        x.append(w['easting'])
#        y.append(w['northing'])
#    x2 = np.zeros(sat_rusty_data['lat'].shape)
#    y2 = np.zeros(sat_rusty_data['lat'].shape)
#    for i in range(len(x2)): 
#            
#        x_y = utm.from_latlon(sat_rusty_data['lat'][i],sat_rusty_data['lon'][i])
#        x2[i] = x_y[0]
#        y2[i] = x_y[1]
#    
#    
#    u_wind, v_wind, umag_wind, indx = simple_sfbay_windterp(x,y,u,v,x2,y2,sat_rusty_data['depth'])
#    wind_sat_rusty_data = sat_rusty_data
#    wind_sat_rusty_data['u_wind'] = u_wind
#    wind_sat_rusty_data['v_wind'] = v_wind
#    wind_sat_rusty_data['U10'] = umag_wind

#                   
##            # ADJUST THE INDX
##            indx = indx[0::100]
#    
#    for k in wind_sat_rusty_data.keys():
#        wind_sat_rusty_data[k] = wind_sat_rusty_data[k][indx]
#    theta = np.arctan2(u_wind[indx], v_wind[indx]) - np.pi
#                      
#    u10_0 = umag_wind[indx]
#    x_0 = x2[indx]
#    y_0 = y2[indx]
#    T = np.zeros(x_0.shape)
#    Hs = np.zeros(x_0.shape)
#    
#    wind_sat_rusty_data['theta'] = umag_wind
#    for i in range(0,len(x_0)):
#        T[i],Hs[i] = fetch_model.getwaveheight(x_0[i], y_0[i], theta[i], u10_0[i])
#    
##        if np.mod(i,400) == 0:
##            print(i/len(x_0)) 
#    wind_sat_rusty_data['easting'] = x_0
#    wind_sat_rusty_data['northin'] = y_0
#    wind_sat_rusty_data['T'] = T
#    wind_sat_rusty_data['a0'] = Hs
#    indx2 = ~ np.isnan(Hs)
#    
#    for k in wind_sat_rusty_data.keys():
#        wind_sat_rusty_data[k] = wind_sat_rusty_data[k][indx2]
#    with open('/Volumes/Stella/wind_sat_and_rusty/wind_sat_rusty_data_'+save_date+'.dk','wb') as g:
#        pickle.dump(wind_sat_rusty_data,g)
#
#    stress_sat_data = calculate_stress(wind_sat_rusty_data)
#    with open('/Volumes/Stella/stress_sat/stress_sat_data_'+save_date+'.dk','wb') as g:
#        pickle.dump(stress_sat_data,g)     
#
#    return filename      
#        
#inputs = os.listdir(sat_rusty_directory)
#results = Parallel(n_jobs=8)(delayed(stress_worker)(filename) for filename in inputs)


   
#                    
##                
full_data = {}
full_data['date'] = np.array([])              #%0%
for filename in os.listdir(stress_sat_directory):
    stress_sat_data = {}
    print(stress_sat_directory+filename)
    with open(stress_sat_directory+filename,'rb') as f:
        stress_sat_data = pickle.load(f)
        indx = ~ (np.isnan(stress_sat_data['tau']) | np.isinf(stress_sat_data['tau']))
        for k in stress_sat_data.keys():
            if not k in full_data.keys():
                full_data[k] = np.array([])
            full_data[k] = np.append(full_data[k], stress_sat_data[k][indx])
        full_data['date'] = np.append(full_data['date'],[datetime.strptime(filename[16:31],'%Y%m%d_%H%M%S')]*sum(indx))
        f.close()
        
#with open('/Users/jadelson/Dropbox/phdResearch/AllOptical/sfbayrsproper/full_dataset.dk','wb') as f:
#    pickle.dump(full_data, f)
#    