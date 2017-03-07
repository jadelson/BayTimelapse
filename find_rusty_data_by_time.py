#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 29 17:23:56 2017

@author: jadelson
"""

import netCDF4 as nc
import numpy as np
import utm
from datetime import datetime
import os

def find_rusty_data_by_time(test_time): 
    path = '/Volumes/Stella/raw_rusty_tides/' 
    #files2 =  os.listdir(path2)
    base_time = datetime(2014,9,1,8,0) 
    ad_1 = datetime(1,1,1,8,0)

    
    days_since = test_time - ad_1
    file_access = days_since.days + 1
    
    rusty_time_since = test_time - base_time
    seconds_since = rusty_time_since.total_seconds()

    if not  os.path.isfile(path + 'data' + str(file_access) + '.nc'):
        return None
    
    model = nc.Dataset(path + 'data' + str(file_access) + '.nc')
    
    timeindex = np.argmin(np.abs(model['time'][:] - seconds_since))

    node_easting = model['global_node_x'][:]
    node_northing = model['global_node_y'][:]
    depth = model['depth'][:]
    eta = model['surface'][:,timeindex]
    u = model['cell_east_velocity'][:,:,timeindex]
    v = model['cell_north_velocity'][:,:,timeindex]
    face_nodes = model['face_nodes'][:]
    
     
    m = depth.shape
    lat = np.empty((m))
    lon = np.empty((m))
    
    
    for i in range(0,m[0]):
    	northingsum = 0
    	eastingsum = 0
    	for f in face_nodes[i,:]:
    		northingsum = northingsum + node_northing[f]
    		eastingsum = eastingsum + node_easting[f]
    	north = northingsum/3
    	east = eastingsum/3
    		
    	templat, templon = utm.to_latlon(east,north,10,'S')
    	lat[i]= templat
    	lon[i] = templon
    
    u = u.transpose()
    v = v.transpose()
    rusty_data = {'lat': lat.flatten(),'lon':lon.flatten(), 'u':u.flatten(), 'v':v.flatten(), 'depth':depth.flatten(), 'eta':eta}

    return rusty_data
#scipy.io.savemat('rusty.mat', mdict={'lat': lat,'lon':lon, 'u':u[:], 'v':v[:], 'depth':depth, 'eta':eta})



