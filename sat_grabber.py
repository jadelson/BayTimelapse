#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 29 22:36:45 2017

@author: jadelson
"""

import pickle
import numpy as np
import utm
from get_coord import get_coord
import netCDF4 as nc
#stress keys = ['k_x', 'node_northing', 'T', 'Cd', 'tau_w_y', 'tau_m_x', 'theta', 'tau_m_y', 'omega', 'easting', 'face_nodes', 'tau', 'k', 'U10', 'fw', 'v_vel', 'Rew', 'a0', 'depth', 'node_easting', 'u_vel', 'northing', 'k_y', 'tau_y', 'tau_w_x', 'A', 'surface']
# Getting back the objects:

def write_sat_data(satnc, rusty_data, date):
#sat_file = '/Volumes/Stella/landsat_order/landsat8_convert/LC80440332014189LGN00_L2.nc'

    #satnc = nc.Dataset(sat_file)
    sat_data = {}
    sat_keys = [ 'LAT', 'LON',  'RHOW_443', 'RHOW_483', 'RHOW_561', 'RHOW_655', 'RHOW_865']
    final_data = {}
    for k in sat_keys:
        final_data[k] = []
    for k in rusty_data.keys():
        final_data[k] = []
        
    for s in sat_keys:
        sat_data[s] = np.asarray(satnc[s][:,:])
        
    if np.min(sat_data['RHOW_655']>= 9.96921e35):
        return None
    
    maxlat = np.max(sat_data['LAT'])
    minlat = np.min(sat_data['LAT'])
    maxlon = np.max(sat_data['LON'])
    minlon = np.min(sat_data['LON'])
    m = rusty_data['lat'].shape
    n = sat_data['LAT'].shape
    A = None
    lat = rusty_data['lat']
    lon = rusty_data['lon']
    for i in range(0,m[0]):
        if np.mod(i,1000)==0:
            print(i*1./m[0])

        if lat[i] > maxlat or lat[i] < minlat or lon[i] > maxlon or lon[i] < minlon:
            continue
        x,y, A = get_coord(sat_data,lat[i],lon[i],n[0],n[1],A)
        if x >= n[0] or y >= n[1]:
            continue
#        print(sat_data['RHOW_655'][x,y])
        if sat_data['RHOW_655'][x,y] >= 9.96921e35 or sat_data['RHOW_655'][x,y] == np.nan:
            continue
        
        linedata = {k:sat_data[k][x,y] for k in sat_keys}
        for k in rusty_data.keys():
                linedata[k] = rusty_data[k][i]
    
        for k in linedata.keys():
            final_data[k].append(linedata[k])
    
    
    
    for k in final_data.keys():
        final_data[k] = np.asarray(final_data[k])
        
    with open('/Volumes/Stella/sat_and_rusty/sat_and_rusty_data_'+date+'.dk','wb') as f:
        pickle.dump(final_data,f)

    return final_data
#    









