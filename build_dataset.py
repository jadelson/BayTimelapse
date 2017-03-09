#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 29 19:05:44 2017

@author: jadelson
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

# Set directory locations
sat_directory = '/Volumes/Stella/landsat_order/landsat8_convert/'
sat_tide_directory = '/Volumes/Stella/sat_and_rusty/'
wind_sat_tide_directory = '/Volumes/Stella/wind_sat_and_rusty/'
stress_sat_directory = '/Volumes/Stella/stress_sat/'
wind_data_file = '/Volumes/Stella/weather_data/weather_data_2014-2017.dk'
fetch_file = 'cached_fetched.dk'

## Here are some functions!

# A function that writes a combined satellite and tidal dataset
def write_sat_data(satnc, tide_data, date):
    """
    Combine satellite and tide data then write to file
    
    :param satnc: Acolite corrected satellite remote sensing netcdf file
    :param tide_data: dict of tidal data at each point in domain
    :param date: date of data
    """
    sat_data = {}
    sat_keys = [ 'LAT', 'LON',  'RHOW_443', 'RHOW_483', 'RHOW_561', 'RHOW_655', 'RHOW_865'] #Landsat 8
    final_data = {}
    for k in sat_keys:
        final_data[k] = []
    for k in tide_data.keys():
        final_data[k] = []
        
    for s in sat_keys:
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
        if sat_data['RHOW_655'][x,y] >= 9.96921e35 or np.isnan(sat_data['RHOW_655'][x,y]):
            continue
        
        linedata = {k:sat_data[k][x,y] for k in sat_keys}
        for k in tide_data.keys():
                linedata[k] = tide_data[k][i]
    
        for k in linedata.keys():
            final_data[k].append(linedata[k])

    for k in final_data.keys():
        final_data[k] = np.asarray(final_data[k])
        
    with open(sat_tide_directory+'sat_and_rusty_data_'+date+'.dk','wb') as f:
        pickle.dump(final_data,f)
        print('Saved: ' + sat_tide_directory+'sat_and_rusty_data_'+date+'.dk')
    return final_data



# This function extracts tidal data for a given time
def find_tide_data_by_time(test_time): 
    path = '/Volumes/Stella/raw_rusty_tides/' 
    base_time = datetime(2014,9,1,8,0) 
    ad_1 = datetime(1,1,1,8,0)

    
    days_since = test_time - ad_1
    file_access = days_since.days + 1
    
    tide_time_since = test_time - base_time
    seconds_since = tide_time_since.total_seconds()

    if not  os.path.isfile(path + 'data' + str(file_access) + '.nc'):
        print('Cannot find file: data' + str(file_access) + '.nc')
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
    tide_data = {'lat': lat.flatten(),'lon':lon.flatten(), 'u':u.flatten(), 'v':v.flatten(), 'depth':depth.flatten(), 'eta':eta}

    return tide_data


# This is a simple wind interpolation sceme for SF Bay that breaks the bay into 
# four seperate regions each with identical wind conditions based on the
# nearest NOAA station
def simple_sfbay_windterp(x_stations, y_stations, u_stations, v_stations, \
                          x_grid, y_grid, depth_grid):
    
    U = []
    for uu in u_stations:
        U.append(uu.ravel()[0])
    U = np.asarray(U)   
    
    V = []
    for vv in v_stations:
        V.append(vv.ravel()[0])
    V = np.asarray(V)   
   
#    UU = U[~np.isnan(U)] 
#    VV = V[~np.isnan(V)]
    UU = U
    VV = V
    
    X = np.asarray(x_stations)
    Y = np.asarray(y_stations)

    lat1 = 4.185e6
    lat2 = 4.203e6
    lon4 = 5.8e5
    lon5 = 5.42e5 #ocean removal
    lat5 = 4.1879e6 #ocean removal
    lon6 = 5.49513e5 #ocean removal
    lat6 = 4.14842e6 #ocean removal
    lon7 = 6.05e5 #delta removal
    
    lon_mont_min = 576030.82726326736 #Montezuma slough removal
    lon_mont_max = 583010.8610822059
    lat_mont_min = 4220561.5234375
    lat_mont_max = 4246293.9453125

    x2 = np.ravel(x_grid)
    y2 = np.ravel(y_grid)
    
    u2 = np.zeros(x2.shape)
    v2 = np.zeros(x2.shape)
    #wind locations
    full_ind = [None]*4
    full_ind[0] = (y2 < lat1) & (x2 >= lon5)
    full_ind[1] = (y2 >= lat1) & (y2 < lat2) & (x2 < lon4) & (x2 >= lon5)
    full_ind[2] = (y2 >= lat2) & (x2 < lon4) & (x2 >= lon5)
    full_ind[3] = (y2 >= lat1) & (x2 >= lon4) & (x2 >= lon5)
    
    base_ind = [None]*4
    base_ind[0] = (Y < lat1)
    base_ind[1] = (Y >= lat1) & (Y < lat2) & (X < lon4)

    base_ind[2] = (Y >= lat2) & (X < lon4) 

    base_ind[3] = (Y >= lat1) & (X >= lon4) 
    
    
    for i in range(4):
        u2[full_ind[i]] = UU[base_ind[i]]
        v2[full_ind[i]] = VV[base_ind[i]]

    umag = np.hypot(u2, v2)

    indx = (depth_grid < 50)
    indx = indx & (~np.isnan(umag))
#    indx = indx &(umag<20 )
    indx = indx & ((lat6-lat5)*(x2 - lon5) <= (lon6-lon5)*(y2 - lat5))
    indx = indx & (x2 < lon7) 
    indx = indx & ((x2 > lon_mont_min) & (x2 < lon_mont_max) & (y2 > lat_mont_min) & (y2 < lat_mont_max))
    
    u_wind = u2
    v_wind = v2
    umag_wind = umag
    valid_location_indx = indx
    return u_wind, v_wind, umag_wind, valid_location_indx

# Calculate the stress given a dataset of depth, surface level, wave period,
# wave height, and tidal currents.
def calculate_stress(data):
    # Physical constants
    kappa = 0.41
    rho = 1023
    nu = 10**-6
    z0 = 10**-4
    g = 9.81
     
    m = data['a0'].shape  
    data['H'] = np.empty(m)
    data['k'] = np.empty(m)
    data['A'] = np.empty(m)
    data['k_x'] = np.empty(m)
    data['k_y'] = np.empty(m)
    data['fw'] = np.empty(m)
    data['omega'] = np.empty(m)
    data['Rew'] = np.empty(m)
    data['Cd'] = np.empty(m)
    data['tau_m_x'] = np.empty(m)
    data['tau_w_x'] = np.empty(m)
    data['tau_m_y'] = np.empty(m)
    data['tau_w_y'] = np.empty(m)
    data['tau_x'] = np.empty(m)
    data['tau_y'] = np.empty(m)
    data['tau_w'] = np.empty(m)
    data['tau_m'] = np.empty(m)

    data['tau'] = np.empty(m)
    
    data['omega'] = 2*np.pi/data['T']

    depth = np.array(data['depth'])
    surface = np.array(data['eta'])
    u_vel_1 = np.array(data['u'])
    v_vel_1 = np.array(data['v'])

    for i in range(0,m[0]):

        
        H = depth[i] + surface[i]
        Cd = kappa**2*(np.log(H/z0) + z0/H - 1)**-2
        u = u_vel_1[i]
        v = v_vel_1[i]
        mag_u = np.sqrt(u**2 + v**2)
        tau_m_x = rho*Cd*mag_u*u
        tau_m_y = rho*Cd*mag_u*v
            
        omega = data['omega'][i]
        
        a0 = data['a0'][i]
        angle = data['theta'][i]   
        

        def dispersion_relation(k):
            return np.abs(omega**2 - g*k*np.tanh(k*H))
        k = minimize(dispersion_relation,omega/(np.sqrt(g*H)),
                options={'disp': False})
        mag_k = np.abs(k.x)
        
        kx = mag_k*np.cos(angle)
        ky = -mag_k*np.sin(angle)
    
        A = a0/(np.sinh(mag_k*H))
        Rew = A**2*omega/nu
        fw = 2*Rew**(-0.5)
        mag_tau_w = 0.5*rho*fw*(A*omega)**2
        tau_w_x = mag_tau_w*kx/mag_k
        tau_w_y = mag_tau_w*ky/mag_k
    
        tau_x = tau_w_x + tau_m_x
        tau_y = tau_w_y + tau_m_y
        tau = np.sqrt(tau_x**2 + tau_y**2)    
 
        
        data['H'][i] = H
        data['k'][i] = mag_k
        data['k_x'][i] = kx
        data['k_y'][i] = ky
        data['A'][i] = A
        data['fw'][i] = fw      
        data['Rew'][i] = Rew
        data['Cd'][i] = Cd 
        data['tau_m_x'][i] = tau_m_x
        data['tau_w_x'][i] = tau_w_x
        data['tau_m_y'][i] = tau_m_y
        data['tau_w_y'][i] = tau_w_y
        data['tau_w'][i] = mag_tau_w
        data['tau_m'][i] = np.sqrt(tau_m_x*tau_m_x + tau_m_y*tau_m_y)
        data['tau_x'][i] = tau_x
        data['tau_y'][i] = tau_y
        data['tau'][i] = tau
    
    return data



# Calculate stress from wave and tide models at each point in the tidal model
def stress_worker(filename, fetch_model):
    if not filename.endswith(".dk"): 
        return 'FAIL ' + filename
    print(filename)
    # Extract date
    date_string = filename[19:27]
    file_date = datetime.strptime(date_string,'%Y%m%d')
    satf = [s for s in os.listdir(sat_directory) if file_date.strftime('%Y%j') in s]
    satf = satf[0]
    model = nc.Dataset(sat_directory+satf)
    sat_datetime = datetime.strptime(date_string + ' ' + model.TIME[0:8],'%Y%m%d %H %M %S')
    save_date = sat_datetime.strftime('%Y%m%d_%H%M%S')
    epoch = datetime.utcfromtimestamp(0)
    time = (sat_datetime - epoch).total_seconds()/60
    with open(sat_tide_directory+filename, 'rb') as f:
        sat_tide_data = pickle.load(f)

    x2 = np.zeros(sat_tide_data['lat'].shape)
    y2 = np.zeros(sat_tide_data['lat'].shape)
    print(len(x2))
    for i in range(len(x2)): 
            
        x_y = utm.from_latlon(sat_tide_data['lat'][i],sat_tide_data['lon'][i])
        x2[i] = x_y[0]
        y2[i] = x_y[1]
    
    winds = fetch_model.data
    
    x = []
    y = []
    u = []
    v = []

    for w,k in zip(winds.values(),winds.keys()):
        if not (k == '994034-99999' or k == '724940-23234'or k == '994036-99999 ' 
                or k == '998477-99999' or k == '724955-93227'):#'720646-99999' or k == '994016-99999':
            continue       
        temp = w['Uwindf'](time)
        u.append(temp)
        temp = w['Vwindf'](time)
        v.append(temp)
        x.append(w['easting'])
        y.append(w['northing']) 
    
    u_wind, v_wind, umag_wind, indx = simple_sfbay_windterp(x,y,u,v,x2,y2,sat_tide_data['depth'])
    wind_sat_tide_data = sat_tide_data
    wind_sat_tide_data['u_wind'] = u_wind
    wind_sat_tide_data['v_wind'] = v_wind
    wind_sat_tide_data['U10'] = umag_wind

    
    for k in wind_sat_tide_data.keys():
        wind_sat_tide_data[k] = wind_sat_tide_data[k][indx]
    theta = np.arctan2(u_wind[indx], v_wind[indx]) - np.pi
                      
    u10_0 = umag_wind[indx]
    x_0 = x2[indx]
    y_0 = y2[indx]
    T = np.zeros(x_0.shape)
    Hs = np.zeros(x_0.shape)
    for i in range(len(x_0)):
        Ttemp, Hstemp = fetch_model.getwaveheight(np.array(x_0[i]), np.array(y_0[i]), np.array(theta[i]), np.array(u10_0[i]))
        T[i] = Ttemp[0]
        Hs[i] = Hstemp[0]
    wind_sat_tide_data['theta'] = umag_wind


    wind_sat_tide_data['easting'] = x_0
    wind_sat_tide_data['northin'] = y_0
    wind_sat_tide_data['T'] = T
    wind_sat_tide_data['a0'] = Hs
    indx2 = ~ np.isnan(Hs)
    
    for k in wind_sat_tide_data.keys():
        wind_sat_tide_data[k] = wind_sat_tide_data[k][indx2]
        
#    # Save raw wind data without calculating stress for ease of future use
#    with open('/Volumes/Stella/wind_sat_and_rusty/wind_sat_rusty_data_'+save_date+'.dk','wb') as g:
#        pickle.dump(wind_sat_tide_data,g)
#        
    # Calculate the sat-stress database from tide and wave data then save
    stress_sat_data = calculate_stress(wind_sat_tide_data)
    with open('/Volumes/Stella/stress_sat/stress_sat_data_'+save_date+'.dk','wb') as g:
        pickle.dump(stress_sat_data,g)     

    return filename      


def sat_worker(filename):
    model = nc.Dataset(sat_directory+filename)
    sat_date = datetime.strptime(model.DATE + ' '+ model.TIME[0:8],'%Y%m%d %H %M %S')        
    tide_data = find_tide_data_by_time(sat_date)
    if not tide_data == None:
        sat_tide_data = write_sat_data(model, tide_data, model.DATE+model.HOUR+model.MINUTE+':'+model.SECOND[0:2])
        if not sat_tide_data == None:
            print('Sat: ' + sat_directory+filename)
            print(len(tide_data['lat']), len(sat_tide_data['lat']))
    model.close()
    return(sat_date)
    
    
if __name__ == "__main__":
    # Load tide model data (requires current velocities, water level, and depth)
#    sat_inputs = os.listdir(sat_directory)
#    sat_inputs = [k for k in sat_inputs if k.endswith('.nc')]
#    results = Parallel(n_jobs=3)(delayed(sat_worker)(filename) for filename in sat_inputs)
   
    # Load results from the wave model (requires wave height, wave direction, wave period)
    fetch_model = FetchModel()
    if not os.path.isfile(wind_data_file): 
        wind_station_data = fetch_model.downloadwinds(2014,2017)
        fetch_model.savedata(wind_station_data,wind_data_file)
    else:
        wind_station_data = fetch_model.loadwindfromfile(wind_data_file)
#     
#              
    # we run the stress_worker fucntion for each satellite image we have 
    # downloaded that writes a full stress file dataset for each filename
    inputs = os.listdir(sat_tide_directory)
    results = Parallel(n_jobs=3)(delayed(stress_worker)(filename, fetch_model) for filename in inputs)
#    for i in inputs:
#        fetches = stress_worker(i, fetch_model)
#    stress_worker(filename, fetch_model)
##                     
#    full_data = {}
#    full_data['date'] = np.array([])              
#    for filename in os.listdir(stress_sat_directory):
#        stress_sat_data = {}
#        print(stress_sat_directory+filename)
#        with open(stress_sat_directory+filename,'rb') as f:
#            stress_sat_data = pickle.load(f)
#            indx = ~ (np.isnan(stress_sat_data['tau']) | np.isinf(stress_sat_data['tau']))
#            for k in stress_sat_data.keys():
#                if not k in full_data.keys():
#                    full_data[k] = np.array([])
#                full_data[k] = np.append(full_data[k], stress_sat_data[k][indx])
#            full_data['date'] = np.append(full_data['date'],[datetime.strptime(filename[16:31],'%Y%m%d_%H%M%S')]*sum(indx))
#            f.close()
            
#    with open('/Users/jadelson/Dropbox/phdResearch/AllOptical/timelapse_work/full_dataset.dk','wb') as f:
#        pickle.dump(full_data, f)
    #    