#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 20 12:53:18 2016

@author: jadelson
"""
import numpy as np

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
    
    
    UU = U[~np.isnan(U)]
    
    VV = V[~np.isnan(V)]
    UU = U
    VV = V
    
    X = np.asarray(x_stations)
    Y = np.asarray(y_stations)
#    X = X[~np.isnan(U)]
#    Y = Y[~np.isnan(U)]
#    print(X.shape,Y.shape)
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
    #D = np.column_stack((X,Y,UU,VV))
    #
    #
    #from scipy.interpolate import interp2d
    #
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
    #u2 = np.asarray(u2).ravel()
    #v2 = np.asarray(v2).ravel()
    umag = np.hypot(u2, v2)
    print(umag)
    umag_temp = umag
    umag_temp[np.isnan(umag)] = 100000
    indx = (depth_grid < 50)
    indx = indx & (~np.isnan(umag))
    indx = indx &(umag<20 )
    indx = indx & ((lat6-lat5)*(x2 - lon5) <= (lon6-lon5)*(y2 - lat5))
    indx = indx & (x2 < lon7) 
    indx = indx & ~((x2 > lon_mont_min) & (x2 < lon_mont_max) & (y2 > lat_mont_min) & (y2 < lat_mont_max))
    
    u_wind = u2
    v_wind = v2
    umag_wind = umag
    valid_location_indx = indx
    return u_wind, v_wind, umag_wind, valid_location_indx


if __name__ == "__main__":
     u_wind, v_wind, umag_wind, indx = simple_sfbay_windterp(x,y,u,v,x2,y2,sat_rusty_data['depth'])

