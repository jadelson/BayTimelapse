#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 13:09:05 2017

@author: jadelson
"""
sat = 'l8'
base_dir = '/Volumes/Stella/'
working_dir = './'

if sat == 'l8':
    sat_keys = [ 'LAT', 'LON',  'RHOW_443', 'RHOW_483', 'RHOW_561', 'RHOW_655', 'RHOW_865'] #Landsat 8
    sat_test_key = 'RHOW_655'
    raw_sat_directory = base_dir + 'landsat_order/landsat8_convert/'
    wind_data_file = base_dir + 'weather_data/weather_data_2014-2017.dk'
if sat == 'l7':
    sat_keys = [ 'LAT', 'LON',  'RHOW_479', 'RHOW_561', 'RHOW_661', 'RHOW_835'] #Landsat 7
    sat_test_key = 'RHOW_661'
    raw_sat_directory = base_dir + 'landsat_order/landsat7_convert/'
    wind_data_file = base_dir + 'weather_data/weather_data_2000-2017.dk'


sat_tide_directory = base_dir + 'sat_and_rusty/'
wind_sat_tide_directory = base_dir + 'wind_sat_and_rusty/'
stress_sat_directory = base_dir + 'stress_sat/'
fetch_file = sat + '_savedfetchmodel.dk'
sat_filter_directory = base_dir + 'sat_on_mesh/'




