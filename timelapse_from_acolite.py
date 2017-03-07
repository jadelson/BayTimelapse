#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 13:20:34 2017

@author: jadelson
"""
import pandas as pd
from datetime import datetime, timedelta
import urllib
import json
import bisect

with open('/Users/jadelson/Dropbox/phdResearch/AllOptical/sfbayrsproper/full_dataset.dk','rb') as f:
    full_data = pickle.load(f)
data_length = len(full_data['tau'])   
date_list = []

for v in full_data['date']:
    if not v in date_list:
        date_list.append(v)

tide_df = pd.DataFrame()
begin_date = min(date_list)
end_date = begin_date + timedelta(days = 365)
last_date = max(date_list)
while True:
    begin_time_str = begin_date.strftime('%Y%m%d%%2000:00')
    end_time_str = end_date.strftime('%Y%m%d%%2023:59')
    url = ('https://tidesandcurrents.noaa.gov/api/datagetter?'
        'begin_date=%s'
        '&end_date=%s'
        '&station=9414290'
        '&product=high_low'
        '&datum=mllw'
        '&units=metric'
        '&time_zone=gmt'
        '&application=web_services'
        '&format=json')
    url = url % (begin_time_str,end_time_str)
    print(url)
    try:
        response = urllib.request.urlopen(url)
    except Exception as e :
        print(e)
    x = response.read().decode('utf-8')
    json_df = json.loads(x)
    if 'error' in json_df:
        print(json_df['error'])
        break
    else:
        tide_df_temp = json_df['data']
        tide_df_temp = pd.DataFrame.from_dict(tide_df_temp)
        tide_df_temp['t'] = pd.to_datetime(tide_df_temp['t'],infer_datetime_format=True)
        tide_df_temp = tide_df_temp.set_index('t')
        tide_df = tide_df.append(tide_df_temp)
    if end_date > last_date:
        break
    else:
        begin_date = end_date
        end_date = begin_date + timedelta(days = 365) 
        
        
low_water_times = tide_df.index[tide_df['ty']=='LL']
ll_df = tide_df[tide_df['ty']=='LL']
ll_length = len(ll_df.index)
phases = {}
j = 0
for date in date_list:
    
    i = bisect.bisect_left(ll_df.index,date)
    if i >= ll_length:
        phases[date] = np.nan
        j += 1
        continue
    #np.argmin(np.abs(ll_df.index.to_pydatetime() - index)+(((ll_df.index.to_pydatetime() - index) > timedelta(0))*timedelta(days=3)))
    td =  ll_df.index[i] - date
    phases[date] = td.total_seconds()
    j += 1
    

    
data_by_phase = {}
for i in range(0,data_length):
    date = full_data['date'][i]
    if not phases[date] in data_by_phase.keys():
        data_by_phase[phases[date]] = {}
        for k in full_data.keys():
            data_by_phase[phases[date]][k] = []
    for k in full_data.keys():
        data_by_phase[phases[date]][k].append(full_data[k][i])
    