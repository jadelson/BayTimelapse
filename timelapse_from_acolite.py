#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 13:20:34 2017

@author: jadelson
"""
import pandas as pd
from datetime import timedelta
import urllib
import json
import bisect
import numpy as np
import matplotlib.pyplot as plt
import pickle
import imageio

#%%
#with open('/Users/jadelson/Dropbox/phdResearch/AllOptical/timelapse_work/full_dataset.dk','rb') as f:
#    full_data = pickle.load(f)
data_length = len(full_data['tau'])   
date_list = []

rho = np.asarray(full_data['RHOW_865'])
A = 3031.75 
B = 0
C = 0.2114
full_data['ssc'] = A*rho/(1 - rho/C) + B
    
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

#%%        
        
low_water_times = tide_df.index[(tide_df['ty']=='LL') | (tide_df['ty']=='L ')]
ll_df = tide_df[(tide_df['ty']=='LL') | (tide_df['ty']=='L ')]
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

with open('data_by_phase.dk','wb') as f:
    pickle.dump(data_by_phase, f)
    
#%%

with open('data_by_phase.dk','rb') as f:
    data_by_phase = pickle.load(f)
    
phase_list = np.asarray(list(data_by_phase.keys()))
phase_list = np.sort(phase_list)
#for phase in phase_list:
#    x = data_by_phase[phase]['lon']
#    y = data_by_phase[phase]['lat']
#    ssc = data_by_phase[phase]['ssc']
#    
#    ax = plt.subplot(1,1,1)
#    plt.clf()
#    plt.title(phase)
#    plt.xlim([-122.55, -121.75])
#    plt.ylim([37.4, 38.3])
#    ax.set_autoscaley_on(False)
#    plt.scatter(x,y,c=ssc,s=.1, vmin=0, vmax=150)
#    plt.colorbar()
#    plt.savefig('figures/ssc_'+str(phase)+'.png', bbox_inches='tight')
    
   
#%%
hourmax = 12.421
tmax = hourmax*3600
nbuckets = 8
dt = tmax/nbuckets/2
bins = np.linspace(dt, tmax+dt, nbuckets+1)
digitized = np.digitize(phase_list, bins)
digitized[digitized==nbuckets] = 0
bin_counts = [len(phase_list[digitized == i]) for i in range(0, nbuckets)]
bin_means = [phase_list[digitized == i].mean()/3600 for i in range(0, nbuckets)]
filenames = []

for i in range(0,nbuckets):
    if sum(digitized==i) == 0:
        continue
    data = {}
    phase = str(np.round((bins[i] - dt)/3600, 1))
    name = 'figures/median_ssc_'+phase+'.png'
    for p in phase_list[digitized == i]:
        v = data_by_phase[p]
        for j in range(0,len(v['lat'])):
            la = np.round(v['lat'][j],4)
            lo = np.round(v['lon'][j],4)
            if not (la,lo) in data.keys():
                data[(la,lo)] = []
            data[(la,lo)].append(v['ssc'][j])
    lat = []
    lon = []
    ssc = []
    
    for (la,lo) in data.keys():
        lat.append(la)
        lon.append(lo)
        ssc.append(np.median(data[la,lo]))
    plt.clf()
    plt.title(phase)
    plt.scatter(lon,lat,c=ssc,s=.1, vmin=0, vmax=150)
    plt.colorbar()
    plt.xlim([-122.55, -121.75])
    plt.ylim([37.4, 38.3])
    
    a = plt.axes([.5, .35, .22, .22])
    plt.plot(np.arange(0,hourmax), -np.cos(np.arange(0,hourmax)/hourmax*2*np.pi))
    plt.scatter(float(phase),-np.cos(float(phase)/hourmax*2*np.pi),color='red')    
#    plt.title('Tide')
    plt.xticks([])
    plt.yticks([])
    plt.savefig(name, bbox_inches='tight', format='png', dpi=800)
    filenames.append(name)

#%%
with imageio.get_writer('figures/tide.gif', mode='I') as writer:
    for filename in filenames:
        for k in range(0,5):
            image = imageio.imread(filename)
            writer.append_data(image)
    









