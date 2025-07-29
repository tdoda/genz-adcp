"""
Plot spatial variability of Signature data

"""
# -*- coding: utf-8 -*-
import os
import sys
import numpy as np
import pandas as pd
from datetime import datetime, timezone
sys.path.append(os.path.join(os.getcwd(), '../../functions'))
from general_functions import  read_netCDF
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import cmocean
import bottleneck as bn # For moving average

plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.size'] = 12
plt.rcParams['font.family'] = 'Arial'  # or any installed font

plt.close('all')

#%% Specify filename here:
foldername="../../analysis/3-Combine_Signature_GPS"
filename="Signature_GPS_20250514.nc"

# Parameters:
dval=1.5 # [m] # Depth of echo data:
dt_avg=10 # [s] # Temporal moving average
savefig=True # True to save the figure
period_avg=[datetime(2025,5,14,22,50,0),datetime(2025,5,15,0,15,0)] # Period when the data should be considered
dt=60 # Temporal resolution shown on the figure [sec], to plot the data faster
#%% Load data
Sign_data, Sign_genatt, Sign_varatt, Sign_dim=read_netCDF(os.path.join(foldername,filename))
temp_echo1=Sign_data["temp"][Sign_data["ind_echo1"].astype("int")]

#%% Compute anomalies
periodnum_avg=[dt.replace(tzinfo=timezone.utc).timestamp() for dt in period_avg]
ind_period=np.where((Sign_data["time_echo1"]>periodnum_avg[0])&(Sign_data["time_echo1"]<periodnum_avg[1]))[0]

meanE=np.nanmean(Sign_data["echo1"][:,ind_period],axis=1)
dE=Sign_data["echo1"]-np.array([meanE,]).T

dE_smooth=bn.move_mean(dE, window=int(dt_avg/(Sign_data["time_echo1"][1]-Sign_data["time_echo1"][0])), axis=1, min_count=1)
dE_avg=np.nanmean(dE_smooth,axis=0) # Depth-averaged

dval_echo=int(dt/(Sign_data["time_echo1"][1]-Sign_data["time_echo1"][0])) # Temporal index step
#%% Plot backscattering time-series
fig,ax=plt.subplots(2,1,figsize=(7,8),sharex=True,sharey=True,layout="constrained")

pm1=ax[0].pcolormesh(pd.to_datetime(Sign_data["time_echo1"][::dval_echo],unit="s"),Sign_data["depth_echo"],Sign_data["echo1"][:,::dval_echo])
cb_echo=plt.colorbar(pm1)
cb_echo.set_label('Backscattering [dB]')
ax[0].set_ylabel('Depth [m]')


pm2=ax[1].pcolormesh(pd.to_datetime(Sign_data["time_echo1"][::dval_echo],unit="s"),Sign_data["depth_echo"],dE_smooth[:,::dval_echo],cmap=cmocean.cm.balance,vmin=-10,vmax=10)
cb_echo=plt.colorbar(pm2)
cb_echo.set_label('Backscattering difference [dB]')
ax[1].set_ylabel('Depth [m]')
ax[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax[1].tick_params(axis='x', labelrotation=20)

ax[0].invert_yaxis()

# Export the figure
if savefig:
    fig.savefig('time_series_backscat.png',dpi=400)  
            
#%% Plot spatial temperature and backscattering
fig,ax=plt.subplots(3,1,figsize=(7,8),sharex=True,sharey=True,layout="constrained")

sc1=ax[0].scatter(Sign_data["lon_echo1"],Sign_data["lat_echo1"],c=temp_echo1,cmap=cmocean.cm.thermal,s=10)
cb_temp=plt.colorbar(sc1)
cb_temp.set_label('Temp [°C]')
ax[0].set_aspect('equal')
ax[0].set_ylabel('Lat [°]')
ax[0].set_title('ADCP temperature')


sc2=ax[1].scatter(Sign_data["lon_echo1"],Sign_data["lat_echo1"],c=dE_smooth[np.where(Sign_data["depth_echo"]>dval)[0][0],:],s=10,cmap=cmocean.cm.balance,vmin=-10,vmax=10)
cb_dE=plt.colorbar(sc2)
cb_dE.set_label('Echo anomaly [dB]')
ax[1].set_aspect('equal')
ax[1].set_ylabel('Lat [°]')
ax[1].set_title('Backscattering at {:.1f} m'.format(dval))

sc3=ax[2].scatter(Sign_data["lon_echo1"],Sign_data["lat_echo1"],c=dE_avg,s=10,cmap=cmocean.cm.balance,vmin=-1,vmax=1)
cb_dEavg=plt.colorbar(sc3)
cb_dEavg.set_label('Echo anomaly [dB]')
ax[2].set_aspect('equal')
ax[2].set_xlabel('Long [°]')
ax[2].set_ylabel('Lat [°]')
ax[2].set_title('Depth-averaged backscattering')

ax[0].ticklabel_format(style='plain', useOffset=False)

# Export the figure
if savefig:
    fig.savefig('spatial_temp_backscat.png',dpi=400)

    
    

    
