"""
Combine Signature data and GPS data

"""
# -*- coding: utf-8 -*-
import os
import sys
from mhkit import dolfyn # Develop branch of the package, to install with "pip install git+https://github.com/MHKiT-Software/MHKiT-Python.git@develop"
import netCDF4
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import bottleneck as bn # For moving average
from datetime import datetime
from datetime import timezone
import cmocean
from scipy.interpolate import interp1d
import gpxpy

# plt.rcParams['svg.fonttype'] = 'none'
# plt.rcParams['font.size'] = 12
# plt.rcParams['font.family'] = 'Arial'  # or any installed font

plt.close('all')

#%% Specify field campaign here:
campaign_name="20250514"
#filename="S101024A034_burst_avge.ad2cp"
filename="S101024A033_moving_av.ad2cp"
# filename="S101024A032_fixed_avech.ad2cp"
# filename="S101024A028_fixed_avech.ad2cp"
savefig=False
# period_avg=[datetime(2025,5,14,22,50,0),datetime(2025,5,15,0,15,0)]
#period_avg=[datetime(2025,5,19,0,0,0),datetime(2025,5,21,0,0,0)]

dt_avg=10 # [s]
E_thresh=200
dt_interp=5 # [s]


#%% Load the ADCP data
data_Sign = dolfyn.read(os.path.join("..\..\data\Signature",campaign_name,"Level0",filename))
breakpoint()

#%% Load GPS data
# Open the .gpx file
with open(r"C:\Users\tdoda\OneDrive - Université de Lausanne\Fieldwork\Lake_Geneva\13-20250514_Practice_FieldCampaign\Data\GPS\Piste_GENZ TEST.gpx", 'r') as gpx_file:
    gpx = gpxpy.parse(gpx_file)


# Access data
GPS_points=gpx.tracks[0].segments[0].points
dateGPS=np.array([pt.time for pt in GPS_points]) # UTC
tnumGPS=np.array([dt.replace(tzinfo=timezone.utc).timestamp() for dt in dateGPS])+2*3600
longval=np.array([pt.longitude for pt in GPS_points])
latval=np.array([pt.latitude for pt in GPS_points])

            
                      
#%% Get echo data
data_echo=data_Sign[0]
tdate_echo=data_echo.time_echo.values
tnum_echo=tdate_echo.astype(np.float64)*1e-9 # [s]
depth_echo=data_echo.range_echo.values
E=data_echo.echo.values
temp=data_echo.temp_echo.values

# f=interp1d(tnum_echo,E,axis=1,kind='linear') 
# tnum_interp=np.arange(tnum_echo[0],tnum_echo[-1],dt_interp)
# tdate_interp=tnum_interp.astype("datetime64[s]")
# E=f(tnum_interp)
# tnum_echo=np.copy(tnum_interp)
# tdate_echo=np.copy(tdate_interp)

#%%
long_interp=np.interp(tnum_echo,tnumGPS,longval,left=np.nan,right=np.nan)
lat_interp=np.interp(tnum_echo,tnumGPS,latval,left=np.nan,right=np.nan)



#%% Get velocity data
data_avg=data_Sign[1]
tdate_avg=data_avg.time_avg.values
tnum_avg=tdate_avg.astype(np.float64)*1e-9 # [s]
depth_avg=data_avg.range_avg.values # [m]
VE=data_avg.vel_avg.values[0,:,:] # [m/s]
VN=data_avg.vel_avg.values[1,:,:] # [m/s]
VU1=data_avg.vel_avg.values[2,:,:] # [m/s]
VU2=data_avg.vel_avg.values[3,:,:] # [m/s]

#%% Correct the echo
Ecorr=np.copy(E)
Ecorr[np.abs(E)>E_thresh]=np.nan # Remove outliers

#%% Compute anomalies
periodnum_avg=[dt.replace(tzinfo=timezone.utc).timestamp() for dt in period_avg]
meanE=np.nanmean(Ecorr[:,(tnum_echo>periodnum_avg[0])&(tnum_echo<periodnum_avg[1])],axis=1)
dE=Ecorr-np.array([meanE,]).T
#dE_smooth=pd.DataFrame(dE).rolling(window=int(dt_avg/(tnum_echo[1]-tnum_echo[0])), axis=1, center=True).mean().to_numpy()
dE_smooth=bn.move_mean(dE, window=int(dt_avg/(tnum_echo[1]-tnum_echo[0])), axis=1, min_count=1)
dE_avg=np.nanmean(dE_smooth,axis=0)
#%% Plot Echo (low sampling period)
dt_plot=30 # sampling interval for the figure [s]
dcol_plot=int(dt_plot/(tnum_echo[1]-tnum_echo[0]))

fig,ax=plt.subplots(2,1,figsize=(10,5),sharex=True,sharey=True)
hp1=ax[0].pcolormesh(tdate_echo[np.arange(0,len(tnum_echo),dcol_plot)],depth_echo,Ecorr[:,np.arange(0,len(tnum_echo),dcol_plot)],vmin=0,vmax=100)
cb1=fig.colorbar(hp1)
ax[0].set_ylabel("Depth [m]")
cb1.set_label("Backscattering [dB]")

hp2=ax[1].pcolormesh(tdate_echo[np.arange(0,len(tnum_echo),dcol_plot)],depth_echo,dE_smooth[:,np.arange(0,len(tnum_echo),dcol_plot)],cmap=cmocean.cm.balance,vmin=-20,vmax=20)
cb2=fig.colorbar(hp2)
ax[1].set_ylabel("Depth [m]")
cb2.set_label("Backscattering anomaly [dB]")


ax[0].invert_yaxis()
fig.set_tight_layout(True)

#%% Plot Echo

fig,ax=plt.subplots(2,1,figsize=(10,5),sharex=True,sharey=True)
hp1=ax[0].pcolormesh(tdate_echo,depth_echo,Ecorr,vmin=0,vmax=100)
cb1=fig.colorbar(hp1)
ax[0].set_ylabel("Depth [m]")
cb1.set_label("Backscattering [dB]")

hp2=ax[1].pcolormesh(tdate_echo,depth_echo,dE_smooth,cmap=cmocean.cm.balance,vmin=-20,vmax=20)
cb2=fig.colorbar(hp2)
ax[1].set_ylabel("Depth [m]")
cb2.set_label("Backscattering anomaly [dB]")


ax[0].invert_yaxis()
fig.set_tight_layout(True)

#%% Plot Velocity

fig,ax=plt.subplots(4,1,figsize=(10,8),sharex=True,sharey=True)

hp1=ax[0].pcolormesh(tdate_avg,depth_avg,VE*100,cmap=cmocean.cm.balance,vmin=-100,vmax=100)
cb1=fig.colorbar(hp1)
ax[0].set_ylabel("Depth [m]")
cb1.set_label("East velocity [cm/s]")

hp2=ax[1].pcolormesh(tdate_avg,depth_avg,VN*100,cmap=cmocean.cm.balance,vmin=-100,vmax=100)
cb2=fig.colorbar(hp2)
ax[1].set_ylabel("Depth [m]")
cb2.set_label("East velocity [cm/s]")

hp3=ax[2].pcolormesh(tdate_avg,depth_avg,VU1*100,cmap=cmocean.cm.balance,vmin=-100,vmax=100)
cb3=fig.colorbar(hp3)
ax[2].set_ylabel("Depth [m]")
cb3.set_label("Upward velocity 1 [cm/s]")

hp4=ax[3].pcolormesh(tdate_avg,depth_avg,VU2*100,cmap=cmocean.cm.balance,vmin=-100,vmax=100)
cb4=fig.colorbar(hp4)
ax[3].set_ylabel("Depth [m]")
cb4.set_label("Upward velocity 1 [cm/s]")


ax[0].invert_yaxis()
fig.set_tight_layout(True)

#%% Plot spatial temperature and backscattering

fig,ax=plt.subplots(3,1,figsize=(7,8),sharex=True,sharey=True)

sc1=ax[0].scatter(long_interp,lat_interp,c=temp,cmap=cmocean.cm.thermal,s=10)
cb_temp=plt.colorbar(sc1)
cb_temp.set_label('Temp [°C]')
ax[0].set_aspect('equal')
ax[0].set_ylabel('Lat [°]')
ax[0].set_title('ADCP temperature')


dval=1.5 # [m]
sc2=ax[1].scatter(long_interp,lat_interp,c=dE_smooth[np.where(depth_echo>dval)[0][0],:],s=10,cmap=cmocean.cm.balance,vmin=-10,vmax=10)
cb_dE=plt.colorbar(sc2)
cb_dE.set_label('Echo anomaly [dB]')
ax[1].set_aspect('equal')
ax[1].set_ylabel('Lat [°]')
ax[1].set_title('Backscattering at {:.1f} m'.format(dval))

sc3=ax[2].scatter(long_interp,lat_interp,c=dE_avg,s=10,cmap=cmocean.cm.balance,vmin=-1,vmax=1)
cb_dEavg=plt.colorbar(sc3)
cb_dEavg.set_label('Echo anomaly [dB]')
ax[2].set_aspect('equal')
ax[2].set_xlabel('Long [°]')
ax[2].set_ylabel('Lat [°]')
ax[2].set_title('Depth-averaged backscattering')

ax[0].ticklabel_format(style='plain', useOffset=False)

fig.set_tight_layout(True)

# Export the figure
if savefig:
    fig.savefig('spatial_temp_backscat.png',dpi=400)
    # fig.savefig('spatial_temp_backscat.svg',dpi=400)

#%% Plot temperature-backscattering relationship (when GPS data)
fig,ax=plt.subplots(1,2)

ax[0].plot(temp[~np.isnan(long_interp)],dE_smooth[np.where(depth_echo>dval)[0][0],~np.isnan(long_interp)],'k.')

ax[1].plot(temp[~np.isnan(long_interp)],dE_avg[~np.isnan(long_interp)],'k.')

#%% Plot echo time series
d1=2
d2=10
d3=28

fig,ax=plt.subplots(2,1,figsize=(8,8),sharex=True)

ax[0].plot(tdate_echo,dE_smooth[np.where(depth_echo>d1)[0][0],:])
ax[0].plot(tdate_echo,dE_smooth[np.where(depth_echo>d2)[0][0],:])
#ax[0].plot(tdate_echo,dE_smooth[np.where(depth_echo>d3)[0][0],:])

ax[1].plot(tdate_echo,data_echo.temp_echo)

