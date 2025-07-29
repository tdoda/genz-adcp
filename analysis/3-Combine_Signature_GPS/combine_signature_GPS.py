"""
Combine Signature data and GPS data

"""
# -*- coding: utf-8 -*-
import os
import sys
import numpy as np
import pandas as pd
from datetime import datetime
from datetime import timezone
sys.path.append(os.path.join(os.getcwd(), '../../functions'))
from general_functions import  read_netCDF, export_netCDF

#%% Specify field campaign and filename here:
campaigns_all=["20250514"]
gps_folder="..\..\data\GPS"
signature_folder="..\..\data\Signature"
indices_data=[0] # Indices of the datafiles to process

for index_file in indices_data:
    campaign_name=campaigns_all[index_file]
    
    print("********************")
    print("Campaign {}".format(campaign_name))

    #%% Load the ADCP data
    filepath=os.path.join(signature_folder,campaign_name,"Level2")
    filenames=os.listdir(filepath) # All files
    if len(filenames)>1:
        raise Exception("More than one Signature datafile in campaign {}".format(campaign_name))
    filename_str=filenames[0]
    
    Sign_data, Sign_genatt, Sign_varatt, Sign_dim=read_netCDF(os.path.join(filepath,filename_str))
    
    #%% Load GPS data
    filepath=os.path.join(gps_folder,campaign_name,"Level1")
    filenames=os.listdir(filepath) # All files
    if len(filenames)>1:
        raise Exception("More than one GPS datafile in campaign {}".format(campaign_name))
    filename_str=filenames[0]
    
    GPS_data, GPS_genatt, GPS_varatt, GPS_dim=read_netCDF(os.path.join(filepath,filename_str))
    
    #%% Interpolate longitude and latitude values to the echosounder time data
    tzone_Sign=int(Sign_genatt["campaign_Time Zone device (UTC+)"])
    tzone_GPS=int(GPS_genatt["campaign_Time Zone device (UTC+)"])
    lon_interp=np.interp(Sign_data["time_echo"],GPS_data["time"]+(tzone_Sign-tzone_GPS)*3600,GPS_data["lon"],left=np.nan,right=np.nan)
    lat_interp=np.interp(Sign_data["time_echo"],GPS_data["time"]+(tzone_Sign-tzone_GPS)*3600,GPS_data["lat"],left=np.nan,right=np.nan)
    
    Sign_data["lon_vel"]=np.interp(Sign_data["time_vel"],GPS_data["time"]+(tzone_Sign-tzone_GPS)*3600,GPS_data["lon"],left=np.nan,right=np.nan)
    Sign_varatt["lon_vel"]={"var_name":"lon_vel","unit":GPS_varatt["lon"]["units"],
                            "long_name":"Longitude of the velocity data","dim":("time_vel",)}
    
    Sign_data["lat_vel"]=np.interp(Sign_data["time_vel"],GPS_data["time"]+(tzone_Sign-tzone_GPS)*3600,GPS_data["lat"],left=np.nan,right=np.nan)
    Sign_varatt["lat_vel"]={"var_name":"lat_vel","unit":GPS_varatt["lat"]["units"],
                            "long_name":"Latitude of the velocity data","dim":("time_vel",)}
    
    #%% Separate echograms
    n_echo=int(Sign_genatt["n_echograms"])
    ind_echo=[]
    if n_echo==2:
        ind_echo.append(np.arange(0,len(Sign_data["time_echo"]),2))
        ind_echo.append(np.arange(1,len(Sign_data["time_echo"]),2))
    elif n_echo==1:
        ind_echo.append(np.arange(0,len(Sign_data["time_echo"]),1))
    else:
        raise Exception("Incorrect number of echograms")
        
    for ke in range(n_echo):
 
        Sign_data["time_echo"+str(ke+1)]=Sign_data["time_echo"][ind_echo[ke]]
        Sign_varatt["time_echo"+str(ke+1)]={"var_name":"time_echo"+str(ke+1),
                                            "unit":Sign_varatt["time_echo"]["units"],
                                            "long_name":"Time of the echogram "+str(ke+1),
                                            "dim":("time_echo"+str(ke+1),)}
        Sign_dim["time_echo"+str(ke+1)]={"dim_name":"time_echo"+str(ke+1),"dim_size":None}
        
        Sign_data["time_echo"+str(ke+1)+"_qual"]=Sign_data["time_echo_qual"][ind_echo[ke]]
        Sign_varatt["time_echo"+str(ke+1)+"_qual"]={"var_name":"time_echo"+str(ke+1)+"_qual",
                                            "unit":Sign_varatt["time_echo_qual"]["units"],
                                            "long_name":"time_echo"+str(ke+1)+"_qual",
                                            "dim":("time_echo"+str(ke+1),)}
        
        Sign_data["echo"+str(ke+1)]=Sign_data["echo_HR"][:,ind_echo[ke]]
        Sign_varatt["echo"+str(ke+1)]={"var_name":"echo"+str(ke+1),
                                            "unit":Sign_varatt["echo_HR"]["units"],
                                            "long_name":"Echogram "+str(ke+1),
                                            "dim":("depth_echo","time_echo"+str(ke+1))}
        
        Sign_data["lon_echo"+str(ke+1)]=lon_interp[ind_echo[ke]]
        Sign_varatt["lon_echo"+str(ke+1)]={"var_name":"lon_echo"+str(ke+1),
                                            "unit":GPS_varatt["lon"]["units"],
                                            "long_name":"Longitude of echogram "+str(ke+1),
                                            "dim":("time_echo"+str(ke+1),)}
        
        
        Sign_data["lat_echo"+str(ke+1)]=lat_interp[ind_echo[ke]]
        Sign_varatt["lat_echo"+str(ke+1)]={"var_name":"lat_echo"+str(ke+1),
                                            "unit":GPS_varatt["lat"]["units"],
                                            "long_name":"Latitude of echogram "+str(ke+1),
                                            "dim":("time_echo"+str(ke+1),)}
        
        Sign_data["ind_echo"+str(ke+1)]=ind_echo[ke]
        Sign_varatt["ind_echo"+str(ke+1)]={"var_name":"ind_echo"+str(ke+1),
                                            "unit":"-",
                                            "long_name":"Indices of echogram "+str(ke+1)+" in complete burst echo data",
                                            "dim":("time_echo"+str(ke+1),)}
    
        
    #del Sign_data["time_echo"]; del Sign_varatt["time_echo"]; del Sign_dim["time_echo"]
    #del Sign_data["time_echo_qual"]; del Sign_varatt["time_echo_qual"]
    del Sign_data["echo_HR"]; del Sign_varatt["echo_HR"]
    
            
    
    #%% Export to netCDF
    export_netCDF("Signature_GPS_"+campaign_name+".nc",Sign_genatt,Sign_dim,Sign_varatt,Sign_data)
    
    

    
