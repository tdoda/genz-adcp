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
from general_functions import  read_netCDF

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
    
    Sign_data, Sign_genatt, Sign_varatt=read_netCDF(os.path.join(filepath,filename_str))
    
    #%% Load GPS data
    filepath=os.path.join(gps_folder,campaign_name,"Level1")
    filenames=os.listdir(filepath) # All files
    if len(filenames)>1:
        raise Exception("More than one GPS datafile in campaign {}".format(campaign_name))
    filename_str=filenames[0]
    
    GPS_data, GPS_genatt, GPS_varatt=read_netCDF(os.path.join(filepath,filename_str))
    
    #%% Interpolate longitude and latitude values to the echosounder time data
    tzone_Sign=int(Sign_genatt["campaign_Time Zone device (UTC+)"])
    tzone_GPS=int(GPS_genatt["campaign_Time Zone device (UTC+)"])
    long_interp=np.interp(Sign_data["time_echo"],GPS_data["time"]+(tzone_Sign-tzone_GPS)*3600,GPS_data["lon"],left=np.nan,right=np.nan)
    lat_interp=np.interp(Sign_data["time_echo"],GPS_data["time"]+(tzone_Sign-tzone_GPS)*3600,GPS_data["lat"],left=np.nan,right=np.nan)
