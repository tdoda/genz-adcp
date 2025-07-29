"""
Read GPS data and export it to netCDF

"""
# -*- coding: utf-8 -*-
import os
import sys
import json
import yaml
import numpy as np
from gps_device import gps

#%% Specify field campaign and filename here:
campaigns_all=["20250514"]
files_all=["Piste_GENZ TEST.gpx"]
indices_data=[0] # Indices of the datafiles to process
gps_folder="..\..\data\GPS"

for index_file in indices_data:
    campaign_name=campaigns_all[index_file]
    filename=files_all[index_file]
    filename_noext=filename.replace(".gpx","")
    print("********************")
    print("Exporting file {}: {}".format(index_file,filename))

    #%% Load directories
    with open("input_python.yaml", "r") as f:
        directories = yaml.load(f, Loader=yaml.FullLoader) # Path to the data directories L0, L1, L2
       
    for key in directories.keys(): # Update path based on the field campaign
        directories[key]=os.path.join(gps_folder,campaign_name,directories[key])
        if not os.path.exists(directories[key]):
            os.makedirs(directories[key]) # Create the directories if they do not exist

    #%% Load metadata
    path_meta=os.path.join(gps_folder,campaign_name,"Level0",filename_noext+".meta")
    if os.path.isfile(path_meta):
        with open(path_meta, 'r') as f:
            parameter_dict = json.load(f) # Parameters related to the settings of the ADCP for each deployment period 
    else:
        raise Exception("Metadata file {} is missing!".format(filename_noext+".meta"))
        
    #%% Load the GPS data
    g = gps()
    # Add metadata information as general attributes:
    g.add_metadata(parameter_dict)
    if g.read_data(os.path.join(gps_folder,campaign_name,"Level0",filename)):  
        g.export(directories["Level1_dir"], "L1", output_period="file", time_label="time", overwrite_file=True) # Create Level 1 file 
    else:
        print("Could not export the data")




