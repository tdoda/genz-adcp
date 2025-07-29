"""
Read Signature data and export it to netCDF

"""
# -*- coding: utf-8 -*-
import os
import sys
import json
import yaml
import numpy as np
from adcp_signature import ADCP_Signature

#%% Specify field campaign and filename here:
campaigns_all=["20250403" ,"20250411","20250430","20250514"  ,"20250519","20250520"]
files_all=["S101024A027_float2_noav.ad2cp","S101024A028_fixed_avech.ad2cp","S101024A032_fixed_avech.ad2cp","S101024A033_moving_av.ad2cp","S101024A034_burst_avge.ad2cp","S101024A035_burst_avge.ad2cp"]
indices_data=np.arange(1,6) # Indices of the datafiles to process
# indices_data=[1]
signature_folder="..\..\data\Signature"

for index_file in indices_data:
    campaign_name=campaigns_all[index_file]
    filename=files_all[index_file]
    filename_noext=filename.replace(".ad2cp","")
    print("********************")
    print("Exporting file {}: {}".format(index_file,filename))

    #%% Load directories
    with open("input_python.yaml", "r") as f:
        directories = yaml.load(f, Loader=yaml.FullLoader) # Path to the data directories L0, L1, L2
       
    for key in directories.keys(): # Update path based on the field campaign
        directories[key]=os.path.join(signature_folder,campaign_name,directories[key])
        if not os.path.exists(directories[key]):
            os.makedirs(directories[key]) # Create the directories if they do not exist

    #%% Load metadata
    path_meta=os.path.join(signature_folder,campaign_name,"Level0",filename_noext+".meta")
    if os.path.isfile(path_meta):
        with open(path_meta, 'r') as f:
            parameter_dict = json.load(f) # Parameters related to the settings of the ADCP for each deployment period 
    else:
        raise Exception("Metadata file {} is missing!".format(filename_noext+".meta"))
        
    #%% Load the ADCP data
    a = ADCP_Signature()
    # Add metadata information as general attributes:
    a.add_metadata(parameter_dict)
    if a.read_data(os.path.join(signature_folder,campaign_name,"Level0",filename), transducer_depth=parameter_dict["transducer_depth (m)"], 
                burst=parameter_dict["burst"], concurrent=parameter_dict["concurrent"], up=parameter_dict["up"],start_date=parameter_dict["start_date"],end_date=parameter_dict["end_date"]):
        
        a.quality_flags("quality_assurance.json",concurrent=parameter_dict["concurrent"]) # Flag the data based on quality checks  
        a.export(directories["Level1_dir"], "L1", output_period="file", time_label="time_echo", overwrite_file=True) # Create Level 1 file
        a.mask_data() # Replace flagged data by nan 
        a.export(directories["Level2_dir"], "L2", output_period="file", time_label="time_echo", overwrite_file=True) # Create Level 2 file
    else:
        print("Could not export the data")




