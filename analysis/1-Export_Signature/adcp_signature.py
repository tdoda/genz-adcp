# -*- coding: utf-8 -*-
import os
import sys
import netCDF4
import numpy as np
import xarray as xr
#from mhkit import dolfyn # Develop branch of the package, to install with "pip install git+https://github.com/MHKiT-Software/MHKiT-Python.git@develop"
import dolfyn_modified as dolfyn # Modified version of the dolfyn package to deal with different types of Signature datasets
import distutils.util
from functions import json_converter, log
from datetime import datetime
import json
from envass_modified import qualityassurance
from dateutil.relativedelta import relativedelta
sys.path.append(os.path.join(os.path.dirname(__file__), '../../functions'))
from general_functions import GenericInstrument
from quality_checks_adcp import *


class ADCP_Signature(GenericInstrument):
    def __init__(self, *args, **kwargs):
        super(ADCP_Signature, self).__init__(*args, **kwargs)
        self.general_attributes = {
            "institution": "Eawag",
            "source": "ADCP Signature",
            "references": "Floating Signature tomy.doda@unil.ch",
            "history": "See history on Renku",
            "conventions": "CF 1.7",
            "comment": "Data from floating Nortek Signature ADCP in Lake Geneva",
            "title": "Signature backscattering and velocities",
            "Instrument": "Nortek Signature ADCP",
        }
        
        self.dimensions = {
            'time_echo': {'dim_name': 'time_echo', 'dim_size': None},
            'time_vel': {'dim_name': 'time_vel', 'dim_size': None},
            'depth_echo': {'dim_name': 'depth_echo', 'dim_size': None},
            'depth_vel': {'dim_name': 'depth_vel', 'dim_size': None}
        }
        
        self.variables = {
            'time_echo': {'var_name': 'time_echo', 'dim': ('time_echo',), 'unit': 'seconds since 1970-01-01 00:00:00', 'long_name': 'Time of the burst echosounder measurements'},
            'time_vel': {'var_name': 'time_vel', 'dim': ('time_vel',), 'unit': 'seconds since 1970-01-01 00:00:00', 'long_name': 'Time of the burst velocity data'},
            'depth_echo': {'var_name': 'depth_echo', 'dim': ('depth_echo',), 'unit': 'm', 'long_name': 'Depth of the echosounder data'},
            'depth_vel': {'var_name': 'depth_vel', 'dim': ('depth_vel',), 'unit': 'm', 'long_name': 'Depth of the velocity data'},
            'u': {'var_name': 'u', 'dim': ('depth_vel', 'time_vel'), 'unit': 'm s-1', 'long_name': 'eastern velocity'},
            'v': {'var_name': 'v', 'dim': ('depth_vel', 'time_vel'), 'unit': 'm s-1', 'long_name': 'northern velocty'},
            'w': {'var_name': 'w', 'dim': ('depth_vel', 'time_vel'), 'unit': 'm s-1', 'long_name': 'upward velocty'},
            'temp': {'var_name': 'temp', 'dim': ('time_echo',), 'unit': 'degC', 'long_name': 'temperature', },
            'press': {'var_name': 'press', 'dim': ('time_echo',), 'unit': 'dbar', 'long_name': 'pressure', },
            'echo_HR': {'var_name': 'echo_HR', 'dim': ('depth_echo', 'time_echo'), 'unit': 'dB', 'long_name': 'Echosounder data'},
            'accE': {'var_name': 'accE', 'dim': ('time_echo',), 'unit': 'm.s-2', 'long_name': 'Acceleration along the East direction'},
            'accN': {'var_name': 'accN', 'dim': ('time_echo',), 'unit': 'm.s-2', 'long_name': 'Acceleration along the North direction'},
            'accU': {'var_name': 'accU', 'dim': ('time_echo',), 'unit': 'm.s-2', 'long_name': 'Acceleration along the upward direction'},
            'echo1': {'var_name': 'echo1', 'dim': ('depth_vel', 'time_vel'), 'unit': 'dB', 'long_name': 'Beam 1 echo'},
            'echo2': {'var_name': 'echo2', 'dim': ('depth_vel', 'time_vel'), 'unit': 'dB', 'long_name': 'Beam 2 echo'},
            'echo3': {'var_name': 'echo3', 'dim': ('depth_vel', 'time_vel'), 'unit': 'dB', 'long_name': 'Beam 3 echo'},
            'echo4': {'var_name': 'echo4', 'dim': ('depth_vel', 'time_vel'), 'unit': 'dB', 'long_name': 'Beam 4 echo'},
            'corr1': {'var_name': 'corr1', 'dim': ('depth_vel', 'time_vel'), 'unit': '-', 'long_name': 'Beam 1 correlation', },
            'corr2': {'var_name': 'corr2', 'dim': ('depth_vel', 'time_vel'), 'unit': '-', 'long_name': 'Beam 2 correlation', },
            'corr3': {'var_name': 'corr3', 'dim': ('depth_vel', 'time_vel'), 'unit': '-', 'long_name': 'Beam 3 correlation', },
            'corr4': {'var_name': 'corr4', 'dim': ('depth_vel', 'time_vel'), 'unit': '-', 'long_name': 'Beam 4 correlation', },
            'battery': {'var_name': 'battery', 'dim': ('time_echo',), 'unit': '-', 'long_name': 'Battery level'},
            'heading': {'var_name': 'heading', 'dim': ('time_vel',), 'unit': 'deg', 'long_name': 'Heading for the velocity data'},
            'roll': {'var_name': 'roll', 'dim': ('time_vel',), 'unit': 'deg', 'long_name': 'Roll for the velocity data'},
            'pitch': {'var_name': 'pitch', 'dim': ('time_vel',), 'unit': 'deg', 'long_name': 'Pitch for the velocity data'},
            'heading_echo': {'var_name': 'heading_echo', 'dim': ('time_echo',), 'unit': 'deg', 'long_name': 'Heading for the echosounder data'},
            'roll_echo': {'var_name': 'roll_echo', 'dim': ('time_echo',), 'unit': 'deg', 'long_name': 'Roll for the echosounder data'},
            'pitch_echo': {'var_name': 'pitch_echo', 'dim': ('time_echo',), 'unit': 'deg', 'long_name': 'Pitch for the echosounder data'},
        }

        self.derived_variables = {
        }

        self.data = {}

    def read_data(self, file, transducer_depth,  burst=True, concurrent=True, up=False, **kwargs):
        """
        Read the ADCP data and store it in an ADCP object.

        Parameters:
            file (str): path and ADCP filename (.ad2cp file)
            transducer_depth (float): depth of the ADCP [m]
            burst (bool): =True if burst data, =False if averaged data only
            concurrent (bool): =True if concurrent plans (burst + average), =False if continuous data only (burst only or averaged only)
            up (bool): =True is the ADCP is upward-looking, =False if the ADCP is downward-looking
            
        Additional arguments (**kwargs):
            start_date (str, %Y%m%d %H:%M format): starting time of the period to extract
            end_date (str, %Y%m%d %H:%M format): end time of the period to extract
                
        Returns:
            True if the data was correctly read, False otherwise
        """
        
        log("Parsing data from {}.".format(file))
        
        if not isinstance(burst, bool):
            burst = bool(distutils.util.strtobool(burst))
            
        if not isinstance(concurrent, bool):
            concurrent = bool(distutils.util.strtobool(concurrent))
            
        
        if not burst:
            print("Averaged data cannot be read by the current version of the script")
            return False
            
        if not concurrent:
            print("Single plan cannot be read by the current version of the script")
            return False

        try:
            data_Sign = dolfyn.read(file) # Returns a tuple of two datasets if the file is a burst data file: data_Sign[0] is the burst echo data, data_Sign[1] is the burst velocity data (even if name of variables are "*_avg")
            if len(data_Sign)!=2:
                raise Exception("Wrong data format")
            date_echo= data_Sign[0].time_echo.data.astype('datetime64[s]')
            time_echo= data_Sign[0].time_echo.data.astype('datetime64[ns]').astype(np.float64)*1e-9
            date_vel= data_Sign[1].time_avg.data.astype('datetime64[s]')
            time_vel= data_Sign[1].time_avg.data.astype('datetime64[ns]').astype(np.float64)*1e-9

            if not isinstance(up, bool):
                up = bool(distutils.util.strtobool(up))

            # Define the measurement period from specified start and end dates:
            date_subset = []
            if "start_date" in kwargs and kwargs["start_date"]: # Not empty
                start_date = datetime.strptime(kwargs["start_date"], "%Y%m%d %H:%M")
                date_subset.append(start_date)
            else:
                date_subset.append(date_echo[0])
                
            if "end_date" in kwargs and kwargs["end_date"]: # Not empty
                end_date = datetime.strptime(kwargs["end_date"], "%Y%m%d %H:%M")
                date_subset.append(end_date)
            else:
                date_subset.append(date_echo[-1])

            time_subset = np.array(date_subset).astype('int')
            idx_subset_echo = (time_subset[0] <= time_echo) & (time_echo <= time_subset[1])
            idx_subset_vel = (time_subset[0] <= time_vel) & (time_vel <= time_subset[1])
            
            # Create a new dataset with the subset data:
            dlfn_subset=xr.Dataset(coords={key: value.values for key, value in data_Sign[0].coords.items() if key not in ["time","time_echo"]})
            for key,value in data_Sign[1].coords.items():
                if key not in ["time","time_avg","range_avg"] and key not in list(data_Sign[0].coords):
                    dlfn_subset=dlfn_subset.expand_dims({key:value})
            dlfn_subset=dlfn_subset.expand_dims({"time_echo":time_echo[idx_subset_echo]})
            dlfn_subset=dlfn_subset.expand_dims({"time_vel":time_vel[idx_subset_vel]})
            dlfn_subset=dlfn_subset.expand_dims({"range_vel":data_Sign[1].range_avg.data})
            #dlfn_subset=dlfn_subset.assign(timestamp=(["time"],time[idx_subset]))
            for var in data_Sign[0].variables:
                dimvar=list(data_Sign[0][var].dims)
                if "time_echo" in dimvar:
                    print(var)
                    dlfn_subset=dlfn_subset.assign({var:(dimvar,data_Sign[0][var].isel(time_echo=idx_subset_echo).values)})
            for var in data_Sign[1].variables:
                if var!="time": # Dimensions time and time_avg are the same in data_Sign[1]
                    print(var)
                    dimvar=list(data_Sign[1][var].dims)
                    dimvar=["range_vel" if d=="range_avg" else d for d in dimvar]
                    dimvar=["time_vel" if d=="time_avg" else d for d in dimvar]
                    
                    if "_avg" in var:
                        if 'vel' not in var:
                            var_new=var.replace('_avg','_vel')  
                        else:
                            var_new=var.replace('_avg','')  
                        
                    else:
                        if 'vel' not in var:
                            var_new=var+'_vel'
                        else:
                            var_new=var.copy()
                    if "time_vel" in dimvar:
                            dlfn_subset=dlfn_subset.assign({var_new:(dimvar,data_Sign[1][var].isel(time_avg=idx_subset_vel).values)})
                    elif "time" in dimvar:
                            dimvar=["time_vel" if d=="time" else d for d in dimvar]
                            dlfn_subset=dlfn_subset.assign({var_new:(dimvar,data_Sign[1][var].isel(time=idx_subset_vel).values)})
            
            self.general_attributes["Er"] = str(np.nanmin(dlfn_subset.amp_vel.values.astype(float)))
            self.general_attributes["up"] = str(up)
            #self.general_attributes["transducer_depth"] = transducer_depth
            #self.general_attributes["xmit_length"] = float(data_Sign[0].attrs["transmit_pulse_m"]) # transmit pulse length [m] used for backscattering calculation 
            self.general_attributes["beam_angle"] = str(data_Sign[1].attrs["beam_angle"]) 
            self.general_attributes["blank_dist"] = str(data_Sign[1].attrs["blank_dist_avg"])
            #self.general_attributes["beam_freq"] = float(dlfn_data.attrs["freq"])
            
            if up:
                z0_echo = transducer_depth - dlfn_subset.range_echo.values
                z0_vel = transducer_depth - dlfn_subset.range_vel.values
            else:
                z0_echo = transducer_depth + dlfn_subset.range_echo.values
                z0_vel = transducer_depth + dlfn_subset.range_vel.values
           
            self.data["time_echo"] = time_echo[idx_subset_echo]
            self.data["time_vel"] = time_vel[idx_subset_vel]
            self.data["depth_echo"] = z0_echo
            self.data["depth_vel"] = z0_vel
            
            #self.data["r"] = z0/np.cos(self.general_attributes["beam_angle"]*np.pi/180.) # radial distance [m]
            #self.data["eu"] = dlfn_subset.vel[3, :, :].values # Velocity error

            self.data["corr"] = dlfn_subset.corr_vel.values.astype(float)/100. # Between 0 and 1
            self.data["corr1"] = self.data["corr"][0, :, :]
            self.data["corr2"] = self.data["corr"][1, :, :]
            self.data["corr3"] = self.data["corr"][2, :, :]
            self.data["corr4"] = self.data["corr"][3, :, :]
            self.data["heading_echo"] = dlfn_subset.heading_echo.values
            self.data["roll_echo"] = dlfn_subset.roll_echo.values
            self.data["pitch_echo"] = dlfn_subset.pitch_echo.values
            self.data["heading"] = dlfn_subset.heading_vel.values
            self.data["roll"] = dlfn_subset.roll_vel.values
            self.data["pitch"] = dlfn_subset.pitch_vel.values
            self.data["u"] = dlfn_subset.vel[0, :, :].values # Eastward velocity [m/s]
            self.data["v"] = dlfn_subset.vel[1, :, :].values # Northward velocity [m/s]
            self.data["w"] = dlfn_subset.vel[2, :, :].values # Upward velocity [m/s]
            self.data["echo_HR"] = dlfn_subset.echo.values
            self.data["accE"] = dlfn_subset.accel_echo.values[0,:]
            self.data["accN"] = dlfn_subset.accel_echo.values[1,:]
            self.data["accU"] = dlfn_subset.accel_echo.values[2,:]
            
            echo = dlfn_subset.amp_vel.values.astype(float)
            self.data["echo1"] = echo[0, :, :]
            self.data["echo2"] = echo[1, :, :]
            self.data["echo3"] = echo[2, :, :]
            self.data["echo4"] = echo[3, :, :]
            self.data["battery"]=dlfn_subset.batt_echo.values
            self.data["temp"]=dlfn_subset.temp_echo.values
            self.data["press"]=dlfn_subset.pressure_echo.values

            return True
        except:
            log("Failed to process {}.".format(file))
            return False

    def quality_flags(self, file_path = 'quality_assurance.json', simple=True, concurrent=True):
        
        if not isinstance(concurrent, bool):
            concurrent = bool(distutils.util.strtobool(concurrent))
            
        if not concurrent:
            raise Exception("Single plan cannot be quality checked by the current version of the script")    
        

        log("Performing quality assurance")
        log("1. ADCP-specific quality checks",indent=1) # Additional ADCP tests on the velocity matrix, increases qa with a base-2 approach (check#2 returns 0 or 2, chech#3 returns 0 or 4, etc.)
        qa_vel=init_flag_adcp(np.array(self.data["u"])) # Initial qa array (zero values)
        qa_echo=init_flag_adcp(np.array(self.data["echo_HR"])) # Initial qa array (zero values)
        # if self.general_attributes['up']=='True': # Upward looking: surface detection
        #     qa_adcp=qa_adcp_interface_top(qa_adcp,self.data["depth_vel"],self.general_attributes['transducer_depth'],beam_angle=25)
        # else: # Downward looking: sediment detection
        #     qa_adcp=qa_adcp_interface_bottom(qa_adcp,self.data["depth_vel"],self.general_attributes['transducer_depth'],self.general_attributes['bottom_depth'],beam_angle=25)
            
        qa_vel=qa_adcp_corr(qa_vel,self.data["corr1"],self.data["corr2"],self.data["corr3"],self.data["corr4"],corr_threshold=70)
        qa_vel=qa_adcp_tilt(qa_vel,self.data["roll"],self.data["pitch"],tilt_threshold=15)
        qa_echo=qa_adcp_tilt(qa_echo,self.data["roll_echo"],self.data["pitch_echo"],tilt_threshold=15)
        #qa_vel=qa_adcp_corrstd(qa_vel,self.data["corr1"],self.data["corr2"],self.data["corr3"],self.data["corr4"],std_threshold=0.02)
        #qa_adcp=qa_adcp_echodiff(qa_adcp,self.data["echo1"],self.data["echo2"],self.data["echo3"],self.data["echo4"],diff_threshold=30)
        
        log("2. envass quality checks",indent=1) # Corresponds to quality check #1: qa is 0 (all good) or 1 (flagged)
        quality_assurance_dict = json_converter(json.load(open(file_path))) # Load parameters related to simple and advanced quality checks
        
        for key, value in self.variables.copy().items():  
            if (key in quality_assurance_dict):  
                print(key)
                name = key + "_qual" # Create a new variable _qual to flag the data
                self.variables[name] = {'var_name': name, 'dim': value["dim"],
                                        'unit': '0 = nothing to report, 1 = more investigation',
                                        'long_name': name, }
                data_array = self.data[key]
                if not isinstance(data_array, np.ndarray):
                    data_array = np.array(data_array)
                    
                shape_var=data_array.shape
                if shape_var==qa_vel.shape:
                    prior_qa=qa_vel
                    timeval=self.data["time_vel"]
                elif shape_var==qa_echo.shape:
                    prior_qa=qa_echo
                    timeval=self.data["time_echo"]
                else:
                    if 'time_echo' in value["dim"]:
                        timeval=self.data["time_echo"]
                        prior_qa=np.zeros(shape_var)
                    elif 'time_vel' in value["dim"]:
                        timeval=self.data["time_vel"]
                        prior_qa=np.zeros(shape_var)
                    else:
                        raise Exception("Wrong dimension")
                if simple: # Simple quality check only
                    self.data[name] = prior_qa+qualityassurance(data_array, timeval, **quality_assurance_dict[key]["simple"])
                else:
                    quality_assurance_all = dict(quality_assurance_dict[key]["simple"], **quality_assurance_dict[key]["advanced"])
                    self.data[name] = prior_qa+qualityassurance(data_array, timeval, **quality_assurance_all)
               


    def derive_variables(self, rotate_velocity):
        log("Computing derived variables.", indent=1)
        self.variables.update(self.derived_variables)
        
        # Add additional variables here
        
    def add_metadata(self, meta_dict):
        """
        Add the metadata as general attributes.

        Parameters:
            meta_dict (dictionary): metadata from json file.
                
        Returns: nothing
        """
        for key, value in meta_dict.items():
            if isinstance(value, str): # Add it directly
                self.general_attributes[key]=value
            elif isinstance(value, dict): # Add each element of the dictionary
                for key2, value2 in value.items():
                    if isinstance(value2, str): # Add it directly
                        self.general_attributes[key+"_"+key2]=value2
                    else:
                        try:
                            self.general_attributes[key+"_"+key2]=str(value2)
                        except:
                            continue
            else:
                try:
                    self.general_attributes[key]=str(value)
                except:
                    continue
                                
                        
        