# -*- coding: utf-8 -*-
import os
import sys
import netCDF4
import numpy as np
from datetime import datetime, timezone
import gpxpy
sys.path.append(os.path.join(os.path.dirname(__file__), '../../functions'))
from general_functions import GenericInstrument


class gps(GenericInstrument):
    def __init__(self, *args, **kwargs):
        super(gps, self).__init__(*args, **kwargs)
        self.general_attributes = {
            "institution": "Unil",
            "source": "Garmin GPS",
            "references": "",
            "history": "See history on Renku",
            "conventions": "CF 1.7",
            "comment": "GPS track of the Lém boat and floating Signature ADCP",
            "title": "GPS track",
            "Instrument": "Garmin GPS eTrex 20x",
        }
        
        self.dimensions = {
            'time': {'dim_name': 'time', 'dim_size': None},
        }
        
        self.variables = {
            'time': {'var_name': 'time', 'dim': ('time',), 'unit': 'seconds since 1970-01-01 00:00:00', 'long_name': 'Time of the GPS data'},
            'lon': {'var_name': 'lon', 'dim': ('time',), 'unit': 'degrees East', 'long_name': 'Longitude'},
            'lat': {'var_name': 'lat', 'dim': ('time',), 'unit': 'degrees North', 'long_name': 'Latitude'},
        }

        self.derived_variables = {
        }

        self.data = {}

    def read_data(self, file):
        """
        Read the GPS data and store it in a gps object.

        Parameters:
            file (str): path and GPS filename (.gpx file)
                           
        Returns:
            True if the data was correctly read, False otherwise
        """
        
        print("Parsing data from {}.".format(file))
        
        try:
            # Open the .gpx file
            with open(file, 'r') as gpx_file:
                gpx = gpxpy.parse(gpx_file)
    
            # Access data
            GPS_points=gpx.tracks[0].segments[0].points
            dateGPS=np.array([pt.time for pt in GPS_points]) # UTC
            self.data["time"]=np.array([dt.replace(tzinfo=timezone.utc).timestamp() for dt in dateGPS])
            self.data["lon"]=np.array([pt.longitude for pt in GPS_points])
            self.data["lat"]=np.array([pt.latitude for pt in GPS_points])
            return True
    
        except:
            print("Failed to process {}.".format(file))
            return False
 
        
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
                                
                        
        