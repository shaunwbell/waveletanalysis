#!/usr/bin/env
"""
 constants.py
 

   Using Anaconda packaged Python
"""
import os, csv, pickle
import datetime

__author__   = 'Shaun Bell'
__email__    = 'shaun.bell@noaa.gov'
__created__  = datetime.datetime(2013, 12, 18)
__modified__ = datetime.datetime(2013, 12, 18)
__version__  = "0.1.0"
__status__   = "Development"



def nc_global_attr():
    """
    Details:
    --------
        List of mandatory global attribute names for NetCDF Files
        
    Path:
    -----
        Found in EPICNetCDF folder
        
    """
    gattr =    ['CREATION_DATE',
                'CRUISE',
                'CAST',
                'INST_TYPE',
                'DATA_TYPE',
                'WATER_MASS',
                'BAROMETER',
                'WIND_DIR',
                'WIND_SPEED',
                'AIR_TEMP',
                'WATER_DEPTH',
                'STATION_NAME']
    
    return gattr

def nc_global_attr_moor():
    """
    Details:
    --------
        List of mandatory global attribute names for NetCDF Files
        
    Path:
    -----
        Found in EPICNetCDF folder
        
    """
    gattr =    ['CREATION_DATE',
                'MOORING',
                'PROJECT',
                'WATER_DEPTH']
    
    return gattr
        
def nc_vars_ctd():
    """
    Details
    -------
        List of mandatory variables to be in the NetCDF File.
        Provided as EPIC-Codes
        
    Path:
    -----
        Found in EPICNetCDF folder
        
    """
    nc_ctd_variables =  ['time',
                         'time2',
                         'dep',
                         'lat',
                         'lon',
                         'T_28',
                         'T2_35',
                         'PAR_916',
                         'S_41',
                         'S_42']
    
    return nc_ctd_variables
                         
                         
def nc_vars_moor():
    """
    Details
    -------
        List of mandatory variables to be in the NetCDF File.
        Provided as EPIC-Codes
        
    Path:
    -----
        Found in EPICNetCDF folder
        
    """
    nc_moor_variables =  ['time',
                         'time2',
                         'dep',
                         'depth',
                         'lat',
                         'latitude',
                         'lon',
                         'longitude',
                         'T_20',
                         'P_1',
                         'S_41',
                         'C_50',
                         'U_320',
                         'V_321']
    return nc_moor_variables
    
    
def mooring_codes():
    """
    Details
    -------
        Using mooring_codes.txt - ingest information for future use
        
    Path:
    -----
        Found in EPICNetCDF folder
        
    """    
    d = {}
    path = os.path.abspath(__file__)
    dir_path = os.path.dirname(path)

    with open(dir_path + '/data/mooring_codes.txt', 'rb') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        for row in csv_reader:
            print row
            d[row[0].strip()] = row[1:]
            
    pickle.dump( d, open( dir_path + '/data/mooring_codes.p', "wb" ) )

def read_mooring_codes():
    path = os.path.abspath(__file__)
    dir_path = os.path.dirname(path)

    if os.path.isfile(dir_path + '/data/mooring_codes.p'):
        moor_codes = pickle.load( open( dir_path + '/data/mooring_codes.p', "rb" ) )
        return moor_codes
    else:
        sys.exit("No success in creating pickle file!")

def read_plotmooring_meta(path_to_metacsv):
    depth = []
    lat = []
    lon = []
    year = []
    mooring = []
    mooringID = []
    with open(path_to_metacsv + "latlon.csv", 'rb') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=",")
        next(csv_reader) #skip header
        for row in csv_reader:
            r0,r1,r2,r3,r4,r5 = row
            depth.append(r0)
            lat.append(r1)
            lon.append(r2)
            year.append(r3)
            mooring.append(r4)
            mooringID.append(r5)
    return (depth,lat,lon,year,mooring, mooringID) 
"""---------------------------------------------------------------------------------"""

def main():
    mooring_codes()
    
if __name__ == "__main__":
    main()