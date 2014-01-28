#!/usr/bin/env 
"""
 ncutilities.py
 

   Using Anaconda packaged Python
"""
import datetime
import numpy as np
from netCDF4 import Dataset 
#user defined
import constants as constants

__author__   = 'Shaun Bell'
__email__    = 'shaun.bell@noaa.gov'
__created__  = datetime.datetime(2013, 12, 20)
__modified__ = datetime.datetime(2013, 12, 20)
__version__  = "0.1.0"
__status__   = "Development"

"""---------------------------------------------------------------------------------"""

def ncopen(ncfile):
    """
    Parameters
    ----------
    TODO
    
    Returns
    -------
    TODO
              
    """
    nchandle = Dataset(ncfile,'r')
    return nchandle
    
def ncclose(nchandle):
    """
    Parameters
    ----------
    TODO
    
    Returns
    -------
    TODO
              
    """
    nchandle.close()

def get_global_atts(nchandle):
    att_names = nchandle.ncattrs()
    return att_names
    
def ncreadfile(nchandle, params):

    data = np.zeros( (nchandle.variables['time'][:].shape[0],len(params)) ) 
    for j, v in enumerate(params): 
        if v in nchandle.variables.keys(): #check for nc variable
            if (v in params[:8]): #treat one-D vars differently - time, time2 depth lat lon
                data[:,j] = nchandle.variables[v][:]
            else:
                data[:,j] = nchandle.variables[v][:,0,0,0]
        else: #if parameter doesn't exist fill the array with zeros
            data[:,j] = data[:,j-1] * 0.0
    return (data)
    

def check_mand_atts(nchandle):
    
    att_names = get_global_atts(nchandle)
    nc_global_attr = constants.nc_global_attr()
    for val in nc_global_attr[1]:
        if not val in att_names:
            return False
        else:
            return True
            
def ncmeta(nchandle, metaparams):
    """get variables from netcdf file, gets all params in list sent eg. ['time','time2',...]
    Returns a dictionary of values:
    cruisePress: barometric reading during cast
    cruiseID:   Cruise Identifier
    castID:     cast number
    ncdata:     array-like list of values"""
    
    if not check_mand_atts(nchandle):
        mooring = 'missing'
        project = 'missing'
        waterd = 'missing'
    else:
        mooring = (nchandle.MOORING).strip('"')
        project = nchandle.PROJECT
        waterd = nchandle.WATER_DEPTH
    
    return {'Mooring':mooring, 'Project':project, 'Water_Depth':waterd}
    
def nc_missing(data, flag=1e35, setfill='Mean'):
    mx = np.ma.masked_values(data, flag)
        
    if setfill == 'Mean':
        return mx.filled(np.ma.mean(mx))
    elif setfill == 'Zero':
        return mx.filled(0)
    else:
        return (mx)
"""---------------------------------------------------------------------------------"""

def main():
    """ TODO """
    
if __name__ == "__main__":
    main() 