#!/usr/bin/env 
"""
 textreadutilities.py
 

   Using Anaconda packaged Python
"""
import datetime
import numpy as np


__author__   = 'Shaun Bell'
__email__    = 'shaun.bell@noaa.gov'
__created__  = datetime.datetime(2014, 01, 02)
__modified__ = datetime.datetime(2014, 01, 02)
__version__  = "0.1.0"
__status__   = "Development"

"""---------------------------------------------------------------------------------"""

def read_txt(infile):
    """ """
    

    season, year, temp, anom = np.loadtxt(infile, dtype=[('f0','|S4'),('f1',int),('f2',float),('f3',float)], skiprows=1, unpack=True)

    return (season, year, temp, anom)
    
