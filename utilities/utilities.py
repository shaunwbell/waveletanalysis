#!/usr/bin/env 
"""
 utilities.py
 

   Using Anaconda packaged Python
"""
import datetime
import matplotlib as mpl
import numpy as np

__author__   = 'Shaun Bell'
__email__    = 'shaun.bell@noaa.gov'
__created__  = datetime.datetime(2013, 12, 20)
__modified__ = datetime.datetime(2013, 12, 20)
__version__  = "0.1.0"
__status__   = "Development"

"""---------------------------------------------------------------------------------"""

def EPICdate2udunits(time1, time2):
    """
    Inputs
    ------
          time1: array_like
                 True Julian day
          time2: array_like
                 Milliseconds from 0000 GMT
    Returns
    -------
          dictionary:
            'timeint': python serial time
            'interval_min': data interval in minutes
    
    Example
    -------
    Python uses days since 0001-01-01 and a gregorian calendar

      
    Reference
    ---------
    PMEL-EPIC Conventions (misprint) says 2400000
    http://www.epic.noaa.gov/epic/eps-manual/epslib_ch5.html#SEC57 says:
    May 23 1968 is 2440000 and July4, 1994 is 2449538
              
    """
    ref_time_py = datetime.datetime.toordinal(datetime.datetime(1968, 5, 23))
    ref_time_epic = 2440000
    
    offset = ref_time_epic - ref_time_py
    
    pytime = [None] * len(time1)
    pytimestr = [None] * len(time1)
    
    for i, val in enumerate(time1):
        pyday = time1[i] - offset 
        pyfrac = time2[i] / (1000. * 60. * 60.* 24.) #milliseconds in a day
        
        pytime[i] = (pyday + pyfrac)
        
    #number of datapoints to average to get a daily average
    interval_min = np.ceil((60 * 60 * 24) / (np.ceil((pytime[1]-pytime[0]) * 60 * 60 * 24))) 
    
    return {'timeint':pytime, 'interval_min':interval_min}
    
def moving_average(x, n, type='simple'):
    """
    compute an n period moving average.

    type is 'simple' | 'exponential'

    """
    #x = np.asarray(x)
    if type=='simple':
        weights = np.ones(n)
    else:
        weights = np.exp(np.linspace(-1., 0., n))

    weights /= weights.sum()


    a =  np.convolve(x, weights, mode='full')[:len(x)]
    a[:n] = a[n]
    return a
    
def subsample(data, sample_size):
    samples = data[0::sample_size]
    return samples
"""---------------------------------------------------------------------------------"""

def main():
    """ TODO """
    
if __name__ == "__main__":
    main() 