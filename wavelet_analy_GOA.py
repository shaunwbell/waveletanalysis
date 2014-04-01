#!/usr/bin/env 
"""
 Based on waveletanalysis.py
 
 For P. Stabeno
 

   Using Anaconda packaged Python
   modifications for confidence intervals based on wave_matlab at
   http://paos.colorado.edu/research/wavelets/
   
### License ###

The MIT License (MIT)

Copyright (c) 2013 Aaron O'Leary (dev@aaren.me)

Permission is hereby granted, free of charge, to any person obtaining a copy of this 
software and associated documentation files (the "Software"), to deal in the Software 
without restriction, including without limitation the rights to use, copy, modify, merge, 
publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons 
to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies 
or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE
OR OTHER DEALINGS IN THE SOFTWARE.
"""

#Standard packages
import os, datetime

#Science packages
import numpy as np
from netCDF4 import Dataset
from scipy import signal

#User Packages
import wave_sig
from general_utilities.wavelets_bams.wavelets import WaveletAnalysis
from utilities import ncutilities as ncutil
from utilities import utilities as util
from utilities import constants 
import wavelet_analy_plot

__author__   = 'Shaun Bell'
__email__    = 'shaun.bell@noaa.gov'
__created__  = datetime.datetime(2014, 01, 13)
__modified__ = datetime.datetime(2014, 01, 13)
__version__  = "0.1.0"
__status__   = "Development"

def from_netcdf(infile, parameter):
    """ Uses ncreadfile_dic which returns a dictionary of all data from netcdf"""
    ncfile = infile

    ###nc readin/out
    nchandle = ncutil.ncopen(ncfile)
    params = ncutil.get_vars(nchandle) #gets all of them
    if not parameter in params:
        ### specified parameter doesn't exist
        return ('no param')
    ncdata = ncutil.ncreadfile_dic(nchandle, params)
    ncutil.ncclose(nchandle)

    ###data massaging
    pytime = util.EPICdate2udunits(ncdata['time'], ncdata['time2'])
    

    dt = 1. / pytime['interval_min'] #data is 4 times daily
    time_base = 'days'
    time = pytime['timeint']
    #time = util.subsample(time, int(pytime.get('interval_min')) / 4)

    
    # signal
    data = ncdata[parameter][:,0,0,0]
    variance = np.var(data)
    
    #normalize
    print 'Variance = %s ' % (variance)
    x = (data - np.mean(data)) / np.sqrt(variance)
    variance = np.var(x)

    return_data = {'data': data, 'anom': x, 'dt': dt, 'time':np.array(time), 'variance':variance, 'time_base':time_base}
    return (return_data)
    
"""---------------------------------other routines-------------------------------------"""

    
def walk_dir(search_dir):
    """Lists of files, paths, dir"""
    allroot = []
    allfile = []
    for root, dirs, files in os.walk(search_dir):
        for file in files:
            if ('_' in file) or ('_' in file):
                allroot.append(os.path.join(root, file))
                allfile.append(file) 
                
    return (allroot, allfile)
    
"""---------------------------------modules--------------------------------------------"""

def acf(series):
    """Determine autocorrelation factors"""
    n = len(series)
    data = np.asarray(series)
    mean = np.mean(data)
    c0 = np.sum((data - mean) ** 2) / float(n)

    def r(h):
        acf_lag = ((data[:n - h] - mean) * (data[h:] - mean)).sum() / float(n) / c0
        return round(acf_lag, 3)
    x = np.arange(n) # Avoiding lag 0 calculation
    acf_coeffs = map(r, x)
    return acf_coeffs  
    
"""------------------------------ Data Selection --------------------------------------"""

#for year in (2001, 2002, 2003, 2004, 2005):
for year in (2002, 2003):

    ### loop through every rcm an7/an9 file in a dir
    stationID = 'globec2'
    (allroot, allfile) = walk_dir('/Users/bell/Data_Local/FOCI/Mooring/%s/%s/' % (year, stationID) )


    for dfile_1 in allroot:
        print dfile_1
        par_1 = ['S_41', 'salinity']
        (data_1) = from_netcdf(dfile_1, par_1[0])

        
        #alternative if secondary param doesn't exist
        if (data_1 == 'no param'):
            continue

        fig_name_base = 'images/' + dfile_1.split('/')[-1].split('.')[0] + '_' + par_1[0]
    
        depth = dfile_1.split('m.nc')[0][-4:]



        """-----------------------------wavelet analysis           ---------------------------"""

        wa = WaveletAnalysis(data_1['anom'], time=data_1['time'], dt=data_1['dt'], dj=0.125)

        # wavelet power spectrum
        power = wa.wavelet_power
        transform = wa.wavelet_transform

        # scales 
        scales = wa.scales

        # associated time vector
        t = wa.time

        # reconstruction of the original data
        rx = wa.reconstruction()

        # determine acor factor for red noise
        acorr = acf(data_1['anom'])
        lagf = (acorr[1]+np.sqrt(acorr[2]))/2
        print 'acorr lagf is %s' % lagf

        # determine significance levels
        (signif, fft_theory) = wave_sig.wave_signif(data_1['anom'],data_1['dt'],scales,lag1=lagf)
        sig95 = np.ones_like(power) * np.array([signif] * len(t)).transpose()
        sig95 = power / sig95         # where ratio > 1, power is significant

        # Global wavelet spectrum & significance levels:
        global_int = data_1['variance']*(np.sum(power, axis=0) ) / data_1['anom'].shape[0]   # time-average over all times
        gs = ((np.sum(power, axis=1) ) / data_1['anom'].shape[0]) / data_1['variance'] #assume var=1
        gswa = wa.global_wavelet_spectrum
        # Global wavelet significance
        (signif_g, fft_theory_g) = wave_sig.global_wave_signif(data_1['anom'],data_1['dt'],scales,lag1=lagf,sigtest=1, dof=len(data_1['anom']))




        """----------------------------- plot setup ------------------------------------------"""
        T, S = np.meshgrid(t, scales)


        """----------- plotting WaveTransform Power with confidence interval contour ----------"""



        """----------------- plotting contours w/global and timeseries ----------"""
        """----------------- zoom in to specified scales ----------"""

        plt, fig = wavelet_analy_plot.plot_wavetransf_time_zoom(data_1['anom'], wa, T, S, sig95, gs, signif_g, 
            data_1['time_base'], scalemin=.1, scalemax=64, ylabel=par_1[1], plot_percentile=True)
        plt.savefig(fig_name_base + '_wave3' + depth + '.png', bbox_inches='tight', dpi = (100))
        plt.close()

        """----------------------- plotting power spectrum FFT --------------------------------"""
        (plt, fig) = wavelet_analy_plot.fft_power_spec(data_1['anom'], data_1['time_base'], Fs=24)

        plt.savefig(fig_name_base + '_FFTspec' + depth + '.png', bbox_inches='tight', dpi = (100))
        plt.close()





