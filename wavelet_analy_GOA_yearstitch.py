#!/usr/bin/env 
"""
 Wavelt Coherence and X-wavelet transform
 
 For P. Stabeno
 

   Using Anaconda packaged Python
   
   Based on
   --------
   
   xwt from -   "Crosswavelet and wavelet coherence software were provided by
   A. Grinsted." (C) Aslak Grinsted 2002-2004
   
   modifications for confidence intervals based on wave_matlab at
   http://paos.colorado.edu/research/wavelets/
   
   wavelet analysis from https://github.com/aaren/wavelets
   
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
import scipy as sp

#User Packages
import wave_sig
from general_utilities.wavelets_bams.wavelets import WaveletAnalysis
from general_utilities import lanzcos
from utilities import ncutilities as ncutil
from utilities import utilities as util
import wavelet_analy_plot_goa

__author__   = 'Shaun Bell'
__email__    = 'shaun.bell@noaa.gov'
__created__  = datetime.datetime(2014, 03, 04)
__modified__ = datetime.datetime(2014, 07, 23)
__version__  = "0.1.0"
__status__   = "Development"
__keywordss__   = 'wavelet', 'GOA', 'stitch moorings', 'publications'

"""---------------------------Read in Data File----------------------------------------"""


def from_netcdf(infile, parameter, detrend=False):
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
    
    #test routine, set default to false
    if detrend: #detrend if flag is true - default false
        x = sp.signal.detrend(data, axis=0)
    else:
        x = data
        
    #normalize
    print 'Variance = %s ' % (variance)
    x = (x - np.mean(x)) / np.sqrt(variance)
    variance = np.var(x)

    return_data = {'data': data, 'anom': x, 'dt': dt, 'time':np.array(time), 'variance':variance, 'time_base':time_base}
    return (return_data)    

def from_2netcdf(infile, parameter):
    """ For moorings with a, and b deployments
    Ingest and concatenate both files
    """
    ncfile = infile[0]

    ###nc readin/out
    nchandle = ncutil.ncopen(ncfile)
    params = ncutil.get_vars(nchandle) #gets all of them
    if not parameter in params:
        ### specified parameter doesn't exist
        return ('no param')
    ncdata1 = ncutil.ncreadfile_dic(nchandle, params)
    ncutil.ncclose(nchandle)

    ncfile = infile[1]

    ###nc readin/out
    nchandle = ncutil.ncopen(ncfile)
    params = ncutil.get_vars(nchandle) #gets all of them
    if not parameter in params:
        ### specified parameter doesn't exist
        return ('no param')
    ncdata2 = ncutil.ncreadfile_dic(nchandle, params)
    ncutil.ncclose(nchandle)


    ###data massaging
    pytime1 = util.EPICdate2udunits(ncdata1['time'], ncdata1['time2'])
    pytime2 = util.EPICdate2udunits(ncdata2['time'], ncdata2['time2'])
    

    dt = 1. / pytime1['interval_min'] #data is 4 times daily
    time_base = 'days'
    ### fill time gap
    blanktime = np.array([])
    for tg in np.arange(np.max(pytime1['timeint'])+dt,np.min(pytime2['timeint']),dt):
        blanktime = np.hstack((blanktime,tg))

    time = np.hstack((pytime1['timeint'], blanktime))
    time = np.hstack((time,pytime2['timeint']))
    
    
    # signal
    data_interp = np.interp(blanktime,[np.max(pytime1['timeint']),np.min(pytime2['timeint'])],[ncdata1[parameter][-1,0,0,0], ncdata2[parameter][0,0,0,0]])

    data = np.hstack((ncdata1[parameter][:,0,0,0], data_interp))
    data = np.hstack((data, ncdata2[parameter][:,0,0,0]))

    variance = np.var(data)
    
    #normalize
    print 'Variance = %s ' % (variance)
    x = (data - np.mean(data)) / np.sqrt(variance)
    variance = np.var(x)

    return_data = {'data': data, 'anom': x, 'dt': dt, 'time':np.array(time), 'variance':variance, 'time_base':time_base}
    return (return_data)
"""---------------------------------stats routines-------------------------------------"""

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

def ar1_powspec(ar1,period):
    """
    From Grinsted ar1spectum.m
     %http://www.madsci.org/posts/archives/may97/864012045.Eg.r.html
     %fixed typo in numerical recipes
    """
    
    freq = 1. / period
    P = (1. - ar1 ** 2.) / (np.abs(1. - ar1 * np.exp(-2. * np.pi * freq * 1j))) ** 2. #complex
    
    return P          
    
"""-------------------------- Analysis Routine ----------------------------------------"""
"""-------------------------- ---------------- ----------------------------------------"""

for stationID in ('globec1',): #for multiple stations
    for year in (2001,): #for multiple years
        file1 = '/Users/bell/Data_Local/FOCI/Mooring/%s/%s/01gb1a_sc_0031m.nc' % (year, stationID)
        file2 = '/Users/bell/Data_Local/FOCI/Mooring/%s/%s/01gb1b_sc_0017m.nc' % (year, stationID)

        par_1 = 'S_41'
        data_1 = from_2netcdf([file1,file2], par_1)
        #data_1 = from_netcdf(file1, par_1)
        
        
        """-----------------------------wavelet analysis           ---------------------------"""

        """35hr lanzcos filter data - watch boundaries"""
        filter_data = True
        if filter_data:
            data_1c = {}
            data_1c['anom'] = lanzcos.lanzcos35(data_1['anom'],data_1['dt'],Cf=35.)
            data_1c['time'] = data_1['time']
        else:
            data_1c = data_1
    
        wa1c = WaveletAnalysis(data_1c['anom'], time=data_1c['time'], dt=data_1['dt'], dj=0.125)


        # determine acor factor for red noise stream 1
        # uses entire dataperiod for corr 
        
        #TODO: current implementation uses any preprocessed changes in ingest module (detrended, standardized, etc)
        acorr_1 = acf(data_1['anom'])
        lagf_1 = (acorr_1[1]+np.sqrt(acorr_1[2]))/2
        print 'acorr lagf for datastream 1 is %s' % lagf_1

        scalemin = 1
        scalemax = 32
        scale_ind = ((wa1c.scales >= scalemin) & (wa1c.scales <= scalemax))


        # determine significance levels for stream 1
        (signif_1, fft_theory_1) = wave_sig.wave_signif(data_1c['anom'],wa1c.dt,wa1c.scales[scale_ind],lag1=lagf_1)
        sig95_1 = np.ones_like(wa1c.wavelet_power[scale_ind,:]) * np.array([signif_1] * len(wa1c.time)).transpose()
        sig95_1 = wa1c.wavelet_power[scale_ind,:] / sig95_1         # where ratio > 1, power is significant


        # Global wavelet spectrum & significance levels stream 1:
        global_int_1 = data_1['variance']*(np.sum(wa1c.wavelet_power, axis=0) ) / data_1c['anom'].shape[0]   # time-average over all times
        gs_1 = ((np.sum(wa1c.wavelet_power, axis=1) ) / data_1c['anom'].shape[0]) / data_1['variance'] #assume var=1
        gswa_1 = wa1c.global_wavelet_spectrum
        # Global wavelet significance
        (signif_g_1, fft_theory_g_1) = wave_sig.global_wave_signif(data_1c['anom'],wa1c.dt,wa1c.scales,lag1=lagf_1,sigtest=1, dof=len(data_1c['anom']))


        """----------------------------- plot setup ------------------------------------------"""


        """----------- plotting WaveTransform Power with confidence interval contour ----------"""


        fig_name_base = 'images/' + file1.split('/')[-1].split('.')[0] + '_' + "_".join(file2.split('/')[-1].split('_')[1:3]).split('.')[0] + '_'


        """----------------- zoom in to specified scales ----------"""
        T1, S1 = np.meshgrid(data_1c['time'], wa1c.scales[scale_ind])
        plt, fig = wavelet_analy_plot_goa.plot_wavetransf_time_zoom(data_1c['anom'], wa1c, T1, S1, sig95_1,\
         gs_1, signif_g_1, data_1['time_base'], data_ind=scale_ind, scalemin=scalemin, scalemax=scalemax, ylabel=par_1[1], plot_percentile=True)
        plt.savefig(fig_name_base + par_1  + '.pdf', bbox_inches='tight', dpi = (300)) #conver to eps with imagemagick
        plt.close()


