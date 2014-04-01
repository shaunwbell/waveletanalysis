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
from netCDF4 import Dataset


#User Packages
import wave_sig
from general_utilities.wavelets_bams.wavelets import WaveletAnalysis
from utilities import ncutilities as ncutil
from utilities import utilities as util
from utilities import constants 
import wavelet_analy_plot
import wave_xwt

__author__   = 'Shaun Bell'
__email__    = 'shaun.bell@noaa.gov'
__created__  = datetime.datetime(2014, 03, 04)
__modified__ = datetime.datetime(2014, 03, 04)
__version__  = "0.1.0"
__status__   = "Development"

"""---------------------------Read in Data File----------------------------------------"""

def example(infile):
    """ 
    Data file from "Crosswavelet and wavelet coherence software were provided by
   A. Grinsted." (C) Aslak Grinsted 2002-2004
    """
    ###text readin (two columns, column 1 is date, column 2 is value)
    dir_path = os.path.dirname(os.path.abspath(__file__))
    full_path = dir_path +infile 

    data = np.loadtxt(full_path, unpack=True)

    ###data massaging
    dt = 1. #data is yearly
    time_base = 'years'
    time = data[0]
    # convert to date ordinal if necessary or activate autoxlabel
    
    # signal
    x = data[1]
    variance = np.var(x)
    
    #normalize
    print 'Variance = %s ' % (variance)
    x = (x - np.mean(x)) / np.sqrt(variance)
    variance = np.var(x)

    return_data = {'data': data[1], 'anom': x, 'dt': dt, 'time':np.array(time), 'variance':variance, 'time_base':time_base}
    return (return_data)
    
"""---------------------------------other routines-------------------------------------"""

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


### ingest data
(data_1) = example('/data/jao.txt')
(data_2) = example('/data/jbaltic.txt') #reconfigured pdf in example, not done here

"""-------------------------- Find and Match Common intervals -------------------------"""
### find and match common time interval
lbound =  np.max([data_1['time'].min(), data_2['time'].min()])
ubound = np.min([data_1['time'].max(), data_2['time'].max()])
wa2_ind = ((data_2['time'] <= ubound ) & (data_2['time'] >= lbound)).nonzero()[0] 
wa1_ind = ((data_1['time'] <= ubound ) & (data_1['time'] >= lbound)).nonzero()[0] 

data_1c, data_2c = {}, {}
data_2c['time'] = data_2['time'][wa2_ind,:]
data_1c['time'] = data_1['time'][wa1_ind,:]
data_2c['anom'] = data_2['anom'][wa2_ind,:]
data_1c['anom'] = data_1['anom'][wa1_ind,:]

"""-----------------------------wavelet analysis           ---------------------------"""

wa1c = WaveletAnalysis(data_1c['anom'], time=data_1c['time'], dt=data_1['dt'], dj=0.125)
wa2c = WaveletAnalysis(data_2c['anom'], time=data_2c['time'], dt=data_2['dt'], dj=0.125)


# determine acor factor for red noise stream 1
acorr_1 = acf(data_1c['anom'])
lagf_1 = (acorr_1[1]) #simpler lag
print 'acorr lagf for datastream 1 is %s' % lagf_1

# determine acor factor for red noise stream 2
acorr_2 = acf(data_2c['anom'])
lagf_2 = (acorr_2[1]) #simpler lag
print 'acorr lagf for datastream 2 is %s' % lagf_2

# determine significance levels for stream 1
(signif_1, fft_theory_1) = wave_sig.wave_signif(data_1c['anom'],wa1c.dt,wa1c.scales,lag1=lagf_1)
sig95_1 = np.ones_like(wa1c.wavelet_power) * np.array([signif_1] * len(wa1c.time)).transpose()
sig95_1 = wa1c.wavelet_power / sig95_1         # where ratio > 1, power is significant

# determine significance levels for stream 1
(signif_2, fft_theory_2) = wave_sig.wave_signif(data_2c['anom'],wa2c.dt,wa2c.scales,lag1=lagf_2)
sig95_2 = np.ones_like(wa2c.wavelet_power) * np.array([signif_2] * len(wa2c.time)).transpose()
sig95_2 = wa2c.wavelet_power / sig95_2         # where ratio > 1, power is significant

# Global wavelet spectrum & significance levels stream 1:
global_int_1 = data_1['variance']*(np.sum(wa1c.wavelet_power, axis=0) ) / data_1c['anom'].shape[0]   # time-average over all times
gs_1 = ((np.sum(wa1c.wavelet_power, axis=1) ) / data_1c['anom'].shape[0]) / data_1['variance'] #assume var=1
gswa_1 = wa1c.global_wavelet_spectrum
# Global wavelet significance
(signif_g_1, fft_theory_g_1) = wave_sig.global_wave_signif(data_1c['anom'],wa1c.dt,wa1c.scales,lag1=lagf_1,sigtest=1, dof=len(data_1c['anom']))

# Global wavelet spectrum & significance levels stream 2:
global_int_2 = data_2['variance']*(np.sum(wa2c.wavelet_power, axis=0) ) / data_2c['anom'].shape[0]   # time-average over all times
gs_2 = ((np.sum(wa2c.wavelet_power, axis=1) ) / data_2c['anom'].shape[0]) / data_2['variance'] #assume var=1
gswa_2 = wa2c.global_wavelet_spectrum
# Global wavelet significance
(signif_g_2, fft_theory_g_2) = wave_sig.global_wave_signif(data_2c['anom'],wa2c.dt,wa2c.scales,lag1=lagf_2,sigtest=1, dof=len(data_2c['anom']))

"""----------------------------- plot setup ------------------------------------------"""
T1, S1 = np.meshgrid(data_1c['time'], wa1c.scales)
T2, S2 = np.meshgrid(data_2c['time'], wa2c.scales)


"""----------- plotting WaveTransform Power with confidence interval contour ----------"""


fig_name_base = 'images/' + 'xwt_example'


"""----------------- plotting contours w/global and timeseries ----------"""
"""----------------- zoom in to specified scales ----------"""

plt, fig = wavelet_analy_plot.plot_wavetransf_time_zoom(data_1c['anom'], wa1c, T1, S1, sig95_1,\
 gs_1, signif_g_1, data_1['time_base'], scalemin=2, scalemax=64, ylabel='jao', plot_percentile=True)
plt.savefig(fig_name_base + 'jao' + '.png', bbox_inches='tight', dpi = (100))
plt.close()

plt, fig = wavelet_analy_plot.plot_wavetransf_time_zoom(data_2c['anom'], wa2c, T2, S2, sig95_2,\
 gs_2, signif_g_2, data_2['time_base'], scalemin=2, scalemax=64, ylabel='jbaltic', plot_percentile=True)
plt.savefig(fig_name_base + 'jbaltic' + '.png', bbox_inches='tight', dpi = (100))
plt.close()    




"""-----------------------------cross wavelet analysis      ---------------------------"""
### since the time scales are matched and the dt is the same, the wavelet analysis for
#   dataset 1 or two can be used and they should be equivalent for coi, scales and fourier
#   period.

# if data periods are significantly different between dataset, more care should be taken when
# truncating their times and utilizing their time windows.

(xwt, axwt)= wave_xwt.xwt(wa1c,wa2c)

# xwt significance
### values come from 95% CI for sqrt of the product of two chi2 distributions
V=2
Zv=3.9999

P1 = ar1_powspec(lagf_1, wa1c.fourier_periods / wa1c.dt)
P2 = ar1_powspec(lagf_2, wa2c.fourier_periods / wa1c.dt)
sig95 = data_1['variance'] * data_2['variance'] * np.sqrt(P1 * P2) * Zv / V
sig95 = np.ones_like(wa2c.wavelet_power) * np.array([sig95] * len(wa2c.time)).transpose()
sig95 = np.abs(xwt) / sig95

(plt, fig) = wavelet_analy_plot.plot_xwt_wavetransf(np.abs(xwt), data_2c['time'], \
    wa2c, T2, S2, sig95, axwt, data_2['time_base'], scalemin=2, scalemax=64, plot_percentile=True)
plt.savefig(fig_name_base + 'jao_jbaltic_' + 'xwt.png', bbox_inches='tight', dpi = (100))
plt.close()    
    
"""-----------------------------cross wavelet coherence      ---------------------------"""

(sw1, sw2, sw_xwt) = wave_xwt.wtc(wa1c,wa2c, xwt)
    
Rsq = (np.abs((sw_xwt)**2) / (sw1 * sw2))
whereAreNaNs = np.isnan(Rsq);
Rsq[whereAreNaNs] = 0;
            
sig95 = sig95 * 0 

(plt, fig) = wavelet_analy_plot.plot_xwt_wavetransf(Rsq, data_2c['time'], \
    wa2c, T2, S2, sig95, axwt, data_2['time_base'], scalemin=2, scalemax=64, plot_percentile=False)
plt.savefig(fig_name_base + 'jao_jbaltic_' + 'coherence.png', bbox_inches='tight', dpi = (100))
plt.close()  