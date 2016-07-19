#!/usr/bin/env 
"""
 Based on waveletanalysis.py
 
 For Adam
 
 One station, One Year - many depth levels
 

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
import datetime

#Science packages
import numpy as np

#User Packages
import wave_sig
from general_utilities.wavelets_bams.wavelets import WaveletAnalysis
import wavelet_analy_plot_hourly as wavelet_analy_plot

__author__   = 'Shaun Bell'
__email__    = 'shaun.bell@noaa.gov'
__created__  = datetime.datetime(2014, 01, 13)
__modified__ = datetime.datetime(2014, 04, 16)
__version__  = "0.1.0"
__status__   = "Development"

"""----------------------methods for reading raw data----------------------------------"""

def ADCP_2D(datafilein, timefilein, depthfilein,level, matlab_timeoffset=366.):
    """ 
    datafile used yyyy_BV.txt - n rows, m columns for time/depth
    timefile used yyyy_date.txt - n rows, 1 columns for time
    depthfile used is yyyy_depth.txt - nrows, 1 column for depth
    TAPS dataset from Adam
    """

    data = np.loadtxt(datafilein, unpack=True)
    time = np.loadtxt(timefilein, unpack=True) - matlab_timeoffset
    deptharray = np.loadtxt(depthfilein, unpack=True)

    ###data massaging
    dt = 1. #data is hourly
    time_base = 'hours'
    depth = -1. * deptharray

    # signal
    xx = data[:,level]

    #find and fill nan with smallest value
    data_nan_ind = np.where([np.isnan(xx)])[1]
    xx[data_nan_ind] = np.nanmin(xx)

    variance = np.var(xx)


    #normalize
    print 'Variance = %s ' % (variance)
    x = (xx - np.mean(xx)) / np.sqrt(variance)
    variance = np.var(x)

    sl = len(x)
    #time = (np.arange(0,sl ,1) * dt/24. ) + startdate #date_object = datetime.datetime.strptime('04 23 2009', '%m %d %Y')

    return (data, x,dt,np.array(time) * 24., variance, time_base, depth) 

def ADCP_2D_subset(datafilein, timefilein, depthfilein, t_min, t_max, level, matlab_timeoffset=366.):
    """ 
    datafile used yyyy_BV.txt - n rows, m columns for time/depth
    timefile used yyyy_date.txt - n rows, 1 columns for time
    depthfile used is yyyy_depth.txt - nrows, 1 column for depth
    TAPS dataset from Adam
    Specify bounding times as python ordinaltimes
    """
    data = np.loadtxt(datafilein, unpack=True)
    time = np.loadtxt(timefilein, unpack=True) - matlab_timeoffset
    t_ind = (time<t_max) & (time>=t_min)
    time = time[t_ind]
    
    deptharray = np.loadtxt(depthfilein, unpack=True)
    
    ###data massaging
    dt = 1. #data is hourly
    time_base = 'hours'
    depth = -1. * deptharray

    # signal
    xx = data[t_ind,level]

    variance = np.var(xx)

   
    #normalize
    print 'Variance = %s ' % (variance)
    x = (xx - np.mean(xx)) / np.sqrt(variance)
    variance = np.var(x)

    sl = len(x)
    #time = (np.arange(0,sl ,1) * dt/24. ) + startdate #date_object = datetime.datetime.strptime('04 23 2009', '%m %d %Y')

    return (data, x,dt,np.array(time) * 24., variance, time_base, depth)     
    
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
for level in range(0,11,1):
    print level
    mooring = '12ckp2a'
    (raw_data, x, dt, time, variance, time_base, depth) = \
    ADCP_2D('/Users/bell/Data_Local/from_phyllis/WaveletAnalysisBS2/wcp_data/'+mooring+'/'+mooring+'_ABS.txt', \
    '/Users/bell/Data_Local/from_phyllis/WaveletAnalysisBS2/wcp_data/'+mooring+'/'+mooring+'_dates.txt', \
    '/Users/bell/Data_Local/from_phyllis/WaveletAnalysisBS2/wcp_data/'+mooring+'/'+mooring+'_depth.txt', level=level) 
    fig_name_base = 'images/'+mooring+'_ADCP'

    """-----------------------------wavelet analysis           ---------------------------"""

    wa = WaveletAnalysis(x, time=time, dt=dt, dj=0.125)

    # wavelet power spectrum
    power = wa.wavelet_power
    transform = wa.wavelet_transform

    # scales 
    scales = wa.scales

    # associated time vector
    t = wa.time / 24.

    # reconstruction of the original data
    rx = wa.reconstruction()

    # determine acor factor for red noise
    acorr = acf(x)
    lagf = (acorr[1]+np.sqrt(acorr[2]))/2
    print 'acorr lagf is %s' % lagf

    # determine significance levels
    (signif, fft_theory) = wave_sig.wave_signif(x,dt,scales,lag1=lagf)
    sig95 = np.ones_like(power) * np.array([signif] * len(t)).transpose()
    sig95 = power / sig95         # where ratio > 1, power is significant

    # Global wavelet spectrum & significance levels:
    global_int = variance*(np.sum(power, axis=0) ) / x.shape[0]   # time-average over all times
    gs = ((np.sum(power, axis=1) ) / x.shape[0]) / variance #assume var=1
    gswa = wa.global_wavelet_spectrum
    # Global wavelet significance
    (signif_g, fft_theory_g) = wave_sig.global_wave_signif(x,dt,scales,lag1=lagf,sigtest=1, dof=len(x))




    """----------------------------- plot setup ------------------------------------------"""
    T, S = np.meshgrid(t, scales)


    """----------- plotting WaveTransform Power with confidence interval contour ----------"""

    plt, fig = wavelet_analy_plot.plot_wavetransf(wa, T, S, sig95, time_base, plot_percentile=True)

    plt.savefig((fig_name_base + '_wave_' + str(depth[level]).replace('.0','m') + '.png'), bbox_inches='tight', dpi = (100))
    plt.close()

    """----------------- plotting contours w/global and timeseries ----------"""

    plt, fig = wavelet_analy_plot.plot_wavetransf_time_zoom(x, wa, T, S, sig95, gs, signif_g,
             time_base, scalemin=2, scalemax=256, ylabel='Echo Intens.', plot_percentile=True)
    plt.savefig(fig_name_base + '_wave2_' + str(depth[level]).replace('.0','m') + '.png', bbox_inches='tight', dpi = (100))
    plt.close()

    """----------------- plotting contours w/global and timeseries ----------"""
    """----------------- zoom in on time series ----------"""
    
    t_min = datetime.datetime.strptime('2011-02-15','%Y-%m-%d').toordinal() 
    t_max = datetime.datetime.strptime('2011-11-15','%Y-%m-%d').toordinal() 
    months = 'ALL' #JJA
    plt, fig = wavelet_analy_plot.plot_wavetransf_time_zoom_tbound(x, wa, T, S, sig95, gs, signif_g,
             time_base, t_min, t_max, scalemin=2, scalemax=256, ylabel='Echo Intens.', plot_percentile=True)
    plt.savefig(fig_name_base + '_wave_' + months + str(depth[level]).replace('.0','m') + '.png', bbox_inches='tight', dpi = (100))
    plt.close()
    
    """----------------------- plotting power spectrum FFT --------------------------------"""
    (plt, fig) = wavelet_analy_plot.fft_power_spec(x, time_base, Fs=1)

    plt.savefig(fig_name_base + '_FFTspec_' + str(depth[level]).replace('.0','m') + '.png', bbox_inches='tight', dpi = (100))
    plt.close()

    """
    # Do FFT analysis of array
    sp = np.fft.fft(x)
    # Getting the related frequencies
    freq = np.fft.fftfreq(t.shape[-1], d=.25)
    pyy = sp*np.conj(sp)
    """

    """----------------- plot scale averaged power timeseries -----------------------------"""
    #Need to know what scales to plot
    scales_bin = [8, 30]

    #find indices
    indices = np.where((scales >8) & (scales < 30))

    #scales[indices]
    scale_ave1 = power[indices].mean(axis=0)

    (plt, fig) = wavelet_analy_plot.scale_ave_timeseries(scale_ave1, wa.time / 24. , scales_bin)
    plt.savefig(fig_name_base + '_scaleave_' + str(depth[level]).replace('.0','m') + '.png', bbox_inches='tight', dpi = (100))
    plt.close()

    """-------------------------- plot orig data (2d) -------------------------------------"""

(plt, fig) = wavelet_analy_plot.plot2dvar(raw_data, time / 24., depth)
plt.savefig((fig_name_base + '_rawdata_' + str(depth[level]).replace('.0','m') + '.png'), bbox_inches='tight', dpi = (100))
plt.close()



