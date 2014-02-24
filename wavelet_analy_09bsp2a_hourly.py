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


#User Packages
import wave_sig
from general_utilities.wavelets_bams.wavelets import WaveletAnalysis
from utilities import ncutilities as ncutil
from utilities import utilities as util
from utilities import constants 
import wavelet_analy_plot_hourly as wavelet_analy_plot

__author__   = 'Shaun Bell'
__email__    = 'shaun.bell@noaa.gov'
__created__  = datetime.datetime(2014, 01, 13)
__modified__ = datetime.datetime(2014, 01, 13)
__version__  = "0.1.0"
__status__   = "Development"


"""----------------------Many different methods for reading raw data-------------------"""


def reanaly_sfc_press(file_in):
    """ NetCDF files from Reanalysis webportal"""


    ncfile = file_in
    
    ###nc readin/out
    nchandle = Dataset(ncfile,'r')
    params = ['time','lat', 'lon', 'pres']
    time = nchandle.variables[params[0]][:]
    lat = nchandle.variables[params[1]][:]
    lon = nchandle.variables[params[2]][:]
    data = nchandle.variables[params[3]][:,:]
    ncutil.ncclose(nchandle)
    
    site_loc = [52.5, 190.0] #52.8, -169.3 for M6 Samalga mooring
    ind_lat = np.where(lat == site_loc[0])[0][0]
    ind_lon = np.where(lon == site_loc[1])[0][0]
    
    xx = data[:,ind_lat,ind_lon] / 100 #pa to hPa/mbar

    # convert time to serial date
    base_date = datetime.datetime.fromordinal(1) #time is in hours since 1-1-1
    time_delta = datetime.timedelta(hours = 1)
                                                                      # there is a -2 correction factor in the following
                                                                      # conversion to calculate the date (unknown why its needed)
    time = [base_date + (int(t - 48) * time_delta) for t in time]       #convert to integer for datetime calculation
    time = [t.toordinal() for t in time ]

    variance = np.var(xx)
    #normalize
    print 'Variance = %s ' % (variance)
    x = (xx - np.mean(xx)) / np.sqrt(variance)
    variance = np.var(x)

    return (xx, x,dt,np.array(time), variance, time_base) 
    
def reanaly_sfc_press_multifile(files_in):
    """ NetCDF files from Reanalysis webportal"""

    time = []
    xx = []
    for ncfile in files_in:
    
        ###nc readin/out
        nchandle = Dataset(ncfile,'r')
        params = ['time','lat', 'lon', 'pres']
        time = np.hstack((time , nchandle.variables[params[0]][:]))
        lat = nchandle.variables[params[1]][:]
        lon = nchandle.variables[params[2]][:]
        data = nchandle.variables[params[3]][:,:]
        ncutil.ncclose(nchandle)
    
        site_loc = [52.5, 190.0] #52.8, -169.3 for M6 Samalga mooring
        ind_lat = np.where(lat == site_loc[0])[0][0]
        ind_lon = np.where(lon == site_loc[1])[0][0]
    
        xx = np.hstack((xx, data[:,ind_lat,ind_lon] / 100)) #pa to hPa/mbar
    dt = 1. #data is daily
    time_base = 'days'
    
    # convert time to serial date
    base_date = datetime.datetime.fromordinal(1) #time is in hours since 1-1-1
    time_delta = datetime.timedelta(hours = 1)
                                                                      # there is a -2 correction factor in the following
                                                                      # conversion to calculate the date (unknown why its needed)
    time = [base_date + (int(t - 48) * time_delta) for t in time]       #convert to integer for datetime calculation

    time = [t.toordinal() for t in time ]

    variance = np.var(xx)
    #normalize
    print 'Variance = %s ' % (variance)
    x = (xx - np.mean(xx)) / np.sqrt(variance)
    variance = np.var(x)

    return (xx, x,dt,np.array(time), variance, time_base) 

def mooring_2dvar(ncfile, level):
    """Standard EcoFOCI Mooring .nc files with two dimensional parameters as a function of time
    (such as ein - echo intensity
    Timestep of data is assumed to be in fractions of a day"""

    ###nc readin/out
    nchandle = Dataset(ncfile,'r')
    params = ['time', 'time2', 'depth','latitude', 'longitude', 'AGC1_1221']
    time = nchandle.variables[params[0]][:]
    time2 = nchandle.variables[params[1]][:]
    lat = nchandle.variables[params[3]][:]
    lon = nchandle.variables[params[4]][:]
    depth = nchandle.variables[params[2]][:]
    ncdata = nchandle.variables[params[5]][:,:,0,0]
    nchandle.close()
    
    dt = 1. #data is hourly
    time_base = 'hours'
    
    pytime = util.EPICdate2udunits(time, time2)
    
    xx = ncdata[:,level]
    
    dt = 24. * (1. / pytime['interval_min']) #data is 4 times daily
    print dt
    time = pytime['timeint']

    variance = np.var(xx)
    #normalize
    print 'Variance = %s ' % (variance)
    x = (xx - np.mean(xx)) / np.sqrt(variance)
    variance = np.var(x)

    return (ncdata, x,dt,np.array(time) * 24., variance, time_base, depth) 

def mooring_1dvar(file_in):
    """Standard EcoFOCI Mooring .nc files with one dimensional parameters as a function of time
    such as temperature at a particular level
    Timestep of data is assumed to be in fractions of a day"""
    #dir_path = os.path.dirname(os.path.abspath(__file__))    
    ncfile = file_in
    
    ###nc readin/out
    nchandle = ncutil.ncopen(ncfile)
    params = constants.nc_vars_moor()
    ncdata = ncutil.ncreadfile(nchandle, params)
    ncutil.ncclose(nchandle)

    ###data massaging
    #time_all = ncdata[:,0] + ncdata[:,1]
    xx = ncutil.nc_missing(ncdata[:,9], flag=1e35, setfill='Zero')
    pytime = util.EPICdate2udunits(ncdata[:,0], ncdata[:,1])
    

    dt = 1. / pytime['interval_min'] #data is 4 times daily
    time_base = 'days'
    time = pytime['timeint']
    #time = util.subsample(time, int(pytime.get('interval_min')) / 4)

    variance = np.var(xx)
    #normalize
    print 'Variance = %s ' % (variance)
    x = (xx - np.mean(xx)) / np.sqrt(variance)
    variance = np.var(x)

    return (xx, x,dt,np.array(time), variance, time_base)  


def ADCP_ABS(datafilein, depthfilein,level,startdate=733530.):
    """ 
    datafile used yyyy_ADCP_ABS.txt - 12 rows, n columns for time
    ADCP backscatter dataset from Adam
    """

    data = np.loadtxt(datafilein, unpack=True)
    depth = np.loadtxt(depthfilein, unpack=True)
    
    ###data massaging
    dt = 1. #data is hourly
    time_base = 'hours'
    depth = depth * -1.

    # signal

    xx = data[:,level]
    variance = np.var(xx)
    
    #normalize
    print 'Depth is %s: ' % (depth[level])
    print 'Variance = %s ' % (variance)
    x = (xx - np.mean(xx)) / np.sqrt(variance)
    variance = np.var(x)

    sl = len(x)
    time = (np.arange(0,sl ,1) * dt /24.) + startdate #date_object = datetime.datetime.strptime('04 23 2009', '%m %d %Y')

    return (data, x,dt,np.array(time) * 24., variance, time_base, depth) 

def ADCP_TAPS_BV(datafilein, timefilein, matlab_timeoffset=366.):
    """ 
    datafile used yyyy_BV.txt - n rows, 1 columns for time
    timefile used yyyy_date.txt - n rows, 1 columns for time
    
    TAPS dataset from Adam
    """

    data = np.loadtxt(datafilein, unpack=True)
    time = np.loadtxt(timefilein, unpack=True) - matlab_timeoffset
    
    ###data massaging
    dt = 1. #data is hourly
    time_base = 'hours'
    depth = 18

    # signal

    xx = data
    variance = np.var(xx)
    
    #normalize
    print 'Depth is %s: ' % (depth)
    print 'Variance = %s ' % (variance)
    x = (xx - np.mean(xx)) / np.sqrt(variance)
    variance = np.var(x)

    sl = len(x)
    #time = (np.arange(0,sl ,1) * dt /24.) + startdate #date_object = datetime.datetime.strptime('04 23 2009', '%m %d %Y')

    return (data, x,dt,np.array(time) * 24., variance, time_base, depth) 

def ADCP_ABS_filtered(datafilein, startdate=733530.):
    """ 
    datafile used yyyy_ADCP_ABS_Lanzcos.txt - n rows, 1 columns for time
    ADCP backscatter dataset from Adam
    """

    data = np.loadtxt(datafilein, unpack=True)
    
    ###data massaging
    dt = 1. #data is hourly
    time_base = 'hours'
    depth = 17

    # signal

    xx = data
    variance = np.var(xx)
    
    #normalize
    print 'Depth is %s: ' % (depth)
    print 'Variance = %s ' % (variance)
    x = (xx - np.mean(xx)) / np.sqrt(variance)
    variance = np.var(x)

    sl = len(x)
    time = (np.arange(0,sl ,1) * dt /24.) + startdate #date_object = datetime.datetime.strptime('04 23 2009', '%m %d %Y')

    return (data, x,dt,np.array(time) * 24., variance, time_base, depth) 


def example():
    """ 
    Data file from http://paos.colorado.edu/research/wavelets/software.html
    Used to validate program retrievals.
    Should be consistant with the website as well as Torrence and Campo 1997 BAMS article
    """
    ###text readin
    dir_path = os.path.dirname(os.path.abspath(__file__))
    infile = '/data/sst_nino3.dat'    
    SSTfile = dir_path +infile 

    anom = np.loadtxt(SSTfile, unpack=True)

    ###data massaging
    dt = 1. / 4. #data is 4 times yearly
    time_base = 'years'
    sl = len(anom)
    time = (np.arange(0,sl ,1) * dt) +1871.

    # signal
    x = anom
    variance = np.var(x)
    
    #normalize
    print 'Variance = %s ' % (variance)
    x = (x - np.mean(x)) / np.sqrt(variance)
    variance = np.var(x)

    return (anom, x,dt,np.array(time), variance, time_base)  

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

for level in range(0,1,1): #if data is 2-dimensional (multiple depths)
    """Example:
    (raw_data, x, dt, time, variance, time_base) = example()
    fig_name_base = 'images/Nino3SST'
    """
    """
    #using ein data
    print level
    #(raw_data, x, dt, time, variance, time_base, depth) = mooring_2dvar('/Users/bell/Data_Local/from_phyllis/wcp_data/09bsp2a/09bsp2a_wcp_ein.nc', level=level) 
    #fig_name_base = 'images/bsp2a_ein'
    #(raw_data, x, dt, time, variance, time_base, depth) = mooring_2dvar('/Users/bell/Data_Local/from_phyllis/wcp_data/10bsp2a/10bsp2a_wcp_ein.nc', level=level) 
    #fig_name_base = 'images/bsp2a_ein'
    """
    """
    #using backscatter data
    (raw_data, x, dt, time, variance, time_base, depth) = \
    ADCP_ABS('/Users/bell/Data_Local/from_phyllis/Phyllis Data/2009_ADCP_ABS.txt', \
    '/Users/bell/Data_Local/from_phyllis/Phyllis Data/2009_ADCP_depth.txt', level=level, startdate=733530.) 
    fig_name_base = 'images/bsp2a_ADCP'

    (raw_data, x, dt, time, variance, time_base, depth) = \
    ADCP_ABS('/Users/bell/Data_Local/from_phyllis/Phyllis Data/2007_ADCP_ABS.txt', \
    '/Users/bell/Data_Local/from_phyllis/Phyllis Data/2007_ADCP_depth.txt', level=level, startdate=732792.) 
    fig_name_base = 'images/bsp2a_ADCP'

    (raw_data, x, dt, time, variance, time_base, depth) = \
    ADCP_ABS('/Users/bell/Data_Local/from_phyllis/Phyllis Data/2010_ADCP_ABS.txt', \
    '/Users/bell/Data_Local/from_phyllis/Phyllis Data/2010_ADCP_depth.txt', level=level, startdate=733890.) 
    fig_name_base = 'images/bsp2a_ADCP'
    """
    """
    #using Lanzcos filtered data 
    (raw_data, x, dt, time, variance, time_base, depth) = \
    ADCP_ABS_filtered('/Users/bell/Data_Local/from_phyllis/Phyllis Data/2009_ADCP_ABS_Lanzcos_17m.txt', startdate=732800.) 
    fig_name_base = 'images/bsp2a_ADCP_Lanzcos'    
    depth = [17]
    """
    """
    (raw_data, x, dt, time, variance, time_base, depth) = \
    ADCP_TAPS_BV('/Users/bell/Data_Local/from_phyllis/Phyllis Data/2007_BV.txt', \
    '/Users/bell/Data_Local/from_phyllis/Phyllis Data/2007_date.txt') 
    fig_name_base = 'images/bsp2a_TAPS'    
    depth = [17]
    """
    (raw_data, x, dt, time, variance, time_base, depth) = \
    ADCP_TAPS_BV('/Users/bell/Data_Local/from_phyllis/Phyllis Data/rebspwavelet2/2007E2_BV.txt', \
    '/Users/bell/Data_Local/from_phyllis/Phyllis Data/rebspwavelet2/2007E2_date.txt') 
    fig_name_base = 'images/2007_bsp2a_rebsp2'    
    depth = [17]    
    
    
    
    """ Begin Routines Below"""
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

    plt.savefig((fig_name_base + '_wave' + str(depth[level]).replace('.0','m') + '.png'), bbox_inches='tight', dpi = (100))
    plt.close()

    """----------------- plotting contours w/global and timeseries ----------"""

    plt, fig = wavelet_analy_plot.plot_wavetransf_time_zoom(x, wa, T, S, sig95, gs, signif_g,
             time_base, scalemin=2, scalemax=256, ylabel='Echo Intens.', plot_percentile=True)
    plt.savefig(fig_name_base + '_wave2' + str(depth[level]).replace('.0','m') + '.png', bbox_inches='tight', dpi = (100))
    plt.close()

    """----------------- plotting contours w/global and timeseries ----------"""
    """----------------- zoom in to specified scales ----------"""

    plt, fig = wavelet_analy_plot.plot_wavetransf_time_zoom(x, wa, T, S, sig95, gs, signif_g, 
        time_base, scalemin=4, scalemax=36, ylabel='Echo Intens.', plot_percentile=True)
    plt.savefig(fig_name_base + '_wave3' + str(depth[level]).replace('.0','m') + '.png', bbox_inches='tight', dpi = (100))
    plt.close()

    """----------------------- plotting power spectrum FFT --------------------------------"""
    (plt, fig) = wavelet_analy_plot.fft_power_spec(x, time_base, Fs=1)

    plt.savefig(fig_name_base + '_FFTspec' + str(depth[level]).replace('.0','m') + '.png', bbox_inches='tight', dpi = (100))
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
    plt.savefig(fig_name_base + '_scaleave' + str(depth[level]).replace('.0','m') + '.png', bbox_inches='tight', dpi = (100))
    plt.close()

    """-------------------------- plot orig data (2d) -------------------------------------"""

(plt, fig) = wavelet_analy_plot.plot2dvar(raw_data, time / 24., depth)
plt.savefig((fig_name_base + '_rawdata' + str(depth[level]).replace('.0','m') + '.png'), bbox_inches='tight', dpi = (100))
plt.close()

