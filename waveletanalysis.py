#!/usr/bin/env 
"""
 textreadutilities.py
 

   Using Anaconda packaged Python
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from subroutines.wavelets_bams.wavelets import WaveletAnalysis

def example():

    # given a signal x(t)
    x = np.random.randn(1000)
    # and a sample spacing
    dt = 0.1

    wa = WaveletAnalysis(x, dt=dt)

    # wavelet power spectrum
    power = wa.wavelet_power

    # scales 
    scales = wa.scales

    # associated time vector
    t = wa.time

    # reconstruction of the original data
    rx = wa.reconstruction()

def nino3test():

    from utilities import textreadutilities as tutil

    ###text readin
    dir_path = os.path.dirname(os.path.abspath(__file__))
    infile = '/data/nino34seasonalanom.txt'    
    SSTfile = dir_path +infile 

    (season, year, temp, anom) = tutil.read_txt(SSTfile)

    ###data massaging
    dt = 1.0 / 12.0
    time_base = 'years'
    sl = len(season)
    timet = range(0,12,1) * int(sl / 12)
    time = np.array(timet) /12. + year

    # signal
    x = anom
    variance = np.var(x)
    
    return (x,dt,np.array(time), variance, time_base)
    
    
def shortninoexample():
    """ data file from http://paos.colorado.edu/research/wavelets/software.html
    Used to validate program retrievals.
    Should be consistant with the website as well as Torrence and Campo 1997 BAMS article"""
    from utilities import textreadutilities as tutil

    ###text readin
    dir_path = os.path.dirname(os.path.abspath(__file__))
    infile = '/data/sst_nino3.dat'    
    SSTfile = dir_path +infile 

    anom = np.loadtxt(SSTfile, unpack=True)

    ###data massaging
    dt = 0.25
    time_base = 'years'
    sl = len(anom)
    time = (np.arange(0,sl ,1) * dt) +1871.

    # signal
    x = anom
    variance = np.var(x)
    
    #normalize
    x = (x - np.mean(x)) / np.sqrt(variance)
    variance = np.var(x)

    return (x,dt,np.array(time), variance, time_base)

def testnino3long():
    
    from utilities import ncutilities as ncutil
    from utilities import utilities as util
    from utilities import constants 

    dir_path = os.path.dirname(os.path.abspath(__file__))    
    ncfile = dir_path + '/data/nino3data.cdf'

    ###nc readin/out
    nchandle = ncutil.ncopen(ncfile)
    params = ['T', 'NINO3']

    ncdata = np.zeros( (nchandle.variables['T'][:].shape[0],len(params)) ) 
    for j, v in enumerate(params): 
        if v in nchandle.variables.keys(): #check for nc variable
            ncdata[:,j] = nchandle.variables[v][:]
        else: #if parameter doesn't exist fill the array with zeros
            ncdata[:,j] = data[:,j-1] * 0.0

    ncutil.ncclose(nchandle)


    ###data massaging
    dt = 1.0 / 12.0
    x = ncdata[:,1]
    time = ncdata[:,0]
    t0=715510 ###Jan 1, 1960
    time = t0 + time * 30. #a bit of a cheat... assuming 30days per month
    variance = np.var(x)
    time_base = 'years'

    variance = np.var(x)
    #normalize
    x = (x - np.mean(x)) / np.sqrt(variance)
    variance = np.var(x)

    return (x,dt,np.array(time), variance, time_base)

def uv_current_example():

    from utilities import ncutilities as ncutil
    from utilities import utilities as util
    from utilities import constants 

    dir_path = os.path.dirname(os.path.abspath(__file__))    
    ncfile = dir_path + '/data/test.nc'

    ###nc readin/out
    nchandle = ncutil.ncopen(ncfile)
    params = constants.nc_vars_moor()
    ncdata = ncutil.ncreadfile(nchandle, params)
    ncutil.ncclose(nchandle)

    ###data massaging
    time_all = ncdata[:,0] + ncdata[:,1]
    xx = ncutil.nc_missing(ncdata[:,12], flag=1e35, setfill='Zero')
    pytime = util.EPICdate2udunits(ncdata[:,0], ncdata[:,1])
    
    x = util.moving_average(xx, int(pytime.get('interval_min')) / 4, type='simple')
    x = util.subsample(x, int(pytime.get('interval_min')) / 4) #one/fourth day subsample of averaged data
    dt = 1.  / 4.
    time_base = 'days'
    time = pytime['timeint']
    time = util.subsample(time, int(pytime.get('interval_min')) / 4)

    variance = np.var(x)
    #normalize
    x = (x - np.mean(x)) / np.sqrt(variance)
    variance = np.var(x)

    return (x,dt,np.array(time), variance, time_base)


if __name__ == "__main__":
    (x, dt, time, variance, time_base) = shortninoexample() 

import wave_sig
    
wa = WaveletAnalysis(x, time=time, dt=dt, dj=0.25)

# wavelet power spectrum
power = wa.wavelet_power

# scales 
scales = wa.scales

# associated time vector
t = wa.time

# reconstruction of the original data
rx = wa.reconstruction()

# determine significance levels

(signif, fft_theory) = wave_sig.wave_signif(x,dt,scales, lag1=0.72)
sig95 = np.ones_like(power) * np.array([signif] * len(t)).transpose()
sig95 = power / sig95         # where ratio > 1, power is significant


# plotting time series
fig = plt.figure(1)
ax1 = plt.subplot(1,1,1)
ax1.grid(True)
p1 = ax1.plot(time,x,'r', time, rx, 'b', time,x-rx, 'k')
ax1.autoscale_view(tight=True)
#ax1.xaxis.set_major_formatter(DateFormatter('%Y-%m-%d'))
#fig.autofmt_xdate()

f = plt.gcf()
DefaultSize = f.get_size_inches()
f.set_size_inches( (DefaultSize[0], DefaultSize[1]) )
plt.show()


# plotting contours
fig, ax = plt.subplots()
T, S = np.meshgrid(t, scales)
levels= [1, 2, 5, 10, 20]
csf =plt.contourf(T, S, power, levels, cmap='Blues')
levels = [-99, 1] # values greater than 1 are significant
plt.contour(T, S, sig95,levels)
ax.set_yscale('log')
ax.grid(True)

# put the ticks at powers of 2 in the scale
ticks = np.unique(2 ** np.floor(np.log2(scales)))[1:]
ax.yaxis.set_ticks(ticks)
ax.yaxis.set_ticklabels(ticks.astype(str))
ax.set_ylim(64, 0.5)

# second y scale with equivalent fourier periods to scales
# except with the ticks at the powers of 2
ax_fourier = ax.twinx()
ax_fourier.set_yscale('log')
# match the fourier ticks to the scale ticks
ax_fourier.set_yticks(ticks)
ax_fourier.set_yticklabels(ticks.astype(str))
ax_fourier.set_ylabel('fourier period (%s)' % time_base )
fourier_lim = [wa.fourier_period(i) for i in ax.get_ylim()]
ax_fourier.set_ylim(fourier_lim)

#ax.xaxis.set_major_formatter(DateFormatter('%Y-%m-%d'))
#fig.autofmt_xdate()

# shade the region between the edge and coi
C, S = wa.coi
ax.fill_between(x=C, y1=S, y2=scales.max(), color='gray', alpha=0.3)
ax.set_xlim(t.min(), t.max())
        
plt.show()

