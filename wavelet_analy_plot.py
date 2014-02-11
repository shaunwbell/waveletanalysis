#!/usr/bin/env 
"""
 Wavelet_analy_plots.py
 to be used with
 Wavelet_analy.py
 
 For P. Stabeno
 

   Using Anaconda packaged Python
   modifications for confidence intervals based on wave_matlab at
   http://paos.colorado.edu/research/wavelets/
"""

#Standard packages
import os

#Science packages
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
import numpy as np
from scipy import stats

#Other Packages
import brewer2mpl as b2m #colorbrewer maps

"""----------------------------- plot setup ------------------------------------------"""

#set up color map with brewer2mpl, not mandatory as matplotlib as the colorbrewer maps but
# provides an easy way to limit the number of colors in a map
bmap = b2m.get_map('Blues', 'Sequential', 5, reverse=False)
bmap = bmap.mpl_colors

"""----------------------------- plots  ----------------------------------------------"""

def plot_wavetransf(wa, T, S, sig95, time_base, plot_percentile=False):
    """plotting WaveTransform Power with confidence interval contour"""

    fig = plt.figure(1)
    ax = plt.subplot(1,1,1)
    if plot_percentile:
        #use following to contour at "percentiles variances" when using non-normalized data to match web output
        csf =plt.contourf(T, S, wa.wavelet_power, levels=[ 0, stats.scoreatpercentile(wa.wavelet_power, 25), stats.scoreatpercentile(wa.wavelet_power, 50),
                                           stats.scoreatpercentile(wa.wavelet_power, 75), stats.scoreatpercentile(wa.wavelet_power, 95), 
                                           stats.scoreatpercentile(wa.wavelet_power, 100)], colors=bmap)
    else:
        #use following to contour at "normalized variances" BAMS
        csf =plt.contourf(T, S, wa.wavelet_power, levels=[ 0, 1,2,5,10], colors=bmap)
    cbar = plt.colorbar(pad=.1, shrink=.5, format='%.4f', extend='both') #move and shrink colorbar
    levels = [-99, 1] # values greater than 1 are significant
    plt.contour(T, S, sig95,levels, colors='black', linewidth=5)
    ax.set_yscale('log')
    ax.grid(True)

    # put the ticks at powers of 2 in the scale
    ticks = np.unique(2 ** np.floor(np.log2(wa.scales)))[1:]
    ax.yaxis.set_ticks(ticks)
    ax.yaxis.set_ticklabels(ticks.astype(str))
    ax.set_ylim(256, 0.5)
    ax.set_ylabel('scales')

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

    ax.xaxis.set_major_formatter(DateFormatter('%Y-%m-%d'))
    #fig.autofmt_xdate()

    # shade the region between the edge and coi
    C, S = wa.coi
    ax.fill_between(x=C, y1=S, y2=wa.scales.max(), color='gray', alpha=0.5)
    ax.set_xlim(wa.time.min(), wa.time.max())
        
    #plt.show()
    DefaultSize = fig.get_size_inches()
    fig.set_size_inches( (DefaultSize[0]*2, DefaultSize[1]) )
    
    return (plt, fig)
    
def plot_wavetransf_time(x, wa, T, S, sig95, gs, signif_g, time_base, ylabel='Pressure (mb)', plot_percentile=False):
    """plotting contours w/global and timeseries"""

    fig = plt.figure(2)
    ax = plt.subplot2grid((3, 4), (1, 0), colspan=3, rowspan=2)
    # use following with unnormalized data to match web output
    if plot_percentile:
        #use following to contour at "percentiles variances" when using non-normalized data to match web output
        csf =plt.contourf(T, S, wa.wavelet_power, levels=[ 0, stats.scoreatpercentile(wa.wavelet_power, 25), stats.scoreatpercentile(wa.wavelet_power, 50),
                                           stats.scoreatpercentile(wa.wavelet_power, 75), stats.scoreatpercentile(wa.wavelet_power, 95), 
                                           stats.scoreatpercentile(wa.wavelet_power, 100)], colors=bmap)
    else:
        #use following to contour at "normalized variances" BAMS
        csf =plt.contourf(T, S, wa.wavelet_power, levels=[ 0, 1,2,5,10], colors=bmap)
    levels = [-99, 1] # values greater than 1 are significant
    plt.contour(T, S, sig95,levels, colors='black', linewidth=5)
    ax.set_yscale('log')
    ax.set_ylabel('Scales')
    ax.grid(True)

    # put the ticks at powers of 2 in the scale
    ticks = np.unique(2 ** np.floor(np.log2(wa.scales)))[1:]
    ax.yaxis.set_ticks(ticks)
    ax.yaxis.set_ticklabels(ticks.astype(str))
    ax.set_ylim(256, 0.5)

    ax.xaxis.set_major_formatter(DateFormatter('%Y-%m-%d'))
    #fig.autofmt_xdate()

    # shade the region between the edge and coi
    C, S = wa.coi
    ax.fill_between(x=C, y1=S, y2=wa.scales.max(), color='gray', alpha=0.5)
    ax.set_xlim(wa.time.min(), wa.time.max())

    ax = plt.subplot2grid((3, 4), (0, 0), colspan=3, rowspan=1)
    p1 = ax.plot(wa.time,x,'r', wa.time, wa.reconstruction(), 'b')
    ax.set_xlim([wa.time.min(), wa.time.max()])
    ax.set_xticklabels([])
    ax.grid(True)
    ax.set_ylabel(ylabel)
    ax2= ax.twinx()
    p2 = ax2.plot(wa.time,x-wa.reconstruction(), 'k')
    ax2.set_xlim([wa.time.min(), wa.time.max()])
    ax2.set_yticklabels([])

    ax = plt.subplot2grid((3, 4), (1, 3), colspan=1, rowspan=2)
    p1 = ax.plot(gs,wa.scales, signif_g, wa.scales, 'k--')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.grid(True)

    # put the ticks at powers of 2 in the scale
    ticks = np.unique(2 ** np.floor(np.log2(wa.scales)))[1:]
    ax.yaxis.set_ticks(ticks)
    ax.set_yticklabels([])
    #ax.yaxis.set_ticklabels(ticks.astype(str))
    ax.set_ylim(256, 0.5)

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

    #plt.show()
    DefaultSize = fig.get_size_inches()
    fig.set_size_inches( (DefaultSize[0]*2, DefaultSize[1]) )
    
    return(plt, fig)
    
def scaleogram(wa):
    """scaleogram"""

    fig = plt.figure(3)
    ax = plt.subplot(1,1,1)
    csf = plt.imshow(wa.wavelet_power,  
            cmap='jet', aspect='auto', origin='lower', extent=[wa.time.min(),wa.time.max(),wa.scales.min(),wa.scales.max()])
    ax.set_ylim(ax.get_ylim()[::-1]) #this reverses the yaxis (i.e. deep at the bottom)
    cbar = plt.colorbar()
    
    ax.xaxis.set_major_formatter(DateFormatter('%Y-%m-%d'))
    fig.autofmt_xdate()

    return(plt, fig)
    
def fft_power_spec(x, time_base):
    """Plot original data (normalized) FFT power spectral density"""

    fig = plt.figure(4)
    ax1 = plt.subplot(1,1,1)
    ax1.psd(x, NFFT=256, pad_to=None, noverlap=0, Fs=4)
    ax1.set_ylabel('Power Spectral Density dB/(cycles/' + time_base + ')')
    ax1.set_xlabel('Frequency (cycles/' + time_base + ')')
    ax1.set_xscale('log')
    #ax4.set_ylabel('')
    #plt.title('overlap')
    
    return(plt, fig)
    
def scale_ave_timeseries(scale_ave, time, scales_bin):
    """plotting WaveTransform scale averaged power for selected bins"""

    fig = plt.figure(5)
    ax1 = plt.subplot(1,1,1)
    p1 = ax1.plot(time,scale_ave,'k')
    ax1.set_xlim([time.min(), time.max()])
    ax1.grid(True)
    ax1.set_ylabel('Average power between scales ' + str(scales_bin[0]) + ' and ' + str(scales_bin[1]) + '')
    ax1.set_xlabel('Time (UTC)')
    ax1.xaxis.set_major_formatter(DateFormatter('%Y-%m-%d'))

    #plt.show()
    DefaultSize = fig.get_size_inches()
    fig.set_size_inches( (DefaultSize[0]*2, DefaultSize[1]) )

    return(plt, fig)

def scale_ave_timeseries2D(scale_ave1, scale_ave2, time1, time2, scales_bin):
    """plotting WaveTransform scale averaged power for selected bins"""

    fig = plt.figure(5)
    ax1 = plt.subplot(1,1,1)
    p1 = ax1.plot(time1,scale_ave1,'k', time2,scale_ave2,'r--')
    ax1.set_xlim([time1.min(), time1.max()])
    ax1.grid(True)
    ax1.set_ylabel('Average power between scales ' + str(scales_bin[0]) + ' and ' + str(scales_bin[1]) + '')
    ax1.set_xlabel('Time (UTC)')
    ax1.xaxis.set_major_formatter(DateFormatter('%Y-%m-%d'))

    #plt.show()
    DefaultSize = fig.get_size_inches()
    fig.set_size_inches( (DefaultSize[0]*2, DefaultSize[1]) )

    return(plt, fig)