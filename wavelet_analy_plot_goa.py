#!/usr/bin/env 
"""
 wavelet_analy_plots_goa.py
 to be used with
 Wavelet_analy_GOA_yearstitch.py
 
 For P. Stabeno
 

   Using Anaconda packaged Python
   modifications for confidence intervals based on wave_matlab at
   http://paos.colorado.edu/research/wavelets/
"""

#Standard packages
import os

#Science packages
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter, MonthLocator
import matplotlib.ticker as plticker
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


    

def plot_wavetransf_time_zoom(x, wa, T, S, sig95, gs, signif_g, time_base, data_ind=True, scalemin=0, scalemax=6, ylabel='Pressure (mb)', plot_percentile=True):
    """plotting contours w/global and timeseries"""

    fig = plt.figure(22)
    ax = plt.subplot2grid((3, 4), (1, 0), colspan=3, rowspan=2)
    # use following with unnormalized data to match web output
    if plot_percentile:
        #use following to contour at "percentiles variances" when using non-normalized data to match web output
        csf =plt.contourf(T, S, wa.wavelet_power[data_ind,:], levels=[ 0, stats.scoreatpercentile(wa.wavelet_power[data_ind,:], 25), stats.scoreatpercentile(wa.wavelet_power[data_ind,:], 50),
                                           stats.scoreatpercentile(wa.wavelet_power[data_ind,:], 75), stats.scoreatpercentile(wa.wavelet_power[data_ind,:], 95), 
                                           stats.scoreatpercentile(wa.wavelet_power[data_ind,:], 100)], colors=bmap)
    else:
        #use following to contour at "normalized variances" BAMS
        csf =plt.contourf(T, S, wa.wavelet_power[data_ind,:], levels=[ 0, 1,2,5,10], colors=bmap)
    levels = [-99, 1] # values greater than 1 are significant
    plt.contour(T, S, sig95,levels, colors='black', linewidth=5)
    """
    cbar = fig.colorbar(csf, ticks=[ 0, stats.scoreatpercentile(wa.wavelet_power[data_ind,:], 25), stats.scoreatpercentile(wa.wavelet_power[data_ind,:], 50),
                                           stats.scoreatpercentile(wa.wavelet_power[data_ind,:], 75), stats.scoreatpercentile(wa.wavelet_power[data_ind,:], 95), 
                                           stats.scoreatpercentile(wa.wavelet_power[data_ind,:], 100)], orientation='horizontal')
    cbar.ax.set_xticklabels(['', '25%', '50%', '75%', '95%', ''])
    """
    ax.set_yscale('log')
    ax.set_ylabel('Scales')
    ax.grid(True)

    # put the ticks at powers of 2 in the scale
    ticks = np.unique(2 ** np.floor(np.log2(wa.scales[data_ind])))[1:]
    ax.yaxis.set_ticks(ticks)
    ax.yaxis.set_ticklabels(ticks.astype(str))
    ax.set_ylim(scalemax, scalemin)

    ax.xaxis.set_major_formatter(DateFormatter('%Y-%m'))
    ax.xaxis.set_major_locator(MonthLocator())
    #fig.autofmt_xdate()

    # shade the region between the edge and coi
    C, S = wa.coi
    ax.fill_between(x=C, y1=S, y2=wa.scales[data_ind].max(), color='gray', alpha=0.5)
    ax.set_xlim(wa.time.min(), wa.time.max())

    #plt.show()
    DefaultSize = fig.get_size_inches()
    fig.set_size_inches( (DefaultSize[0]*2, DefaultSize[1]) )
    
    return(plt, fig)    

