#!/usr/bin/env 
"""
 wave_xwt.py
 
 Purpose:
 --------
 Calculate Cross Wavelet Transform of 2 time series
 
 References:
 -----------
 These routines are a translation / modification of matlab code from Aslak Grinsted
 which also builds on the software from Torrence and Campo based on wave_matlab at
 http://paos.colorado.edu/research/wavelets/
 
 ::
 
  Please acknowledge the use of this software in any publications:
   
  Crosswavelet and wavelet coherence software were provided by A. Grinsted.
 
  (C) Aslak Grinsted_. 2002-2004
 
  http://www.pol.ac.uk/home/research/waveletcoherence/
  .. _Grinsted: http://www.glaciology.net/Home
  
  -------------------------------------------------------------------------
    Copyright (C) 2002-2004, Aslak Grinsted
    This software may be used, copied, or redistributed as long as it is not
    sold and this copyright notice is reproduced on each copy made.  This
    routine is provided as is without any express or implied warranties
    whatsoever.

   Using Anaconda packaged Python
"""

import numpy as np
import scipy as sp
import sys

"""---------------------------------math routines-------------------------------------"""

def angle(Wxy):
    """input: cross wavelet transfer
       output: angle between real and imaginary comp"""

    P = np.arctan2(  np.imag(Wxy), np.real(Wxy) )    
        
    return ( P )


"""---------------------------------xwt/wtc routines-------------------------------------"""

def xwt(wa1,wa2):
    """ Cross wavelet transfprm
        The regions in time frequency space where the two time series show high common power
        
        wa1 - wavelet analysis of timeseries 1 using WaveletAnalysis routine
        wa2 - wavelet analysis of timeseries 2 using WaveletAnalysis routine
        
        Returns:
        --------
        
        Wxy: Cross Wavelet Transform
        aWxy: Phase Angle
        """
    if (wa1.dt != wa2.dt):
        sys.exit("Timesteps must be equivalent")


    
    Wxy=wa1.wavelet_transform * np.conj(wa2.wavelet_transform);
    aWxy = angle(Wxy)

    return(Wxy, aWxy)
    
def wtc(wa1,wa2, xwt):
    """ Wavelet Coherence
        The regions in time frequency space where the two time series co-vary but does not
        necessarily have high power
        
        wa1 - wavelet analysis of timeseries 1 using WaveletAnalysis routine
        wa2 - wavelet analysis of timeseries 2 using WaveletAnalysis routine
    """
    if (wa1.dt != wa2.dt):
        sys.exit("Timesteps must be equivalent")
        
    print "Smoothing Series 1"
    sw1 = smoothwavelet(wa1)
    print "Smoothing Series 2"
    sw2 = smoothwavelet(wa2)
    print "Smoothing XWT"
    sw_xwt = smoothwavelet_xwt(xwt,wa2.dt, wa2.fourier_periods, wa2.dj, wa2.scales)

    return (sw1, sw2, sw_xwt)
        
def smoothwavelet(wa):
    """
     Input - wa object (WaveletAnalysis object from waveletanalysis package by https://github.com/aaren/wavelets
     Output - smoothed wavelet analysis
    
    
     Note: Very Slow!!!!
    
     Translated to Python from the following:
    
     Smoothing as in the np.appendix of Torrence and Webster "Inter decadal changes in the ENSO-Monsoon System" 1998
     
     used in wavelet coherence calculations
     
    
     Only applicable for the Morlet wavelet. 
    
     (C) Aslak Grinsted 2002-2005
    

    -------------------------------------------------------------------------
    Copyright (C) 2002-2005, Aslak Grinsted
    This software may be used, copied, or redistributed as long as it is not
    sold and this copyright notice is reproduced on each copy made.  This
    routine is provided as is without any express or implied warranties
    whatsoever.

    """

    
    n = wa.wavelet_power.shape
    twave = np.zeros_like(wa.wavelet_power)

    #zero-pad to power of 2... Speeds up fft calcs if n is large
    npad=2.**np.ceil(np.log2(n[1]))

    k = np.arange( 1,np.fix(npad/2)+1 )
    k = k * ((2. * np.pi ) / npad )
    k1 = np.array([0.])
    k1 = np.hstack((k1, k))
    k1 = np.hstack((k1, -1 * k[::-1][1:]))

    k2=k1 ** 2.


    snorm = wa.scales / wa.dt
    sinv = 1 / wa.scales
    wa_invpower = (np.ones_like(wa.wavelet_power) * np.array([sinv] * n[1]).transpose()) * wa.wavelet_power

    for ii in range(0,n[0],1):
        F = np.exp(-.5 * (snorm[ii] ** 2.) * k2)
        smooth=np.fft.ifft(F * np.fft.fft(wa_invpower[ii,:],n=np.int(npad)))
        twave[ii] = np.real(smooth[0:n[1]])

    
    dj0      = 0.6
    dj0steps = dj0 / (wa.dj * 2)
    kernal = np.array([dj0steps % 1.])
    kernal = np.vstack((kernal,np.ones((2. * np.round(dj0steps)-1,1))))
    kernal = np.vstack((kernal,dj0steps % 1))
    kernal = kernal / (2. * np.round(dj0steps) - 1. + 2. * dj0steps % 1.)

    swave=sp.signal.convolve2d(twave,kernal,mode='same') 
    
    return (swave)
    
def smoothwavelet_xwt(xwt, dt, period, dj, scales):
    """
     Input - wa object (WaveletAnalysis object from waveletanalysis package by https://github.com/aaren/wavelets
     Output - smoothed wavelet analysis
    
    
     Note: Very Slow!!!!
    
     Translated to Python from the following:
    
     Smoothing as in the np.appendix of Torrence and Webster "Inter decadal changes in the ENSO-Monsoon System" 1998
     
     used in wavelet coherence calculations
     
    
     Only applicable for the Morlet wavelet. 
    
     (C) Aslak Grinsted 2002-2005
    

    -------------------------------------------------------------------------
    Copyright (C) 2002-2005, Aslak Grinsted
    This software may be used, copied, or redistributed as long as it is not
    sold and this copyright notice is reproduced on each copy made.  This
    routine is provided as is without any express or implied warranties
    whatsoever.

    """
    
    n = xwt.shape
    twave = np.zeros_like(xwt)

    #zero-pad to power of 2... Speeds up fft calcs if n is large
    npad=2.**np.ceil(np.log2(n[1]))

    k = np.arange( 1,np.fix(npad/2)+1 )
    k = k * ((2. * np.pi ) / npad )
    k1 = np.array([0.])
    k1 = np.hstack((k1, k))
    k1 = np.hstack((k1, -1 * k[::-1][1:]))

    k2=k1 ** 2.
    
    snorm = scales / dt
    sinv = 1 / scales
    wa_invpower = np.ones_like(xwt) * np.array([sinv] * n[1]).transpose() * xwt
    
    for ii in range(0,n[0],1):
        F = np.exp(-.5 * (snorm[ii] ** 2.) *k2 )
        smooth=np.fft.ifft(F * np.fft.fft(wa_invpower[ii,:],n=np.int(npad)))
        twave[ii] = np.real(smooth[0:n[1]])

    #twave = np.real(twave)
        
    dj0      = 0.6
    dj0steps = dj0 / (dj * 2)
    kernal = np.array([dj0steps % 1.])
    kernal = np.vstack((kernal,np.ones((2. * np.round(dj0steps)-1,1))))
    kernal = np.vstack((kernal,dj0steps % 1))
    kernal = kernal / (2. * np.round(dj0steps) - 1. + 2. * dj0steps % 1.)

    swave=sp.signal.convolve2d(twave,kernal,mode='same') 
    
    

    return (swave)
