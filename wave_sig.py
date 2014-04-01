#!/usr/bin/env 
"""
 wave_sig.py
 
 References:
 -----------
 
 See "http://paos.colorado.edu/research/wavelets/"
 Written January 1998 by C. Torrence
 
 Translated to Python 2014 

   Using Anaconda packaged Python
   from wave_matlab
"""

import numpy as np

def wave_signif(Y,dt,scale1,lag1=0.0,sigtest=0, dof=0, siglvl=0.95):
    """ Inputs:
        -------
            Y = timeseries (or variance?)
            dt = sampling time
            scale1 = vector of scale indices
            """
            
    n1 = len(Y)
    J1 = len(scale1)
    scale = scale1
    s0 = scale.min()
    dj = np.log(scale[1]/scale[0])/np.log(2.)

    ## for morlet 6 only
    #lag1 = 0.72 - from Torrence 1998 for rednoise NinoSST
    
    param = 6.
    k0 = param
    fourier_factor = (4. * np.pi) / (k0 + np.sqrt(2. + k0**2))
    empir = np.array([2.,-1,-1,-1])
    if (k0 == 6.):
        empir[1]=0.776
        empir[2]=2.32
        empir[3]=0.60
        
    if (np.size(Y) == 1):
        variance = Y
    else:
        variance = np.var(Y)
        
    period = scale * fourier_factor
    dofmin = empir[0]     # Degrees of freedom with no smoothing
    Cdelta = empir[1]     # reconstruction factor
    gamma_fac = empir[2]  # time-decorrelation factor
    dj0 = empir[3]        # scale-decorrelation factor
    
    freq = dt / period;   # normalized frequency
    fft_theor = (1. - lag1**2.) / (1. - 2. * lag1 * np.cos(freq * 2. * np.pi) + lag1**2.)  # [Eqn(16)]
    fft_theor = variance * fft_theor  # include time-series variance
    signif = fft_theor

    if (sigtest==0):
        dof = dofmin
        chisquare = chisquare_inv(siglvl,dof, scipy=True) / dof 
        signif = fft_theor * chisquare   # [Eqn(18)]
    else: 
        print "No options for sigtest != 0"
        raise

   
    return(signif,fft_theor)

def global_wave_signif(Y,dt,scale1,lag1=0.0,sigtest=1, dof=0, siglvl=0.95):
    """Time averaged significance for global wavelet averages"""
            
    n1 = len(Y)
    J1 = len(scale1)
    scale = scale1
    s0 = scale.min()
    dj = np.log(scale[1]/scale[0])/np.log(2.)

    ## for morlet 6 only
    #lag1 = 0.72 - from Torrence 1998 for rednoise NinoSST
    
    param = 6.
    k0 = param
    fourier_factor = (4. * np.pi) / (k0 + np.sqrt(2. + k0**2))
    empir = np.array([2.,-1,-1,-1])
    if (k0 == 6.):
        empir[1]=0.776
        empir[2]=2.32
        empir[3]=0.60
        
    if (np.size(Y) == 1):
        variance = Y
    else:
        variance = np.var(Y)
        
    period = scale * fourier_factor
    dofmin = empir[0]     # Degrees of freedom with no smoothing
    Cdelta = empir[1]     # reconstruction factor
    gamma_fac = empir[2]  # time-decorrelation factor
    dj0 = empir[3]        # scale-decorrelation factor
    
    freq = dt / period;   # normalized frequency
    fft_theor = (1. - lag1**2.) / (1. - 2. * lag1 * np.cos(freq * 2. * np.pi) + lag1**2.)  # [Eqn(16)]
    fft_theor = variance * fft_theor  # include time-series variance
    signif = fft_theor   
    
    if (sigtest==1):
        dof = dofmin * np.sqrt(1. + (dof * dt / gamma_fac / scale1)**2. )
        dof = [dofmin if trunc < dofmin else trunc for trunc in dof ]
        chisquare = chisquare_inv(siglvl,dof, scipy=True) / dof 
        signif = fft_theor * chisquare
    else: 
        print "No options for sigtest != 1"
        raise
 
    return(signif,fft_theor)
    
    
def chisquare_inv(P,V, scipy=True):
    """ Translated from chisquare_inv.m
        Originally coded by C. Torrence January 1998 
        
    By passing flag scipy = True : use scipy stats.chi2.inverse function
        """
    
    if ((1-P) < 1e-4):
        print "Must use a P <0.9999"
        exit()
    elif ((P==0.95) & (V==2)): # from lookup tables
        X = 5.9915
    
        return (X)
    elif  ((P==0.90) & (V==2)): # from lookup tables
        X = 4.605
    
    if scipy:
        import scipy as sp
        X = sp.stats.chi2.isf((1-P),V)

    return (X) 
    
