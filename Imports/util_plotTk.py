#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     util_plotTk.py                                                    #
#                                                                             #
# PURPOSE:  Plotting functions for the POSSUM pipeline Tk interface.          #
#                                                                             #
# MODIFIED: 19-November-2015 by C. Purcell                                    #
#                                                                             #
# CONTENTS:                                                                   #
#                                                                             #
# xfloat                                                                      #
# xint                                                                        #
# filter_range_indx                                                           #
# tweakAxFormat                                                               #
# format_ticks                                                                #
# plot_I_vs_nu_ax                                                             #
# plot_PQU_vs_nu_ax                                                           #
# plot_rmsIQU_vs_nu_ax                                                        #
# plot_pqu_vs_lamsq_ax                                                        #
# plot_psi_vs_lamsq_ax                                                        #
# plot_q_vs_u_ax                                                              #
# plot_RMSF_ax                                                                #
# gauss                                                                       #
# plot_dirtyFDF_ax                                                            #
# plot_cleanFDF_ax                                                            #
# plot_hist4_ax                                                               #
#                                                                             #
# #-------------------------------------------------------------------------# #
#                                                                             #
# plotSpecIPQU                                                                #
# plotSpecRMS                                                                 #
# plotPolang                                                                  #
# plotFracPol                                                                 #
# plotFracQvsU                                                                #
# plotPolsummary                                                              #
# plotRMSF                                                                    #
# plotDirtyFDF                                                                #
# plotCleanFDF                                                                #
# plotStampI                                                                  #
# plotStampP                                                                  #
# plotSctHstQuery                                                             #
#                                                                             #
#=============================================================================#
#                                                                             #
# The MIT License (MIT)                                                       #
#                                                                             #
# Copyright (c) 2015 Cormac R. Purcell                                        #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the "Software"),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
#=============================================================================#

import os
import sys
import math as m
import numpy as np
import StringIO
import traceback
import astropy.io.fits as pf
import matplotlib as mpl
from matplotlib.ticker import MaxNLocator
from matplotlib.patches import Polygon
from matplotlib.ticker import FuncFormatter
from matplotlib.figure import Figure

from util_plotFITS import plot_fits_map
from util_PPC import xfloat
from normalize import APLpyNormalize

# Alter the default linewidths etc.
mpl.rcParams['lines.linewidth'] = 1.0
mpl.rcParams['axes.linewidth'] = 0.8
mpl.rcParams['xtick.major.size'] = 8.0
mpl.rcParams['xtick.minor.size'] = 4.0
mpl.rcParams['ytick.major.size'] = 8.0
mpl.rcParams['ytick.minor.size'] = 4.0
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.size'] = 12.0

# Constants
C = 2.99792458e8


#-----------------------------------------------------------------------------#
def xint(x, default=None):
    if x is None:
        return default
    return int(x)

 
#-----------------------------------------------------------------------------#
def filter_range_indx(a, dataMin=None, dataMax=None, filterNans=False):
    """Return a boolean array [True, ...] where data falls outside of the
    range [dataMin <= a <= dataMax]."""

    if filterNans:
        iNaN = np.zeros_like(a, dtype="bool")
    else:
        iNaN = a!=a
    if dataMin is None:
        i1 = np.ones_like(a, dtype="bool")
    else:
        i1 = a>=dataMin
        i1+=iNaN
    if dataMax is None:
        i2 = np.ones_like(a, dtype="bool")
    else:
        i2 = a<=dataMax
        i2+=iNaN
    return ~i1+~i2


#-----------------------------------------------------------------------------#
def tweakAxFormat(ax, pad=10, loc='upper right', linewidth=1, ncol=1,
                   bbox_to_anchor=(1.00, 1.00), showLeg=True):
    
    # Axis/tic formatting
    ax.tick_params(pad=pad)
    for line in ax.get_xticklines() + ax.get_yticklines():
        line.set_markeredgewidth(linewidth)
        
    # Legend formatting
    if showLeg:
        leg = ax.legend(numpoints=1, loc=loc, shadow=False,
                        borderaxespad=0.3, ncol=ncol,
                        bbox_to_anchor=bbox_to_anchor)
        for t in leg.get_texts():
            t.set_fontsize('small') 
        leg.get_frame().set_linewidth(0.5)
        leg.get_frame().set_alpha(0.5)
    return ax
    

#-----------------------------------------------------------------------------#
def format_ticks(ax, pad=10, w=1.0):

    ax.tick_params(pad=pad)
    for line in ax.get_xticklines() + ax.get_yticklines():
        line.set_markeredgewidth(w)


#-----------------------------------------------------------------------------#
def plot_I_vs_nu_ax(ax, freqArr_Hz, IArr_mJy, dIArr_mJy=None,
                    freqHirArr_Hz=None, IModArr_mJy=None, axisYright=False,
                    axisXtop=False):
    """Plot the I spectrum and an optional model."""

    # Set the axis positions
    if axisYright:
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")
    if axisXtop:
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position("top")            

    # Default to non-high-resolution inputs
    if freqHirArr_Hz is None:
        freqHirArr_Hz =  freqArr_Hz

    # Plot I versus frequency
    ax.errorbar(x=freqArr_Hz/1e9, y=IArr_mJy, yerr=dIArr_mJy, mfc='none',
                ms=4, fmt='D', ecolor='grey', elinewidth=1.0, capsize=2,
                label='Stokes I')
    if IModArr_mJy is not None:
        ax.plot(freqHirArr_Hz/1e9, IModArr_mJy, color='k', lw=0.5,
                label='I Model')
    #ax.text(0.05, 0.94, 'Stokes I Spectrum', transform=ax.transAxes)
    ax.yaxis.set_major_locator(MaxNLocator(4))
    ax.xaxis.set_major_locator(MaxNLocator(5))

    xRange = (np.nanmax(freqArr_Hz)-np.nanmin(freqArr_Hz))/1e9 
    ax.set_xlim( np.min(freqArr_Hz)/1e9 - xRange*0.05,
                 np.max(freqArr_Hz)/1e9 + xRange*0.05)
    
    # Formatting
    ax.yaxis.set_major_locator(MaxNLocator(4))
    ax.xaxis.set_major_locator(MaxNLocator(4))
    xRange = (np.nanmax(freqArr_Hz)-np.nanmin(freqArr_Hz))/1e9 
    ax.set_xlim( np.min(freqArr_Hz)/1e9 - xRange*0.05,
                 np.max(freqArr_Hz)/1e9 + xRange*0.05)
    ax.set_xlabel('$\\nu$ (GHz)')
    ax.set_ylabel('Flux Density (mJy)')

    # Format tweaks
    ax = tweakAxFormat(ax)
    ax.relim()
    ax.autoscale_view()


#-----------------------------------------------------------------------------#
def plot_PQU_vs_nu_ax(ax, freqArr_Hz, QArr_mJy, UArr_mJy, dQArr_mJy=None,
                      dUArr_mJy=None, freqHirArr_Hz=None, QmodArr=None,
                      UmodArr=None, axisYright=False, axisXtop=False):
    """Plot the P, Q & U spectrum and an optional model. """

    # Set the axis positions
    if axisYright:
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")
    if axisXtop:
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position("top")    
        
    # Default to non-high-resolution inputs
    if freqHirArr_Hz is None:
        freqHirArr_Hz =  freqArr_Hz

    # Calculate P and errors
    PArr_mJy = np.sqrt(np.power(QArr_mJy,2) + np.power(UArr_mJy,2))
    if dQArr_mJy is None or dUArr_mJy is None:
        dPArr_mJy = None
    else:
        dPArr_mJy = np.sqrt(np.power(dQArr_mJy,2) + np.power(dUArr_mJy,2))
        
    # Plot P, Q, U versus frequency
    ax.errorbar(x=freqArr_Hz/1e9, y=QArr_mJy, yerr=dQArr_mJy, mec='b',
                mfc='none', ms=4, fmt='D', color='g', elinewidth=1.0,
                capsize=2, label='Stokes Q')
    ax.errorbar(x=freqArr_Hz/1e9, y=UArr_mJy, yerr=dUArr_mJy, mec='r',
                mfc='none', ms=4, fmt='D', color='r', elinewidth=1.0,
                capsize=2, label='Stokes U')
    ax.errorbar(x=freqArr_Hz/1e9, y=PArr_mJy, yerr=dPArr_mJy, mec='k',
                mfc='none', ms=4, fmt='D', color='k', elinewidth=1.0,
                capsize=2, label='Intensity P')
    
    # Plot the models
    if QmodArr is not None:
        ax.plot(freqHirArr_Hz/1e9, QmodArr, color='b', lw=0.5, label='Model Q')
    if UmodArr is not None:
        ax.plot(freqHirArr_Hz/1e9, UmodArr, color='r', lw=0.5, label='Model U')
    if QmodArr is not None and UmodArr is not None:
        PmodArr = np.sqrt(QmodArr**2.0 + UmodArr**2.0 )
        ax.plot(freqHirArr_Hz/1e9, PmodArr, color='k', lw=0.5, label='Model P')

    # Formatting
    ax.yaxis.set_major_locator(MaxNLocator(4))
    ax.xaxis.set_major_locator(MaxNLocator(4))
    xRange = (np.nanmax(freqArr_Hz)-np.nanmin(freqArr_Hz))/1e9 
    ax.set_xlim( np.min(freqArr_Hz)/1e9 - xRange*0.05,
                 np.max(freqArr_Hz)/1e9 + xRange*0.05)
    ax.set_xlabel('$\\nu$ (GHz)')
    ax.set_ylabel('Flux Density (mJy)')
    ax.axhline(0, color='grey')

    # Format tweaks
    ax = tweakAxFormat(ax)
    ax.relim()
    ax.autoscale_view()


#-----------------------------------------------------------------------------#
def plot_rmsIQU_vs_nu_ax(ax, freqArr_Hz, rmsIArr_mJy,  rmsQArr_mJy,
                         rmsUArr_mJy, axisYright=False, axisXtop=False):
    """Plot the noise spectra in Stokes I, Q & U. """

    # Set the axis positions
    if axisYright:
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")
    if axisXtop:
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position("top")    
        
    # Plot the rms spectra in GHz and mJy
    ax.plot(freqArr_Hz/1e9, rmsIArr_mJy, marker='o', color='k', lw='0.5',
             label='rms I')
    ax.plot(freqArr_Hz/1e9, rmsQArr_mJy, marker='o', color='b', lw='0.5',
             label='rms Q')
    ax.plot(freqArr_Hz/1e9, rmsUArr_mJy, marker='o', color='r', lw='0.5',
             label='rms U')
    #ax.text(0.05, 0.94, 'I, Q & U RMS', transform=ax.transAxes)

    # Formatting
    ax.yaxis.set_major_locator(MaxNLocator(4))
    ax.xaxis.set_major_locator(MaxNLocator(4))
    xRange = (np.nanmax(freqArr_Hz)-np.nanmin(freqArr_Hz))/1e9 
    ax.set_xlim( np.min(freqArr_Hz)/1e9 - xRange*0.05,
                 np.max(freqArr_Hz)/1e9 + xRange*0.05)
    ax.set_xlabel('$\\nu$ (GHz)')
    ax.set_ylabel('Flux Density (mJy bm$^{-1}$)')
    
    # Format tweaks
    ax = tweakAxFormat(ax)
    ax.relim()
    ax.autoscale_view()


#-----------------------------------------------------------------------------#
def plot_pqu_vs_lamsq_ax(ax, lamSqArr_m2, qArr, uArr, dqArr=None, duArr=None,
                         lamSqHirArr_m2=None, qModArr=None, uModArr=None,
                         axisYright=False, axisXtop=False):

    # Set the axis positions
    if axisYright:
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")
    if axisXtop:
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position("top")    
        
    # Default to non-high-resolution inputs
    if lamSqHirArr_m2 is None:
        lamSqHirArr_m2 = lamSqArr_m2

    # Calculate p and errors
    pArr = np.sqrt(qArr**2.0 + uArr**2.0 )
    if dqArr is None or duArr is None:
        dpArr = None
    else:
        dpArr = np.sqrt(dqArr**2.0 + duArr**2.0 )
        
    # Plot p, q, u versus lambda^2
    ax.errorbar(x=lamSqArr_m2, y=qArr, yerr=dqArr, mec='b', mfc='none', ms=4,
                fmt='D', ecolor='b', elinewidth=1.0, capsize=2,
                label='Stokes q')
    ax.errorbar(x=lamSqArr_m2, y=uArr, yerr=duArr, mec='r', mfc='none', ms=4,
                fmt='D', ecolor='r', elinewidth=1.0, capsize=2,
                label='Stokes u')
    ax.errorbar(x=lamSqArr_m2, y=pArr, yerr=dpArr, mec='k', mfc='none', ms=4,
                fmt='D', ecolor='k', elinewidth=1.0, capsize=2,
                label='Intensity p')

    # Plot the models
    if qModArr is not None:
        ax.plot(lamSqHirArr_m2, qModArr, color='b', lw=0.5, label='Model q')
    if uModArr is not None:
        ax.plot(lamSqHirArr_m2, uModArr, color='r', lw=0.5, label='Model u')
    if qModArr is not None and uModArr is not None:
        pModArr = np.sqrt(qModArr**2.0 + uModArr**2.0 )
        ax.plot(lamSqHirArr_m2, pModArr, color='k', lw=0.5, label='Model p')

    # Formatting
    ax.yaxis.set_major_locator(MaxNLocator(4))
    ax.xaxis.set_major_locator(MaxNLocator(4))
    xRange = np.nanmax(lamSqArr_m2)-np.nanmin(lamSqArr_m2)
    ax.set_xlim( np.min(lamSqArr_m2) - xRange*0.05,
                 np.max(lamSqArr_m2) + xRange*0.05)
    ax.set_xlabel('$\\lambda^2$ (m$^2$)')
    ax.set_ylabel('Fractional Polarisation')
    ax.axhline(0, color='grey')

    # Format tweaks
    ax = tweakAxFormat(ax)
    ax.relim()
    ax.autoscale_view()


#-----------------------------------------------------------------------------#
def plot_psi_vs_lamsq_ax(ax, lamSqArr_m2, qArr, uArr, dqArr=None, duArr=None,
                         lamSqHirArr_m2=None, qModArr=None, uModArr=None,
                         axisYright=False, axisXtop=False):

    # Set the axis positions
    if axisYright:
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")
    if axisXtop:
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position("top")    
        
    # Default to non-high-resolution inputs
    if lamSqHirArr_m2 is None:
        lamSqHirArr_m2 = lamSqArr_m2

    # Plot psi versus lambda^2
    pArr = np.sqrt(qArr**2.0 + uArr**2.0 )
    psiArr_deg = np.degrees( np.arctan2(uArr, qArr) / 2.0 )
    if dqArr is None or duArr is None:
        dQUArr = None
        dPsiArr_deg = None
    else:
        dQUArr = np.sqrt(dqArr**2.0 + duArr**2.0)
        dPsiArr_deg = np.degrees( np.sqrt( (qArr * duArr)**2.0 +
                                           (uArr * dqArr)**2.0) /
                                  (2.0 * pArr**2.0) )
    ax.errorbar(x=lamSqArr_m2, y=psiArr_deg, yerr=dPsiArr_deg, mec='k',
                 mfc='none', ms=4, fmt='D', ecolor='k', elinewidth=1.0,
                 capsize=2)
    if qModArr is not None and uModArr is not None:
        psiHirArr_deg = np.degrees( np.arctan2(uModArr, qModArr) / 2.0 )
        ax.plot(lamSqHirArr_m2, psiHirArr_deg, color='k', lw=0.5,
                 label='Model $\psi$')
    ax.set_ylim(-99.9, 99.9)
    ax.axhline(0, color='grey')

    # Formatting
    ax.yaxis.set_major_locator(MaxNLocator(4))
    ax.xaxis.set_major_locator(MaxNLocator(4))
    xRange = np.nanmax(lamSqArr_m2)-np.nanmin(lamSqArr_m2)
    ax.set_xlim( np.min(lamSqArr_m2) - xRange*0.05,
                 np.max(lamSqArr_m2) + xRange*0.05)
    ax.set_xlabel('$\\lambda^2$ (m$^2$)')
    ax.set_ylabel('$\psi$ (degrees)')
    
    # Format tweaks
    ax = tweakAxFormat(ax, showLeg=False)
    ax.relim()
    ax.autoscale_view()
    

#-----------------------------------------------------------------------------#
def plot_q_vs_u_ax(ax, lamSqArr_m2, qArr, uArr, dqArr=None, duArr=None,
                   lamSqHirArr_m2=None, qModArr=None, uModArr=None,
                   axisYright=False, axisXtop=False):

    # Set the axis positions
    if axisYright:
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")
    if axisXtop:
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position("top")    
        
    # Default to non-high-resolution inputs
    if lamSqHirArr_m2 is None:
        lamSqHirArr_m2 = lamSqArr_m2

    # Plot U versus Q 
    ax.errorbar(x=uArr, y=qArr, xerr=duArr, yerr=dqArr, mec='k', mfc='none',
                 ms=4, fmt='D', ecolor='k', elinewidth=1.0, capsize=2)
    if qModArr is not None and uModArr is not None:
        ax.plot(uModArr, qModArr, color='k', lw=0.5, label='Model q & u')
    ax.axhline(0, color='grey')
    ax.axvline(0, color='grey')
    
    # Formatting
    ax.yaxis.set_major_locator(MaxNLocator(4))
    ax.xaxis.set_major_locator(MaxNLocator(4))
    ax.set_xlabel('Stokes u')
    ax.set_ylabel('Stokes q')
    format_ticks(ax, 10, 1.2)
    ax.relim()
    ax.autoscale_view(False,True,True)


#-----------------------------------------------------------------------------#
def plot_RMSF_ax(ax, phiArr, RMSFArr, fwhmRMSF=None, axisYright=False,
                 axisXtop=False):

    # Set the axis positions
    if axisYright:
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")
    if axisXtop:
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position("top")    
        
    # Plot the RMSF
    ax.step(phiArr, RMSFArr.real, where='mid', color='b', lw='0.5',
            label='Real')
    ax.step(phiArr, RMSFArr.imag, where='mid', color='r', lw='0.5',
            label='Imaginary')
    ax.step(phiArr, np.abs(RMSFArr) , where='mid', color='k', lw='1.0',
            label='PI')
    #ax.text(0.05, 0.94, 'RMSF', transform=ax.transAxes)

    # Plot the Gaussian fit
    if fwhmRMSF is not None:
        yGauss = gauss([1.0, 0.0, fwhmRMSF])(phiArr)
        ax.plot(phiArr, yGauss, color='magenta',marker='None',mfc='w',
                mec='g', ms=10, label='Gaussian Fit', lw=2.0, ls='--')
    
    # Scaling
    ax.yaxis.set_major_locator(MaxNLocator(4))
    ax.xaxis.set_major_locator(MaxNLocator(4))
    xRange = np.nanmax(phiArr)-np.nanmin(phiArr)
    ax.set_xlim( np.min(phiArr) - xRange*0.01,
                 np.max(phiArr) + xRange*0.01)
    ax.set_ylabel('Normalised Units')
    ax.set_xlabel('$\phi$ rad m$^{-2}$')
    ax.axhline(0, color='grey')

    # Format tweaks
    ax = tweakAxFormat(ax)
    ax.relim()
    ax.autoscale_view()


#-----------------------------------------------------------------------------
def gauss(p):
    """Return a fucntion to evaluate a Gaussian with parameters
    p = [amp, mean, FWHM]"""
    
    a, b, w = p
    gfactor = 2.0 * m.sqrt(2.0 * m.log(2.0))
    s = w / gfactor
    
    def rfunc(x):
        y = a * np.exp(-(x-b)**2.0 /(2.0 * s**2.0))
        return y
    
    return rfunc
    

#-----------------------------------------------------------------------------#
def plot_dirtyFDF_ax(ax, phiArr, FDFArr_mJy, gaussParm=[], title="Dirty FDF",
                     axisYright=False, axisXtop=False):

    # Set the axis positions
    if axisYright:
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")
    if axisXtop:
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position("top")
        
    # Plot the FDF
    FDFpiArr_mJy = np.sqrt( np.power(FDFArr_mJy.real, 2.0) +
                            np.power(FDFArr_mJy.imag, 2.0) )
    ax.step(phiArr, FDFArr_mJy.real, where='mid', color='b', lw='0.5',
            label='Real')
    ax.step(phiArr, FDFArr_mJy.imag, where='mid', color='r', lw='0.5',
            label='Imaginary')
    ax.step(phiArr, FDFpiArr_mJy, where='mid', color='k', lw='1.0',
            label='PI')
    #ax.text(0.05, 0.94, title, transform=ax.transAxes)

    # Plot the Gaussian peak
    if len(gaussParm)==3:
        # [amp, mean, FWHM]
        phiTrunkArr = np.where(phiArr>=gaussParm[1]-gaussParm[2]/3.0,
                               phiArr, np.nan)
        phiTrunkArr = np.where(phiArr<=gaussParm[1]+gaussParm[2]/3.0,
                               phiTrunkArr, np.nan)
        yGauss = gauss(gaussParm)(phiTrunkArr)
        ax.plot(phiArr, yGauss, color='magenta',marker='None',mfc='w',
                mec='g', ms=10, label='Peak Fit', lw=2.5, ls='-')

    # Scaling
    ax.yaxis.set_major_locator(MaxNLocator(4))
    ax.xaxis.set_major_locator(MaxNLocator(4))
    xRange = np.nanmax(phiArr)-np.nanmin(phiArr)
    ax.set_xlim( np.min(phiArr) - xRange*0.01,
                 np.max(phiArr) + xRange*0.01)
    ax.set_ylabel('Flux Density (mJy)')
    ax.set_xlabel('$\phi$ rad m$^{-2}$')
    ax.axhline(0, color='grey')

    # Format tweaks
    ax = tweakAxFormat(ax)
    ax.relim()
    ax.autoscale_view()


#-----------------------------------------------------------------------------#
def plot_cleanFDF_ax(ax, phiArr, cleanFDFArr_mJy, ccFDFArr_mJy=None,
                     dirtyFDFArr_mJy=None, gaussParm=[], title="Clean FDF",
                     cutoff_mJy=None, axisYright=False, axisXtop=False):

    # Set the axis positions
    if axisYright:
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")
    if axisXtop:
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position("top")    
        
    # Plot the dirty FDF in the background
    if not dirtyFDFArr_mJy is None:
        dirtyFDFpiArr_mJy = np.sqrt( np.power(dirtyFDFArr_mJy.real, 2.0) +
                                     np.power(dirtyFDFArr_mJy.imag, 2.0) )
        
        ax.step(phiArr, dirtyFDFpiArr_mJy, where='mid', color='grey',
                 lw='0.5', label='Dirty')

    # Plot the clean FDF
    cleanFDFpiArr_mJy = np.sqrt( np.power(cleanFDFArr_mJy.real, 2.0) +
                                 np.power(cleanFDFArr_mJy.imag, 2.0) )
    ax.step(phiArr, cleanFDFArr_mJy.real, where='mid', color='b', lw='0.5',
            label='Real')
    ax.step(phiArr, cleanFDFArr_mJy.imag, where='mid', color='r', lw='0.5',
            label='Imaginary')
    ax.step(phiArr, cleanFDFpiArr_mJy, where='mid', color='k', lw='1.0',
            label='PI')
    #ax.text(0.05, 0.94, title, transform=ax.transAxes)

    # Plot the CC spectrum
    if not ccFDFArr_mJy is None:
        ax.step(phiArr, ccFDFArr_mJy, where='mid', color='g', lw='0.5',
                label='CC')

    # Plot the Gaussian peak
    if len(gaussParm)==3:
        # [amp, mean, FWHM]
        phiTrunkArr = np.where(phiArr>=gaussParm[1]-gaussParm[2]/3.0,
                               phiArr, np.nan)
        phiTrunkArr = np.where(phiArr<=gaussParm[1]+gaussParm[2]/3.0,
                               phiTrunkArr, np.nan)
        yGauss = gauss(gaussParm)(phiTrunkArr)
        ax.plot(phiArr, yGauss, color='magenta',marker='None',mfc='w',
                mec='g', ms=10, label='Peak Fit', lw=2.5, ls='-')

    # Plot the clean cutoff line
    if not cutoff_mJy is None:
        ax.axhline(cutoff_mJy, color="r", ls='--')

    # Scaling
    ax.yaxis.set_major_locator(MaxNLocator(4))
    ax.xaxis.set_major_locator(MaxNLocator(4))
    xRange = np.nanmax(phiArr)-np.nanmin(phiArr)
    ax.set_xlim( np.min(phiArr) - xRange*0.01,
                 np.max(phiArr) + xRange*0.01)
    ax.set_ylabel('Flux Density (mJy)')
    ax.set_xlabel('$\phi$ rad m$^{-2}$')
    ax.axhline(0, color='grey')

    # Format tweaks
    ax = tweakAxFormat(ax)
    ax.relim()
    ax.autoscale_view()


#-----------------------------------------------------------------------------#
def plot_hist4_ax(ax, popLst, nBins=10, doXlog=False, doYlog=False, styIndx=0,
                  xMin=None, xMax=None, yMax=None, xLabel="", yLabel="",
                  title="", legLabLst=[], legLoc="tr", verbose=False):
    
    # Format of the histogram lines and shading.
    # Up to four histograms are supported and two alternative styles
    edgeColourLst = [['black', 'red', 'blue', 'black'],
                     ['grey','red', 'blue', 'black']]
    fillColourLst = [['#dddddd', 'none', 'none', 'none'],
                     ['none', 'none', 'none', 'none']]
    hatchStyleLst = [['', '/', '\\', ''],
                     ['/', '\\', '', '']]
    histLinewidthLst = [[1.0,1.0,1.0,1.5],
                        [1.0,1.0,1.0,1.5]]
    
    # Translate the legend location code
    if legLoc not in ["tl", "tr", "bl", "br"]: legLoc = "tr"
    locTab = {"tl": "upper left", "tr": "upper right",
              "bl": "lower left", "br": "lower right"}
    legLoc = locTab[legLoc]

    # TODO: Remove extra columns in the recarrays


    # Determine the max and min of the ensemble population
    popEnsemble = np.concatenate(popLst).astype(np.float)
    xMinData = float(np.nanmin(popEnsemble))      
    xMaxData = float(np.nanmax(popEnsemble))

    # All valid data must have the same sign for log plots
    if doXlog:
        if not (xMinData<0) == (xMaxData<0):
            print "\nErr: for log axis all data must have the same sign!"
            return
    sign = np.sign(popEnsemble)[0]
    
    # Calculate the bin edges
    if doXlog:
        logBins = np.linspace(m.log10(abs(xMinData)), m.log10(abs(xMaxData)),
                              int(nBins+1))
        b = np.power(10.0, logBins) * sign
    else:
        b = np.linspace(xMinData, xMaxData, int(nBins+1))

    # Bin the data in each population
    nLst = []
    for p in popLst:
        n, b =  np.histogram(p.astype(np.float), bins = b)
        n = np.array(n, dtype=np.float)
        nLst.append(n)

    # Print the binned values to the screen
    if verbose:
        print "\n#BIN, COUNTS ..."
        binCentreArr = b[:-1]+np.diff(b)/2.0
        for i in range(len(binCentreArr)):
            print binCentreArr[i],
            for j in range(len(nLst)):
                print nLst[j][i],
            print
    
    # Set the Y-axis limits
    nEnsemble = np.concatenate(nLst)
    if doYlog:
        yZeroPt = 0.8
        yMin = yZeroPt
        if yMax is None:
            yMaxData = float(max(nEnsemble))
            yFac = abs(yMaxData/yZeroPt)
            yMax = yMaxData*(1+ m.log10(yFac)*0.3)
    else:
        yZeroPt = 0.0
        yMin = yZeroPt
        if yMax is None:
            yMax = float(max(nEnsemble))*1.2
        
    # Set the X-axis limits, incorporating a single padding bin
    xFac = (len(b)-1)*0.05
    if doXlog:
        sign = np.sign(b)[0]
        logBins = np.log10(b*sign)
        logBinWidth = np.max(np.diff(logBins))
        if xMin is None:
            xMin = 10**(logBins[0] - logBinWidth*xFac)*sign
        if xMax is None:
            xMax = 10**(logBins[-1] + logBinWidth*xFac)*sign
    else:        
        linBinWidth = np.max(np.diff(b))
        if xMin is None:
            xMin = b[0] - linBinWidth*xFac
        if xMax is None:
            xMax = b[-1] + linBinWidth*xFac
        
    # Set the axis formatter for log scale axes
    if doXlog:
        ax.set_xscale('symlog')
        majorFormatterX = FuncFormatter(label_format_exp(5.0))
        ax.xaxis.set_major_formatter(majorFormatterX)
    if doYlog:
        ax.set_yscale('symlog')
        majorFormatterY = FuncFormatter(label_format_exp(3.0))
        ax.yaxis.set_major_formatter(majorFormatterY)
    
    # Create individual histogram polygons. Manually creating histograms gives
    # more control than inbuilt matplotlib function - when originally writing
    # this code the fill styles did not work well.
    for i in range(len(nLst)):
        
        # Set the legend labels
        try:
            legLabel = legLabLst[i]
            if legLabLst[i]=="":
                raise Exception
        except Exception:
            legLabel = "Query %s" % (i+1)

            
        # Create the histograms from line-segments
        polyCoords = mk_hist_poly(b, nLst[i], doYlog, zeroPt=0.7)
        hist = Polygon(polyCoords, closed=True, animated=False, linewidth=2.7,
                       label=legLabel)
        hist.set_linewidth(histLinewidthLst[styIndx][i])
        hist.set_edgecolor(edgeColourLst[styIndx][i])
        hist.set_facecolor(fillColourLst[styIndx][i])
        hist.set_hatch(hatchStyleLst[styIndx][i])
        ax.add_patch(hist)

    # Set the X axis limits
    ax.set_xlim(xMin, xMax)
    ax.set_ylim(yMin, yMax)
    
    # Draw the labels on the plot 
    ax.set_xlabel(xLabel)
    ax.set_ylabel(yLabel)
    ax.set_title(title, size=14)

    # Format tweaks
    tweakAxFormat(ax, showLeg=True, loc=legLoc)


#-----------------------------------------------------------------------------#
def plot_scatter4_ax(ax, popLst, doXlog=False, doYlog=False, zPower=1.0,
                     styIndx=0, xMin=None, xMax=None, yMin=None,
                     yMax=None, zMin=None, zMax=None, xLabel="", yLabel="",
                     zLabel="", title="", legLabLst=[], showCbar=False,
                     show11Line=False, legLoc="tr", verbose=False):

    # Format of the scatter points and shading.
    # Up to four populations are supported and two alternative styles
    edgeColourLst = [['black', 'black', 'black', 'black'],
                     ['black','red', 'green', 'blue']]
    fillColourLst = [['black', 'red', 'green', 'blue'],
                     ['none', 'none', 'none', 'none']]
    symbolLst =     [['o', 's', 'd', '^'],
                     ['o', '+', 's', 'd']]
    symbolSize =    [[45, 30, 30, 30],
                     [45, 80, 30, 30]]

    # Translate the legend location code
    if legLoc not in ["tl", "tr", "bl", "br"]: legLoc = "tr"
    locTab = {"tl": "upper left", "tr": "upper right",
              "bl": "lower left", "br": "lower right"}
    legLoc = locTab[legLoc]

    # Separate out the X Y and Z data
    xLst = []
    yLst = []
    zLst = []
    for i in range(len(popLst)):
        colNames = popLst[i].dtype.names
        nCols = len(colNames)
        xLst.append(popLst[i][colNames[0]])
        yLst.append(popLst[i][colNames[1]])
        if nCols>2:
            yLst.append(popLst[i][colNames[2]])

    # Determine the max and min of the ensemble population
    xEnsemble = np.concatenate(xLst).astype(np.float)
    signX = np.sign(xEnsemble)[0]
    xMinData = float(np.nanmin(xEnsemble))
    xMaxData = float(np.nanmax(xEnsemble))
    yEnsemble = np.concatenate(yLst).astype(np.float)
    signY = np.sign(yEnsemble)[0]
    yMinData = float(np.nanmin(yEnsemble))
    yMaxData = float(np.nanmax(yEnsemble))
    if not zLst==[]:
        zEnsemble = np.concatenate(zLst).astype(np.float)
        signZ = np.sign(zEnsemble)[0]
        zMinData = float(np.nanmin(zEnsemble))
        zMaxData = float(np.nanmax(zEnsemble))

    # All valid data must have the same sign for log plots
    if doXlog:
        if not (xMinData<0) == (xMaxData<0):
            print "\nErr: for log X-axis all data must have the same sign!"
            sys.exit()
    if doYlog:
        if not (yMinData<0) == (yMaxData<0):
            print "\nErr: for log Y-axis all data must have the same sign!"
            sys.exit()
    if zLst is not None and zPower!=1.0:
        if not (zMinData<0) == (zMaxData<0):
            print "\nErr: for log Z-axis all data must have the same sign!"
            sys.exit()

    # Set the plotting ranges (& colour limits)
    if doXlog:
        xFac = abs(xMaxData/xMinData)
        if xMin is None:
            xMin = xMinData/(1+ m.log10(xFac)*0.1)
        if xMax is None:
            xMax = xMaxData*(1+ m.log10(xFac)*0.1)
    else:
        xPad = abs(xMaxData - xMinData) * 0.04
        if xMin is None:
            xMin = xMinData - xPad
        if xMax is None:
            xMax = xMaxData + xPad
    if doYlog:
        yFac = abs(yMaxData/yMinData)
        if yMin is None:
            yMin = yMinData/(1+ m.log10(yFac)*0.1)
        if yMax is None:
            yMax = yMaxData*(1+ m.log10(yFac)*0.1)
    else:
        yPad = abs(yMaxData - yMinData) * 0.05
        if yMin is None:
            yMin = yMinData - yPad
        if yMax is None:
            yMax = yMaxData + yPad

    # Set the z-colour range
    if not zLst==[]:
        if not np.all(np.isnan(zEnsemble)):
            if zMin is None:
                zMin = zMinData
            if zMax is None:
                zMax = zMaxData

    # Set the axis formatter for log scale axes
    if doXlog:
        ax.set_xscale('log')
        majorFormatterX = FuncFormatter(label_format_exp(5.0))
        ax.xaxis.set_major_formatter(majorFormatterX)
    if doYlog:
        ax.set_yscale('log')
        majorFormatterY = FuncFormatter(label_format_exp(3.0))
        ax.yaxis.set_major_formatter(majorFormatterY)
    norm = APLpyNormalize(stretch='power', exponent=zPower,
                                vmin=zMin, vmax=zMax)

    # Plot each set of points in turn
    sc3D = None
    zMap = 'r'
    for i in range(len(xLst)):

        # Map the z axis to the colours
        if not zLst==[]:
            if np.all(np.isnan(zLst[i])):
                zMap = fillColourLst[styIndx][i]
            else:
                zMap = zLst[i]

        # Set the legend labels
        try:
            legLabel = legLabLst[i]
            if legLabLst[i]=="":
                raise Exception
        except Exception:
            legLabel = "Query %s" % (i+1)

        # Add the points to the plot
        sc = ax.scatter(xLst[i], yLst[i],
                        marker = symbolLst[styIndx][i],
                        s = symbolSize[styIndx][i],
                        c = zMap,
                        norm = norm,
                        vmin = zMin,
                        vmax = zMax,
                        linewidths = 0.5,
                        edgecolors=edgeColourLst[styIndx][i],
                        label=legLabel)
        if not zLst==[]:
            if not np.all(np.isnan(zLst[i])):
                sc3D = sc
        
    # Set the X axis limits
    ax.set_xlim(xMin, xMax)
    ax.set_ylim(yMin, yMax)
    
    # Draw the labels on the plot 
    ax.set_xlabel(xLabel)
    ax.set_ylabel(yLabel)
    ax.set_title(title)

    # Format tweaks
    tweakAxFormat(ax, showLeg=True, loc=legLoc)

    



    
# Axis code above
#=============================================================================#
# Figure code below


#-----------------------------------------------------------------------------#
def plotSpecIPQU(dataMan, indx, io='fig'):

    # Get the data
    freqArr_Hz, IArr_Jy, rmsIArr_Jy= dataMan.get_specI_byindx(indx)
    dummy, QArr_Jy, rmsQArr_Jy= dataMan.get_specQ_byindx(indx)
    dummy, UArr_Jy, rmsUArr_Jy= dataMan.get_specU_byindx(indx)
    
    # Get the models to overplot
    freqHirArr_Hz, qModArr, uModArr = dataMan.get_thin_qumodel_byindx(indx, 
                                                             oversample=True)
    freqHirArr_Hz, IModArr_mJy = dataMan.get_modI_byindx(indx, oversample=True)
    QmodArr = qModArr * IModArr_mJy
    UmodArr = uModArr * IModArr_mJy
    
    # Setup the figure
    fig = Figure()
    fig.set_size_inches([8,8])
    
    # Plot I versus nu,
    ax1 = fig.add_subplot(211)
    plot_I_vs_nu_ax(ax=ax1,
                    freqArr_Hz  = freqArr_Hz,
                    IArr_mJy    = IArr_Jy*1e3,
                    dIArr_mJy   = rmsIArr_Jy*1e3,
                    freqHirArr_Hz = freqHirArr_Hz,
                    IModArr_mJy = IModArr_mJy*1e3)
    ax1.set_xlabel('')
    [label.set_visible(False) for label in ax1.get_xticklabels()]

    # Plot Stokes P, Q & U
    ax2 = fig.add_subplot(212, sharex=ax1)
    plot_PQU_vs_nu_ax(ax=ax2,
                      freqArr_Hz = freqArr_Hz,
                      QArr_mJy   = QArr_Jy*1e3,
                      UArr_mJy   = UArr_Jy*1e3,
                      dQArr_mJy  = rmsQArr_Jy*1e3,
                      dUArr_mJy  = rmsUArr_Jy*1e3,
                      freqHirArr_Hz=freqHirArr_Hz,
                      QmodArr=QmodArr*1e3,
                      UmodArr=UmodArr*1e3)
    
    # Write to the pipe
    if io=='string':
        sio = StringIO.StringIO()
        setattr(sio, "name", "foo.jpg")
        fig.savefig(sio, format='jpg' )    
        return sio
    else:
        return fig
    

#-----------------------------------------------------------------------------#
def plotSpecRMS(dataMan, indx, io='fig'):

    # Get the data
    freqArr_Hz, IArr_Jy, rmsIArr_Jy= dataMan.get_specI_byindx(indx)
    dummy, QArr_Jy, rmsQArr_Jy= dataMan.get_specQ_byindx(indx)
    dummy, UArr_Jy, rmsUArr_Jy= dataMan.get_specU_byindx(indx)
    
    # Setup the figure
    fig = Figure()
    fig.set_size_inches([8,8])

    # Plot Stokes I, Q & U
    ax1 = fig.add_subplot(111)
    plot_rmsIQU_vs_nu_ax(ax=ax1,
                         freqArr_Hz  = freqArr_Hz,
                         rmsIArr_mJy = rmsIArr_Jy*1e3,
                         rmsQArr_mJy = rmsQArr_Jy*1e3,
                         rmsUArr_mJy = rmsUArr_Jy*1e3)
    
    # Write to the pipe
    if io=='string':
        sio = StringIO.StringIO()
        setattr(sio, "name", "foo.jpg")
        fig.savefig(sio, format='jpg' )    
        return sio
    else:
        return fig
    
#-----------------------------------------------------------------------------#
def plotPolang(dataMan, indx, io='fig'):

    # Get the data
    freqArr_Hz, IArr_Jy, rmsIArr_Jy= dataMan.get_specI_byindx(indx)
    dummy, modIArr_Jy = dataMan.get_modI_byindx(indx)
    dummy, QArr_Jy, rmsQArr_Jy= dataMan.get_specQ_byindx(indx)
    dummy, UArr_Jy, rmsUArr_Jy= dataMan.get_specU_byindx(indx)
    
    # Calculate fractional polarisation spectra
    qArr = QArr_Jy / modIArr_Jy
    uArr = UArr_Jy / modIArr_Jy
    dqArr = qArr * np.sqrt( (rmsQArr_Jy/QArr_Jy)**2.0 + 
                            (rmsIArr_Jy/IArr_Jy)**2.0 )
    duArr = uArr * np.sqrt( (rmsUArr_Jy/UArr_Jy)**2.0 + 
                            (rmsIArr_Jy/IArr_Jy)**2.0 )
    lamSqArr_m2 = np.power(C/freqArr_Hz, 2.0)

    # Get the models to overplot
    freqHirArr_Hz, qModArr, uModArr = dataMan.get_thin_qumodel_byindx(indx, 
                                                             oversample=True)
    lamSqHirArr_m2 = np.power(C/freqHirArr_Hz, 2.0)

    # Setup the figure
    fig = Figure()
    fig.set_size_inches([8,8])
    
    # Plot psi versus lambda^2
    ax1 = fig.add_subplot(111)
    plot_psi_vs_lamsq_ax(ax=ax1,
                         lamSqArr_m2 = lamSqArr_m2,
                         qArr        = qArr,
                         uArr        = uArr,
                         dqArr       = dqArr,
                         duArr       = duArr,
                         lamSqHirArr_m2 = lamSqHirArr_m2,
                         qModArr     = qModArr,
                         uModArr     = uModArr, 
                         axisYright  = False)
    
    # Write to the pipe
    if io=='string':
        sio = StringIO.StringIO()
        setattr(sio, "name", "foo.jpg")
        fig.savefig(sio, format='jpg' )    
        return sio
    else:
        return fig

#-----------------------------------------------------------------------------#
def plotFracPol(dataMan, indx, io='fig'):

    # Get the data
    freqArr_Hz, IArr_Jy, rmsIArr_Jy= dataMan.get_specI_byindx(indx)
    dummy, modIArr_Jy = dataMan.get_modI_byindx(indx)
    dummy, QArr_Jy, rmsQArr_Jy= dataMan.get_specQ_byindx(indx)
    dummy, UArr_Jy, rmsUArr_Jy= dataMan.get_specU_byindx(indx)
    
    # Calculate fractional polarisation spectra
    qArr = QArr_Jy / modIArr_Jy
    uArr = UArr_Jy / modIArr_Jy
    dqArr = qArr * np.sqrt( (rmsQArr_Jy/QArr_Jy)**2.0 + 
                            (rmsIArr_Jy/IArr_Jy)**2.0 )
    duArr = uArr * np.sqrt( (rmsUArr_Jy/UArr_Jy)**2.0 + 
                            (rmsIArr_Jy/IArr_Jy)**2.0 )
    lamSqArr_m2 = np.power(C/freqArr_Hz, 2.0)

    # Get the models to overplot
    freqHirArr_Hz, qModArr, uModArr = dataMan.get_thin_qumodel_byindx(indx, 
                                                             oversample=True)
    lamSqHirArr_m2 = np.power(C/freqHirArr_Hz, 2.0)

    # Setup the figure
    fig = Figure()
    fig.set_size_inches([8,8])

    # Plot p, q, u versus lambda^2
    ax2 = fig.add_subplot(111)
    plot_pqu_vs_lamsq_ax(ax=ax2,
                         lamSqArr_m2 = lamSqArr_m2,
                         qArr        = qArr,
                         uArr        = uArr,
                         dqArr       = dqArr,
                         duArr       = duArr,
                         lamSqHirArr_m2 = lamSqHirArr_m2,
                         qModArr     = qModArr,
                         uModArr     = uModArr)

    # Write to the pipe
    if io=='string':
        sio = StringIO.StringIO()
        setattr(sio, "name", "foo.jpg")
        fig.savefig(sio, format='jpg' )    
        return sio
    else:
        return fig


#-----------------------------------------------------------------------------#
def plotFracQvsU(dataMan, indx, io='fig'):


    # Get the data
    freqArr_Hz, IArr_Jy, rmsIArr_Jy= dataMan.get_specI_byindx(indx)
    dummy, modIArr_Jy = dataMan.get_modI_byindx(indx)
    dummy, QArr_Jy, rmsQArr_Jy= dataMan.get_specQ_byindx(indx)
    dummy, UArr_Jy, rmsUArr_Jy= dataMan.get_specU_byindx(indx)
    
    # Calculate fractional polarisation spectra
    qArr = QArr_Jy / modIArr_Jy
    uArr = UArr_Jy / modIArr_Jy
    dqArr = qArr * np.sqrt( (rmsQArr_Jy/QArr_Jy)**2.0 + 
                            (rmsIArr_Jy/IArr_Jy)**2.0 )
    duArr = uArr * np.sqrt( (rmsUArr_Jy/UArr_Jy)**2.0 + 
                            (rmsIArr_Jy/IArr_Jy)**2.0 )
    lamSqArr_m2 = np.power(C/freqArr_Hz, 2.0)

    # Get the models to overplot
    freqHirArr_Hz, qModArr, uModArr = dataMan.get_thin_qumodel_byindx(indx, 
                                                             oversample=True)
    lamSqHirArr_m2 = np.power(C/freqHirArr_Hz, 2.0)

    # Setup the figure
    fig = Figure()
    fig.set_size_inches([8,8])

    # Plot U versus Q
    ax1 = fig.add_subplot(111)
    plot_q_vs_u_ax(ax=ax1,
                   lamSqArr_m2 = lamSqArr_m2,
                   qArr        = qArr,
                   uArr        = uArr,
                   dqArr       = dqArr,
                   duArr       = duArr,
                   lamSqHirArr_m2 = lamSqHirArr_m2,
                   qModArr     = qModArr,
                   uModArr     = uModArr, 
                   axisYright  = False)
    
    # Write to the pipe
    if io=='string':
        sio = StringIO.StringIO()
        setattr(sio, "name", "foo.jpg")
        fig.savefig(sio, format='jpg' )    
        return sio
    else:
        return fig

    
#-----------------------------------------------------------------------------#
def plotPolsummary(dataMan, indx, io='fig'):
    
    # Get the data
    freqArr_Hz, IArr_Jy, rmsIArr_Jy= dataMan.get_specI_byindx(indx)
    dummy, modIArr_Jy = dataMan.get_modI_byindx(indx)
    dummy, QArr_Jy, rmsQArr_Jy= dataMan.get_specQ_byindx(indx)
    dummy, UArr_Jy, rmsUArr_Jy= dataMan.get_specU_byindx(indx)
    
    # Calculate fractional polarisation spectra
    qArr = QArr_Jy / modIArr_Jy
    uArr = UArr_Jy / modIArr_Jy
    dqArr = qArr * np.sqrt( (rmsQArr_Jy/QArr_Jy)**2.0 + 
                            (rmsIArr_Jy/IArr_Jy)**2.0 )
    duArr = uArr * np.sqrt( (rmsUArr_Jy/UArr_Jy)**2.0 + 
                            (rmsIArr_Jy/IArr_Jy)**2.0 )
    lamSqArr_m2 = np.power(C/freqArr_Hz, 2.0)

    # Get the models to overplot
    freqHirArr_Hz, qModArr, uModArr = dataMan.get_thin_qumodel_byindx(indx, 
                                                             oversample=True)
    lamSqHirArr_m2 = np.power(C/freqHirArr_Hz, 2.0)

    # Setup the figure
    fig = Figure()
    fig.set_size_inches([8,8])

    # Plot I versus nu,
    ax1 = fig.add_subplot(221)
    plot_I_vs_nu_ax(ax=ax1,
                    freqArr_Hz  = freqArr_Hz,
                    IArr_mJy    = IArr_Jy*1e3,
                    dIArr_mJy   = rmsIArr_Jy*1e3,
                    IModArr_mJy = modIArr_Jy*1e3,
                    axisXtop=True)

    # Plot p, q, u versus lambda^2
    ax2 = fig.add_subplot(223)
    plot_pqu_vs_lamsq_ax(ax=ax2,
                         lamSqArr_m2 = lamSqArr_m2,
                         qArr        = qArr,
                         uArr        = uArr,
                         dqArr       = dqArr,
                         duArr       = duArr,
                         lamSqHirArr_m2 = lamSqHirArr_m2,
                         qModArr     = qModArr,
                         uModArr     = uModArr)
    
    # Plot psi versus lambda^2
    ax3 = fig.add_subplot(222)
    plot_psi_vs_lamsq_ax(ax=ax3,
                         lamSqArr_m2 = lamSqArr_m2,
                         qArr        = qArr,
                         uArr        = uArr,
                         dqArr       = dqArr,
                         duArr       = duArr,
                         lamSqHirArr_m2 = lamSqHirArr_m2,
                         qModArr     = qModArr,
                         uModArr     = uModArr,
                         axisYright=True,
                         axisXtop=True)
    
    # Plot U versus Q
    ax4 = fig.add_subplot(224)
    plot_q_vs_u_ax(ax=ax4,
                   lamSqArr_m2 = lamSqArr_m2,
                   qArr        = qArr,
                   uArr        = uArr,
                   dqArr       = dqArr,
                   duArr       = duArr,
                   lamSqHirArr_m2 = lamSqHirArr_m2,
                   qModArr     = qModArr,
                   uModArr     = uModArr, 
                   axisYright  = True)

    # Write to the pipe
    if io=='string':
        sio = StringIO.StringIO()
        setattr(sio, "name", "foo.jpg")
        fig.savefig(sio, format='jpg' )    
        return sio
    else:
        return fig
    
    
#-----------------------------------------------------------------------------#
def plotRMSF(dataMan, indx, io='fig'):

    # Get the data and Gaussian fit to RMSF
    phiArr, RMSFArr =  dataMan.get_RMSF_byindx(indx)
    pDict = dataMan.get_RMSF_params_byindx(indx)
    
    # Setup the figure
    fig = Figure()
    fig.set_size_inches([8,8])
    
    # Plot the RMSF
    ax1 = fig.add_subplot(111)    
    plot_RMSF_ax(ax=ax1,
                 phiArr  = phiArr,
                 RMSFArr = RMSFArr,
                 fwhmRMSF=pDict["fwhmRMSF"])
    
    # Write to the pipe
    if io=='string':
        sio = StringIO.StringIO()
        setattr(sio, "name", "foo.jpg")
        fig.savefig(sio, format='jpg' )    
        return sio
    else:
        return fig


#-----------------------------------------------------------------------------#
def plotDirtyFDF(dataMan, indx, io='fig'):

    # Get the data
    phiArr, FDFArr_Jy = dataMan.get_dirtyFDF_byindx(indx)

    # Get the peak results
    pDict = dataMan.get_FDF_peak_params_byindx(indx)
    pDict1 = dataMan.get_RMSF_params_byindx(indx)
    gaussParm=[pDict["ampPeakPIfit_Jybm"]*1e3,
               pDict["phiPeakPIfit_rm2"],
               pDict1["fwhmRMSF"]]
    
    # Setup the figure
    fig = Figure()
    fig.set_size_inches([8,8])
    
    # Plot the FDF
    ax1 = fig.add_subplot(111)
    plot_dirtyFDF_ax(ax=ax1,
                     phiArr     = phiArr,
                     FDFArr_mJy = FDFArr_Jy*1e3,
                     gaussParm  = gaussParm,
                     title="Dirty Faraday Dispersion Function")

    # Write to the pipe
    if io=='string':
        sio = StringIO.StringIO()
        setattr(sio, "name", "foo.jpg")
        fig.savefig(sio, format='jpg' )    
        return sio
    else:
        return fig


#-----------------------------------------------------------------------------#
def plotCleanFDF(dataMan, indx, io='fig'):

    # Get the data
    phiArr, dirtyFDFArr_Jy = dataMan.get_dirtyFDF_byindx(indx)
    dummy, cleanFDFArr_Jy = dataMan.get_cleanFDF_byindx(indx)
    dummy, ccFDF_Jy = dataMan.get_ccFDF_byindx(indx)
    
    # Get the peak results
    pDict = dataMan.get_FDF_peak_params_byindx(indx, doClean=True)
    pDict1 = dataMan.get_RMSF_params_byindx(indx)
    gaussParm=[pDict["ampPeakPIfit_Jybm"]*1e3,
               pDict["phiPeakPIfit_rm2"],
               pDict1["fwhmRMSF"]]
    
    # Setup the figure
    fig = Figure()
    fig.set_size_inches([8,8])
    
    # Plot the clean FDF
    ax1 = fig.add_subplot(111)
    plot_cleanFDF_ax(ax=ax1,
                     phiArr          = phiArr,
                     cleanFDFArr_mJy = cleanFDFArr_Jy*1e3,
                     ccFDFArr_mJy    = ccFDF_Jy*1e3,
                     dirtyFDFArr_mJy = dirtyFDFArr_Jy*1e3,
                     gaussParm       = gaussParm,
                     title           = "Clean Faraday Dispersion Function",
                     cutoff_mJy      =  pDict["cleanCutoff_Jybm"]*1e3)

    # Write to the pipe
    if io=='string':
        sio = StringIO.StringIO()
        setattr(sio, "name", "foo.jpg")
        fig.savefig(sio, format='jpg' )    
        return sio
    else:
        return fig
    

#-----------------------------------------------------------------------------#
def plotStampI(dataMan, indx, io='fig'):

    # Get the data & header of the saved postage stamp
    data, head = dataMan.get_stampI_byindx(indx)
    
    # Setup the figure
    fig = Figure()
    fig.set_size_inches([8,8])

    fig = plot_fits_map(data, head, fig=fig)

    # Write to the pipe
    if io=='string':
        sio = StringIO.StringIO()
        setattr(sio, "name", "foo.jpg")
        fig.savefig(sio, format='jpg' )    
        return sio
    else:
        return fig
    

#-----------------------------------------------------------------------------#
def plotStampP(dataMan, indx, io='fig'):

    # Get the data & header of the saved postage stamp
    data, head = dataMan.get_stampP_byindx(indx)

    # Setup the figure
    fig = Figure()
    fig.set_size_inches([8,8])

    fig = plot_fits_map(data, head, fig=fig)

    # Write to the pipe
    if io=='string':
        sio = StringIO.StringIO()
        setattr(sio, "name", "foo.jpg")
        fig.savefig(sio, format='jpg' )    
        return sio
    else:
        return fig
    

#-----------------------------------------------------------------------------#
def plotSctHstQuery(dataMan, plotParm, io='fig'):

    # What type of plot are we creating?
    plotType = plotParm.configDict.get("TYPE", "Histogram")
    
    # Execute each query in turn and store results in list of recarrays
    popLst = []
    names = []
    nCols = 0
    for i in range(len(plotParm.queryLst)-1,-1,-1):
        sql =  plotParm.queryLst[i]
        try:
            resultArr = dataMan.query_database(sql)
            colNames = resultArr.dtype.names
            nCols = max(len(colNames), nCols)
            popLst.append(resultArr)            
        except Exception:
            popLst.append(None)
            print "\nWarn: failed to execute query:"
            print "'%s'\n" % sql
            print traceback.format_exc(), "\n"
    popLst.reverse()
    popLst = popLst[:4]

    # Filter data for limits given in the driving file (default None)
    xMinDataCmd = plotParm.configDict.get("XDATAMIN", None)
    xMaxDataCmd = plotParm.configDict.get("XDATAMAX", None)
    xMinData = xfloat(xMinDataCmd, None)
    xMaxData = xfloat(xMaxDataCmd, None)
    yMinDataCmd = plotParm.configDict.get("YDATAMIN", None)
    yMaxDataCmd = plotParm.configDict.get("YDATAMAX", None)
    yMinData = xfloat(yMinDataCmd, None)
    yMaxData = xfloat(yMaxDataCmd, None)
    zMinDataCmd = plotParm.configDict.get("ZDATAMIN", None)
    zMaxDataCmd = plotParm.configDict.get("ZDATAMAX", None)
    zMinData = xfloat(zMinDataCmd, None)
    zMaxData = xfloat(zMaxDataCmd, None)
    for i in range(len(popLst)):
        msk = filter_range_indx(popLst[i][colNames[0]], xMinData, xMaxData)
        if plotType=="Scatter" and nCols>1:
            msk += filter_range_indx(popLst[i][colNames[1]], yMinData, yMaxData)
        if plotType=="Scatter" and nCols>2:
            msk += filter_range_indx(popLst[i][colNames[2]], zMinData, zMaxData)
        popLst[i] = popLst[i][~msk]
        
    # Labels from driving file (default column name in DB)
    if plotParm.configDict["XLABEL"]=="":
        plotParm.configDict["XLABEL"] = colNames[0]
    xLabel = plotParm.configDict.get("XLABEL", colNames[0])
    if plotParm.configDict["YLABEL"]=="":
        if plotType=="Scatter":
            plotParm.configDict["YLABEL"] = colNames[1]
        else:
            plotParm.configDict["YLABEL"] = "Count"
    yLabel = plotParm.configDict.get("YLABEL", "Count")
    if plotParm.configDict["ZLABEL"]=="":
        if plotType=="Scatter":
            plotParm.configDict["ZLABEL"] = colNames[2]
    zLabel = plotParm.configDict.get("ZLABEL", "")
    plotTitle = plotParm.configDict.get("TITLE", "")

    # Other driving parameters    
    nBins = xint(plotParm.configDict.get("NBINS", 10))
    doXlog = xint(plotParm.configDict.get("DOLOGX", 0))
    doYlog = xint(plotParm.configDict.get("DOLOGY", 0))
    zPower = xfloat(plotParm.configDict.get("ZPOWER", 1.0))

    # Setup the figure
    fig = Figure()
    fig.set_size_inches([8,8])
    ax = fig.add_subplot(111)

    # Bin the data and create the histogram
    if plotType=="Histogram":
        plot_hist4_ax(ax,
                      popLst = popLst,
                      nBins = nBins,
                      doXlog = doXlog,
                      doYlog = doYlog,
                      styIndx = 0,
                      xMin = None,
                      xMax = None,
                      yMax = None,
                      xLabel = xLabel,
                      yLabel = yLabel,
                      title = plotTitle,
                      legLabLst = plotParm.queryLabLst)
    if plotType=="Scatter":    
        plot_scatter4_ax(ax,
                         popLst = popLst,
                         doXlog = doXlog,
                         doYlog = doYlog,
                         zPower = zPower,
                         styIndx = 0,
                         xMin = None,
                         xMax = None,
                         yMin = None,
                         yMax = None,
                         zMin = None,
                         zMax = None,
                         xLabel = xLabel,
                         yLabel = yLabel,
                         zLabel = zLabel,
                         title = plotTitle,
                         legLabLst = plotParm.queryLabLst,
                         showCbar = False,
                         show11Line = False,
                         legLoc = "tr",
                         verbose = False)


        
    # Write to the pipe
    if io=='string':
        sio = StringIO.StringIO()
        setattr(sio, "name", "foo.jpg")
        fig.savefig(sio, format='jpg' )    
        return sio
    else:
        return fig


    
#-----------------------------------------------------------------------------#
def mk_hist_poly(bins, n, logScaleY=False, zeroPt=0.8, addZeroPt=True):
    """Create the line segments for the a polygon used to draw a histogram"""

    if logScaleY is True:
        for i in range(len(n)):
            if n[i] <= 0.0:
                n[i] = zeroPt 
    else:
        zeroPt = 0.0

    # Starting position
    polyCoordLst = []
    if addZeroPt:
        polyCoordLst.append([bins[0],zeroPt])

    # Form the line segments
    i = 0
    j = 0
    while i <= len(bins)-1:
        if j < len(n):
            polyCoordLst.append([bins[i],n[j]])
        if i == j:
            i += 1
        else:
            j += 1

    # Ground the polygon line and close
    if addZeroPt:
        polyCoordLst.append([bins[-1],zeroPt])
        polyCoordLst.append([bins[0],zeroPt])
    polyCoords = np.array(polyCoordLst)
                        
    return polyCoords


#-----------------------------------------------------------------------------#
def label_format_exp(switchExp=3.0):
    """Return a function to format labels for log axes. Switches to a power
    format for log10(|number|) >= switchExp."""
    
    def rfunc(num, pos=None):
        absNum = 0.0
        sign = ""
        exponent = 0.0
        if num!=0.0:
            absNum = abs(num)
            sign = "-" if int(num/absNum)<0 else ""
            exponent = m.log10(absNum)
        if abs(exponent)>=switchExp:
            return r"$%s10^{%i}$" % (sign, m.log10(absNum))
        else:
            return  r"$%s%g$" % (sign, absNum)
    return rfunc
