#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     util_plotTk.py                                                    #
#                                                                             #
# PURPOSE:  Plotting functions for the POSSUM pipeline Tk interface.          #
#                                                                             #
# MODIFIED: 27-May-2015 by C. Purcell                                         #
#                                                                             #
# CONTENTS:                                                                   #
#                                                                             #
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
#                                                                             #
#=============================================================================#

import os
import sys
import numpy as np
import StringIO
import astropy.io.fits as pf
import matplotlib as mpl
from matplotlib.ticker import MaxNLocator
from matplotlib.figure import Figure

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
                mfc='none', ms=4, fmt='D', color='b', elinewidth=0.3,
                capsize=2, label='Stokes Q', linestyle='-')
    ax.errorbar(x=freqArr_Hz/1e9, y=UArr_mJy, yerr=dUArr_mJy, mec='r',
                mfc='none', ms=4, fmt='D', color='r', elinewidth=0.3,
                capsize=2, label='Stokes U', linestyle='-')
    ax.errorbar(x=freqArr_Hz/1e9, y=PArr_mJy, yerr=dPArr_mJy, mec='k',
                mfc='none', ms=4, fmt='D', color='k', elinewidth=0.3,
                capsize=2, label='Stokes P', linestyle='-')
    
    # Plot the models
    if QmodArr is not None:
        ax.plot(freqArr_Hz/1e9, QmodArr, color='b', lw=0.5, label='Model Q')
    if UmodArr is not None:
        ax.plot(freqArr_Hz/1e9, UmodArr, color='r', lw=0.5, label='Model U')
    if QmodArr is not None and UmodArr is not None:
        PmodArr = np.sqrt(QmodArr**2.0 + UmodArr**2.0 )
        ax.plot(freqArr_Hz/1e9, PmodArr, color='k', lw=0.5, label='Model P')

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
        yGauss = gauss([1.0, 0.0, fwhmRMSF])(phi2Arr)
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
        yGauss = gauss(gaussParm)(phi2Arr)
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
        yGauss = gauss(gaussParm)(phi2Arr)
        ax.plot(phiArr, yGauss, color='magenta',marker='None',mfc='w',
                mec='g', ms=10, label='Peak Fit', lw=2.5, ls='-')

    # Plot the clean cutoff line
    if not cutoff_mJy is None:
        ax.axhline(cutoff_mJy, color="r")

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

    
# Axis code above
#=============================================================================#
# Figure code below


#-----------------------------------------------------------------------------#
def plotSpecIPQU(dataMan, indx, io='fig'):

    # Get the data
    freqArr_Hz, IArr_Jy, rmsIArr_Jy= dataMan.get_specI_byindx(indx)
    dummy, modIArr_Jy = dataMan.get_modI_byindx(indx)
    dummy, QArr_Jy, rmsQArr_Jy= dataMan.get_specQ_byindx(indx)
    dummy, UArr_Jy, rmsUArr_Jy= dataMan.get_specU_byindx(indx)
    
    # Setup the figure
    fig = Figure()
    fig.set_size_inches([8,8])
    
    # Plot I versus nu,
    ax1 = fig.add_subplot(211)
    plot_I_vs_nu_ax(ax=ax1,
                    freqArr_Hz  = freqArr_Hz,
                    IArr_mJy    = IArr_Jy*1e3,
                    dIArr_mJy   = rmsIArr_Jy*1e3,
                    IModArr_mJy = modIArr_Jy*1e3)
    ax1.set_xlabel('')
    [label.set_visible(False) for label in ax1.get_xticklabels()]

    # Plot Stokes P, Q & U
    ax2 = fig.add_subplot(212, sharex=ax1)
    plot_PQU_vs_nu_ax(ax=ax2,
                      freqArr_Hz = freqArr_Hz,
                      QArr_mJy   = QArr_Jy*1e3,
                      UArr_mJy   = UArr_Jy*1e3,
                      dQArr_mJy  = rmsQArr_Jy*1e3,
                      dUArr_mJy  = rmsUArr_Jy*1e3)

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
                         duArr       = duArr) 
        
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
                         duArr       = duArr)
    
    # Plot psi versus lambda^2
    ax3 = fig.add_subplot(222)
    plot_psi_vs_lamsq_ax(ax=ax3,
                         lamSqArr_m2 = lamSqArr_m2,
                         qArr        = qArr,
                         uArr        = uArr,
                         dqArr       = dqArr,
                         duArr       = duArr,
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

    # Get the data
    phiArr, RMSFArr =  dataMan.get_RMSF_byindx(indx)

    # Setup the figure
    fig = Figure()
    fig.set_size_inches([8,8])
    
    # Plot the RMSF
    ax1 = fig.add_subplot(111)    
    plot_RMSF_ax(ax=ax1,
                 phiArr  = phiArr,
                 RMSFArr = RMSFArr)

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
    
    # Setup the figure
    fig = Figure()
    fig.set_size_inches([8,8])
    
    # Plot the FDF
    ax1 = fig.add_subplot(111)
    plot_dirtyFDF_ax(ax=ax1,
                     phiArr     = phiArr,
                     FDFArr_mJy = FDFArr_Jy*1e3,
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
                     gaussParm       = [],
                     title           = "Clean Faraday Dispersion Function")

    # Write to the pipe
    if io=='string':
        sio = StringIO.StringIO()
        setattr(sio, "name", "foo.jpg")
        fig.savefig(sio, format='jpg' )    
        return sio
    else:
        return fig
    
