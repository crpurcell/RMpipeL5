#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     util_plotTk.py                                                    #
#                                                                             #
# PURPOSE:  Plotting functions for the POSSUM pipeline Tk interface.          #
#                                                                             #
# MODIFIED: 19-May-2015 by C. Purcell                                         #
#                                                                             #
# CONTENTS:                                                                   #
#                                                                             #
# tweakAxFormat                                                               #
# plotSpecParms                                                               #
# plot_QUfit_spectra                                                          #
# plotDirtyFDF                                                                #
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
mpl.rcParams['font.size'] = 10.0

# Constants
C = 2.99792458e8


#-----------------------------------------------------------------------------#
def tweakAxFormat(ax, pad=10, loc='upper right', linewidth=1):
    
    # Axis/tic formatting
    ax.tick_params(pad=pad)
    for line in ax.get_xticklines() + ax.get_yticklines():
        line.set_markeredgewidth(linewidth)
        
    # Legend formatting
    leg = ax.legend(numpoints=1, loc=loc, shadow=False,
                     borderaxespad=0.7)
    for t in leg.get_texts():
        t.set_fontsize('small') 
    leg.get_frame().set_linewidth(0.5)

    return ax


#-----------------------------------------------------------------------------#
def plotSpecParms(sessionPath, uniqueName, io='fig'):
    
    specDir = sessionPath + '/OUT'
    
    # Setup the figure
    fig = Figure()
    fig.set_size_inches([8,4])

    # Load the Stokes I spectrum and model
    specIdat = specDir +  '/' + uniqueName +  '_specI.dat'
    IArr = np.loadtxt(specIdat, delimiter=' ', unpack=True)
    
    rmsSpecIDat = specDir +  '/' + uniqueName +  '_rmsSpecI.dat'
    rmsIArr = np.loadtxt(rmsSpecIDat, delimiter=' ', unpack=True)
    
    specImodelDat = specDir +  '/' + uniqueName +  '_specImodel.dat'
    modelIArr = np.loadtxt(specImodelDat, delimiter=' ', unpack=True)

    # Load the Stokes Q and U spectra
    specQdat = specDir +  '/' + uniqueName +  '_specQ.dat'
    QArr = np.loadtxt(specQdat, delimiter=' ', unpack=True)
    
    rmsSpecQdat = specDir +  '/' + uniqueName +  '_rmsSpecQ.dat'
    rmsQArr = np.loadtxt(rmsSpecQdat, delimiter=' ', unpack=True)
    
    specUdat = specDir +  '/' + uniqueName +  '_specU.dat'
    UArr = np.loadtxt(specUdat, delimiter=' ', unpack=True)

    rmsSpecUdat = specDir +  '/' + uniqueName +  '_rmsSpecU.dat'
    rmsUArr = np.loadtxt(rmsSpecUdat, delimiter=' ', unpack=True)
    
    # Calculate fractional polarisation spectra
    qArr = QArr[1] / modelIArr[1]
    uArr = UArr[1] / modelIArr[1]
    dqArr = qArr * np.sqrt( (rmsQArr[1]/QArr[1])**2.0 + 
                            (rmsIArr[1]/IArr[1])**2.0 )
    duArr = uArr * np.sqrt( (rmsUArr[1]/UArr[1])**2.0 + 
                            (rmsIArr[1]/IArr[1])**2.0 )
    lamSqArr_m2 = np.power(C/IArr[0], 2.0)
    plot_QUfit_spectra(fig=fig,
                       lamSqArr_m2 = lamSqArr_m2,
                       IArr_mJy = IArr[1]*1e3,
                       dIArr_mJy = rmsIArr[1]*1e3,
                       qArr = qArr,
                       dqArr = dqArr,
                       uArr = uArr,
                       duArr = duArr,
                       IModArr_mJy =modelIArr[1]*1e3 )

    # Write to the pipe
    if io=='string':
        sio = StringIO.StringIO()
        setattr(sio, "name", "foo.jpg")
        fig.savefig(sio, format='jpg' )    
        return sio
    else:
        return fig

    
#-----------------------------------------------------------------------------#
def plot_QUfit_spectra(fig, lamSqArr_m2, IArr_mJy, dIArr_mJy, qArr, dqArr,
                       uArr, duArr, lamSqHirArr_m2=None, IModArr_mJy=None,
                       qModArr=None, uModArr=None, outName=None):
    
    if lamSqHirArr_m2 is None:
        lamSqHirArr_m2 = lamSqArr_m2
    freqArr_Hz =  C / lamSqArr_m2**0.5
    freqHirArr_Hz =  C / lamSqHirArr_m2**0.5

    # Plot I versus nu,
    ax1 = fig.add_subplot(221)
    ax1.errorbar(x=freqArr_Hz/1e9, y=IArr_mJy, yerr=dIArr_mJy, mfc='none',
                 ms=4, fmt='D', ecolor='grey', elinewidth=1.0, capsize=2)
    if IModArr_mJy is not None:
        ax1.plot(freqHirArr_Hz/1e9, IModArr_mJy, color='k', lw=0.5,
                 label='I Model')
    ax1.text(0.95, 0.90, 'Stokes I Spectrum', transform=ax1.transAxes,
             horizontalalignment='right')
    ax1.yaxis.set_major_locator(MaxNLocator(4))
    ax1.xaxis.set_major_locator(MaxNLocator(5))
    ax1.set_xlabel('$\\nu$ (GHz)')
    ax1.set_ylabel('Flux Density (mJy)')
    format_ticks(ax1, 10, 1.2)
    ax1.relim()
    ax1.autoscale_view(False,True,True)

    # Plot p, q, u versus lambda^2
    ax1 = fig.add_subplot(223)
    ax1.errorbar(x=lamSqArr_m2, y=qArr, yerr=dqArr, mec='b', mfc='none', ms=4,
                 fmt='D', ecolor='b', elinewidth=1.0, capsize=2,
                 label='Stokes q')
    ax1.errorbar(x=lamSqArr_m2, y=uArr, yerr=duArr, mec='r', mfc='none', ms=4,
                 fmt='D', ecolor='r', elinewidth=1.0, capsize=2,
                 label='Stokes u')
    pArr = np.sqrt(qArr**2.0 + uArr**2.0 )
    dpArr = np.sqrt(dqArr**2.0 + duArr**2.0 )
    ax1.errorbar(x=lamSqArr_m2, y=pArr, yerr=dpArr, mec='k', mfc='none', ms=4,
                 fmt='D', ecolor='k', elinewidth=1.0, capsize=2, label='p')
    if qModArr is not None:
        ax1.plot(lamSqHirArr_m2, qModArr, color='b', lw=0.5, label='Model q')
    if uModArr is not None:
        ax1.plot(lamSqHirArr_m2, uModArr, color='r', lw=0.5, label='Model u')
    if qModArr is not None and uModArr is not None:
        PModArr = np.sqrt(qModArr**2.0 + uModArr**2.0 )
        ax1.plot(lamSqHirArr_m2, PModArr, color='k', lw=0.5, label='Model p')
    ax1.yaxis.set_major_locator(MaxNLocator(4))
    ax1.xaxis.set_major_locator(MaxNLocator(4))
    ax1.set_xlabel('$\\lambda^2$ (m$^2$)')
    ax1.set_ylabel('Fractional Polarisation')
    format_ticks(ax1, 10, 1.2)
    ax1.relim()
    ax1.autoscale_view(False,True,True)
      
    # Plot psi versus lambda^2
    ax1 = fig.add_subplot(222)
    ax1.yaxis.tick_right()
    ax1.yaxis.set_label_position("right")
    dQUArr = np.sqrt(dqArr**2.0 + duArr**2.0)
    psiArr_deg = np.degrees( np.arctan2(uArr, qArr) / 2.0 )
    dPsiArr_deg = np.degrees( np.sqrt( (qArr * duArr)**2.0 +
                                       (uArr * dqArr)**2.0) /
                              (2.0 * pArr**2.0) )
    ax1.errorbar(x=lamSqArr_m2, y=psiArr_deg, yerr=dPsiArr_deg, mec='k',
                 mfc='none', ms=4, fmt='D', ecolor='k', elinewidth=1.0,
                 capsize=2, label='$\psi$')
    if qModArr is not None and uModArr is not None:
        psiHirArr_deg = np.degrees( np.arctan2(uModArr, qModArr) / 2.0 )
        ax1.plot(lamSqHirArr_m2, psiHirArr_deg, color='k', lw=0.5,
                 label='Model $\psi$')
    ax1.set_ylim(-99.9, 99.9)
    ax1.yaxis.set_major_locator(MaxNLocator(4))
    ax1.xaxis.set_major_locator(MaxNLocator(4))
    ax1.set_xlabel('$\\lambda^2$ (m$^2$)')
    ax1.set_ylabel('$\psi$ (degrees)')
    format_ticks(ax1, 10, 1.2)
    ax1.relim()
    ax1.autoscale_view(False,True,True)

    # Plot U versus Q
    ax1 = fig.add_subplot(224)
    ax1.yaxis.tick_right()
    ax1.yaxis.set_label_position("right")
    ax1.errorbar(x=uArr, y=qArr, xerr=duArr, yerr=dqArr, mec='k', mfc='none',
                 ms=4, fmt='D', ecolor='k', elinewidth=1.0, capsize=2)
    if qModArr is not None and uModArr is not None:
        ax1.plot(uModArr, qModArr, color='k', lw=0.5, label='Model q & u')
    ax1.axhline(0, color='grey')
    ax1.axvline(0, color='grey')
    ax1.yaxis.set_major_locator(MaxNLocator(4))
    ax1.xaxis.set_major_locator(MaxNLocator(4))
    ax1.set_xlabel('Stokes u')
    ax1.set_ylabel('Stokes q')
    format_ticks(ax1, 10, 1.2)
    ax1.relim()
    ax1.autoscale_view(False,True,True)

    return fig
    
    
#-----------------------------------------------------------------------------#
def plotDirtyFDF(sessionPath, uniqueName, io='fig'):

    specDir = sessionPath + '/OUT'
    
    # Setup the figure
    fig = Figure()
    fig.set_size_inches([8,4])
    
    # Load the RMSF and the dirty FDF
    FDFdat = specDir +  '/' + uniqueName +  '_dirtyFDF.dat'
    phiArr, FDFreal, FDFimag = np.loadtxt(FDFdat, unpack=True)
    FDFpi = np.sqrt( np.power(FDFreal, 2.0) + np.power(FDFimag, 2.0) )
    RMSFdat =  specDir +  '/' + uniqueName +  '_RMSF.dat'
    RMSF = np.loadtxt(RMSFdat, delimiter=' ', unpack=True)

    # Plot the RMSF on top
    ax1 = fig.add_axes([0.12, 0.56, 0.85, 0.40],xticklabels=[])
    ax1.step(RMSF[0], RMSF[1], where='mid', color='b', lw='0.5',
             label='Real')
    ax1.step(RMSF[0], RMSF[2], where='mid', color='r', lw='0.5',
             label='Imaginary')
    ax1.step(RMSF[0], np.sqrt(RMSF[1]**2.0 + RMSF[2]**2.0) , where='mid',
            color='k', lw='1.0', label='PI')
    ax1.text(0.05, 0.84, 'RMSF', transform=ax1.transAxes)
    
    # Format tweaks
    ax1 = tweakAxFormat(ax1)
    
    # Plot the FDF on bottom
    ax2 = fig.add_axes([0.12, 0.12, 0.85, 0.40])
    ax2.step(phiArr, FDFreal*1e3, where='mid', color='b', lw='0.5',
             label='Real')
    ax2.step(phiArr, FDFimag*1e3, where='mid', color='r', lw='0.5',
             label='Imaginary')
    ax2.step(phiArr, FDFpi*1e3, where='mid', color='k', lw='1.0',
             label='PI')
    ax2.text(0.05, 0.84, 'Dirty FDF', transform=ax2.transAxes)

    # Format tweaks
    ax2 = tweakAxFormat(ax2)
    
    # Set axis limits and labels
    ax1.set_xlim(np.min(phiArr), np.max(phiArr))
    ax2.set_xlim(np.min(phiArr), np.max(phiArr))
    ax1.set_ylabel('Normalised Units')
    ax2.set_xlabel('$\phi$ rad m$^{-2}$')
    ax2.set_ylabel('Flux Density (mJy bm$^{-1}$)')

    # Write to the pipe
    if io=='string':
        sio = StringIO.StringIO()
        setattr(sio, "name", "foo.jpg")
        fig.savefig(sio, format='jpg' )    
        return sio
    else:
        return fig
    

#-----------------------------------------------------------------------------#
def format_ticks(ax, pad=10, w=1.0):

    ax.tick_params(pad=pad)
    for line in ax.get_xticklines() + ax.get_yticklines():
        line.set_markeredgewidth(w)
