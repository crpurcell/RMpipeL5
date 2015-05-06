#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     util_plot.py                                                      #
#                                                                             #
# PURPOSE:  Plotting functions for the POSSUM pipeline interface.             #
#                                                                             #
# MODIFIED: 02-Dec-2014 by C. Purcell                                         #
#                                                                             #
# CONTENTS:                                                                   #
#                                                                             #
# tweakAxFormat                                                               #
# plotSpecI                                                                   #
# plotSpecQU                                                                  #
# plotSpecRMS                                                                 #
# plotDirtyFDF                                                                #
# plotCleanFDF                                                                #
#                                                                             #
#                                                                             #
#                                                                             #
#=============================================================================#

import os
import sys
import numpy as np
import StringIO
import astropy.io.fits as pf

from util_plotFITS import *
from util_PPC import read_paths

# chose a non-GUI backend and set a temporary directory
import matplotlib as mpl
#mpl.use( 'Agg' )
import pylab as pl
from matplotlib.ticker import MaxNLocator

# If this file does not exist, default is current directory
inPathFile = 'env_paths.txt'
    
# Read the paths to the directory containing DATA and SESSIONS
pathDict = read_paths(inPathFile)
pipeDir = pathDict['pipeDir']

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
def tweakAxFormat(ax1, pad=10, loc='upper right'):
    
    # Axis/tic formatting
    ax1.tick_params(pad=pad)
    for line in ax1.get_xticklines() + ax1.get_yticklines():
        line.set_markeredgewidth(1)
        
    # Legend formatting
    leg = ax1.legend(numpoints=1, loc=loc, shadow=False,
                     borderaxespad=0.7)
    for t in leg.get_texts():
        t.set_fontsize('small') 
    leg.get_frame().set_linewidth(0.5)

    return ax1


#-----------------------------------------------------------------------------#
def plotSpecIerrs(session, uniqueName, io='string'):
    
    # Setup the figure
    fig = pl.figure()
    fig.set_size_inches([8,4])

    sessionRootDir = '%s/SESSION' % pipeDir
    specDir = sessionRootDir + '/' + session + '/OUT'
    
    # Set the paths to the data
    #specDir = sessionPath + '/OUT'
    
    # Load the Stokes I spectrum and model
    specIdat = specDir +  '/' + uniqueName +  '_specI.dat'
    specIArr = np.loadtxt(specIdat, delimiter=' ', unpack=True)
    specdIdat = specDir +  '/' + uniqueName +  '_rmsSpecI.dat'
    specdIArr = np.loadtxt(specIdat, delimiter=' ', unpack=True)
    modelIdat = specDir +  '/' + uniqueName +  '_specImodel.dat'
    modelIArr = np.loadtxt(modelIdat, delimiter=' ', unpack=True)

    # Plot the spectrum and model, in GHz and mJy
    ax = fig.add_axes([0.12, 0.14, 0.85, 0.80])
    plot_I_errs_ax(ax,
                   freqArr_Hz=specIArr[0],
                   IArr_mJy=specIArr[1]*1e3,
                   dIArr_mJy=specdIArr[1],
                   freqHirArr_Hz=modelIArr[0],
                   IModArr_mJy=modelIArr[1]*1e3)

    # Write to the pipe
    sio = StringIO.StringIO()
    setattr(sio, "name", "foo.jpg")
    fig.savefig(sio, format='jpg' )
    
    return sio


#-----------------------------------------------------------------------------#
def plot_I_errs_ax(ax, freqArr_Hz, IArr_mJy, dIArr_mJy, freqHirArr_Hz=None,
                   IModArr_mJy=None):
    
    """
    Plot the I spectrum and an optional model.
    """

    # Default to non-high-resolution inputs
    if freqHirArr_Hz is None:
        freqHirArr_Hz =  C / lamSqHirArr_m2**0.5

    # Plot I versus frequency
    ax.errorbar(x=freqArr_Hz/1e9, y=IArr_mJy, yerr=dIArr_mJy, mfc='none',
                ms=4, fmt='D', ecolor='grey', elinewidth=1.0, capsize=2,
                label='Data')
    if IModArr_mJy is not None:
        ax.plot(freqHirArr_Hz/1e9, IModArr_mJy, color='k', lw=0.5,
                label='I Model')
    ax.text(0.05, 0.89, 'Stokes I Spectrum', transform=ax.transAxes)
    ax.yaxis.set_major_locator(MaxNLocator(4))
    ax.xaxis.set_major_locator(MaxNLocator(5))

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
def plotSpecPQUerrs(session, uniqueName, io='string'):
    
    # Setup the figure
    fig = pl.figure()
    fig.set_size_inches([8,4])

    sessionRootDir = '%s/SESSION' % pipeDir
    specDir = sessionRootDir + '/' + session + '/OUT'
    
    # Set the paths to the data
    #specDir = sessionPath + '/OUT'
    
    # Load the Stokes Q and U spectra
    specQdat = specDir +  '/' + uniqueName +  '_specQ.dat'
    specQArr = np.loadtxt(specQdat, delimiter=' ', unpack=True)
    specdQdat = specDir +  '/' + uniqueName +  '_rmsSpecQ.dat'
    specdQArr = np.loadtxt(specQdat, delimiter=' ', unpack=True)
    specUdat = specDir +  '/' + uniqueName +  '_specU.dat'
    specUArr = np.loadtxt(specUdat, delimiter=' ', unpack=True)
    specdUdat = specDir +  '/' + uniqueName +  '_rmsSpecU.dat'
    specdUArr = np.loadtxt(specQdat, delimiter=' ', unpack=True)
    specPArr = np.sqrt(np.power(specQArr[1],2) + np.power(specUArr[1],2))
    specdPArr = np.sqrt(np.power(specdQArr[1],2) + np.power(specdUArr[1],2))
    
    # Plot the spectrum and model, in GHz and mJy
    ax = fig.add_axes([0.12, 0.14, 0.85, 0.80])
    plot_PQU_errs_ax(ax,
                     freqArr_Hz=specQArr[0],
                     PArr_mJy=specPArr*1e3,
                     dPArr_mJy=specdPArr*1e3,
                     QArr_mJy=specQArr[1]*1e3,
                     dQArr_mJy=specdQArr[1]*1e3,
                     UArr_mJy=specUArr[1]*1e3,
                     dUArr_mJy=specdUArr[1]*1e3)
                     

    # Write to the pipe
    sio = StringIO.StringIO()
    setattr(sio, "name", "foo.jpg")
    fig.savefig(sio, format='jpg' )
    
    return sio


#-----------------------------------------------------------------------------#
def plot_PQU_errs_ax(ax, freqArr_Hz, PArr_mJy, dPArr_mJy, QArr_mJy, dQArr_mJy,
                   UArr_mJy, dUArr_mJy):
    
    """
    Plot the P, Q & U spectrum and an optional model.
    """
    
    # Plot P, Q, U versus frequency
    ax.errorbar(x=freqArr_Hz/1e9, y=QArr_mJy, yerr=dQArr_mJy, mec='b',
                mfc='none', ms=4, fmt='D', color='b', elinewidth=0.3,
                capsize=2, label='Stokes Q', linestyle='-')
    ax.errorbar(x=freqArr_Hz/1e9, y=UArr_mJy, yerr=UArr_mJy, mec='r',
                mfc='none', ms=4, fmt='D', color='r', elinewidth=0.3,
                capsize=2, label='Stokes U', linestyle='-')
    ax.errorbar(x=freqArr_Hz/1e9, y=PArr_mJy, yerr=dPArr_mJy, mec='k',
                mfc='none', ms=4, fmt='D', color='k', elinewidth=0.3,
                capsize=2, label='Stokes P', linestyle='-')
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
def plotSpecRMS(session, name, io='string'):

    # Setup the figure
    fig = pl.figure()
    fig.set_size_inches([8,4])

    # Set the paths to the data
    sessionRootDir = '%s/SESSION' % pipeDir
    specDir = sessionRootDir + '/' + session + '/OUT'
    
    # Load the Stokes I, Q and U rms spectra
    specIdat = specDir +  '/' + name +  '_rmsSpecI.dat'
    specIArr = np.loadtxt(specIdat, delimiter=' ', unpack=True)
    specQdat = specDir +  '/' + name +  '_rmsSpecQ.dat'
    specQArr = np.loadtxt(specQdat, delimiter=' ', unpack=True)
    specUdat = specDir +  '/' + name +  '_rmsSpecU.dat'
    specUArr = np.loadtxt(specUdat, delimiter=' ', unpack=True)
    specPIArrY = np.sqrt(np.power(specQArr[1],2) + np.power(specUArr[1],2))

    # Plot the spectra in GHz and mJy
    ax1 = fig.add_axes([0.12, 0.14, 0.85, 0.80])
    ax1.plot(specIArr[0]/1e9, specIArr[1]*1e3, marker='+', color='k',
             lw='0.5', label='rms I')
    ax1.plot(specQArr[0]/1e9, specQArr[1]*1e3, marker='+', color='b',
             lw='0.5', label='rms Q')
    ax1.plot(specUArr[0]/1e9, specUArr[1]*1e3, marker='+', color='r',
             lw='0.5', label='rms U')
    ax1.text(0.05, 0.89, 'I, Q & U RMS', transform=ax1.transAxes)
    
    # Format tweaks
    ax1 = tweakAxFormat(ax1)
    
    # Set axis limits and labels
    ax1.set_xlim( np.min(specQArr[0])/1e9, np.max(specQArr[0])/1e9)
    ax1.set_xlabel('$\\nu$ (GHz)')
    ax1.set_ylabel('Flux Density (mJy bm$^{-1}$)')

    # Write to the pipe
    sio = StringIO.StringIO()
    setattr(sio, "name", "foo.jpg")
    fig.savefig(sio, format='jpg' )
    
    return sio

#-----------------------------------------------------------------------------#
def plotDirtyFDF(session, name, io='string'):

    # Setup the figure
    fig = pl.figure()
    fig.set_size_inches([8,6])

    # Set the paths to the data
    sessionRootDir = '%s/SESSION' % pipeDir
    specDir = sessionRootDir + '/' + session + '/OUT'
    
    # Load the RMSF and the dirty FDF
    FDFdat = specDir +  '/' + name +  '_dirtyFDF.dat'
    phiArr, FDFreal, FDFimag = np.loadtxt(FDFdat, unpack=True)
    FDFpi = np.sqrt( np.power(FDFreal, 2.0) + np.power(FDFimag, 2.0) )
    RMSFdat =  specDir +  '/' + name +  '_RMSF.dat'
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
    sio = StringIO.StringIO()
    setattr(sio, "name", "foo.jpg")
    fig.savefig(sio, format='jpg' )
    
    return sio


#-----------------------------------------------------------------------------#
def plotCleanFDF(session, name, io='string'):

    # Setup the figure
    fig = pl.figure()
    fig.set_size_inches([8,6])

    # Set the paths to the data
    sessionRootDir = '%s/SESSION' % pipeDir
    specDir = sessionRootDir + '/' + session + '/OUT'
    
    # Load the RMSF and the clean FDF
    FDFdat = specDir +  '/' + name +  '_dirtyFDF.dat'
    dummmy, dirtyFDFreal, dirtyFDFimag = np.loadtxt(FDFdat, unpack=True)
    dirtyFDFpi = np.sqrt( np.power(dirtyFDFreal, 2.0) +
                          np.power(dirtyFDFimag, 2.0) )
    FDFdat = specDir +  '/' + name +  '_cleanFDF.dat'
    phiArr, cleanFDFreal, cleanFDFimag = np.loadtxt(FDFdat, unpack=True)
    cleanFDFpi = np.sqrt( np.power(cleanFDFreal, 2.0) +
                     np.power(cleanFDFimag, 2.0) )
    FDFccDat = specDir +  '/' + name +  '_cleanFDFPImodel.dat'
    FDFcc = np.loadtxt(FDFccDat,unpack=True)    
    RMSFdat =  specDir +  '/' + name +  '_RMSF.dat'
    RMSFphiArr, RMSFreal, RMSFimag = np.loadtxt(RMSFdat, unpack=True)
    RMSFpi = np.sqrt( np.power(RMSFreal, 2.0) +
                     np.power(RMSFimag, 2.0) )
    
    # Plot the clean FDF on top
    ax1 = fig.add_axes([0.12, 0.56, 0.85, 0.40], xticklabels=[])
    ax1.step(phiArr, dirtyFDFpi*1e3, where='mid', color='grey',
             lw='0.5', label='Dirty FDF')
    ax1.step(phiArr, cleanFDFreal*1e3, where='mid', color='b', lw='0.5',
             label='Real')
    ax1.step(phiArr, cleanFDFimag*1e3, where='mid', color='r', lw='0.5',
             label='Imaginary')
    ax1.step(phiArr, cleanFDFpi*1e3, where='mid', color='k', lw='1.0',
             label='PI')
    ax1.text(0.05, 0.84, 'Clean FDF', transform=ax1.transAxes)

    # Format tweaks
    ax1 = tweakAxFormat(ax1)
    
    # Plot the clean component spectrum on bottom
    ax2 = fig.add_axes([0.12, 0.12, 0.85, 0.40])
    ax2.step(phiArr, FDFcc[1]*1e3, where='mid', color='g', lw='0.5',
             label='CC')
    ax2.text(0.05, 0.84, 'CC Spectrum', transform=ax2.transAxes)

    # Format tweaks
    ax2 = tweakAxFormat(ax2)
    
    # Set axis limits and labels
    ax1.set_xlim(np.min(phiArr), np.max(phiArr))
    ax2.set_xlim(np.min(phiArr), np.max(phiArr))
    ax2.set_ylim(np.max(FDFcc[1]*1e3)*(-0.1), np.max(FDFcc[1]*1e3)*1.1)
    ax2.set_xlabel('$\phi$ rad m$^{-2}$')
    ax2.set_ylabel('Flux Density (mJy bm$^{-1}$)', position=(1,1.02))

    # Write to the pipe
    sio = StringIO.StringIO()
    setattr(sio, "name", "foo.jpg")
    fig.savefig(sio, format='jpg' )
    
    return sio


#-----------------------------------------------------------------------------#
def plotPostageI(session, name, io='string'):

    # Set the paths to the data
    sessionRootDir = '%s/SESSION' % pipeDir
    specDir = sessionRootDir + '/' + session + '/OUT'

    # Load the FITS data
    fitsI = specDir +  '/' + name +  '_I.fits'
    data = pf.getdata(fitsI)
    head = pf.getheader(fitsI)
    
    # Create the figure
    fig = pl.figure(figsize=(6.5, 5.0))
    fig = plot_fits_map(data*1000, head, fig=fig, bunit='mJy/beam')
    
    # Write to the pipe
    sio = StringIO.StringIO()
    setattr(sio, "name", "foo.jpg")
    fig.savefig(sio, format='jpg' , dpi=70)
    
    return sio


#-----------------------------------------------------------------------------#
def plotPostageP(session, name, io='string'):

    # Set the paths to the data
    sessionRootDir = '%s/SESSION' % pipeDir
    specDir = sessionRootDir + '/' + session + '/OUT'

    # Load the FITS data
    fitsP = specDir +  '/' + name +  '_P.fits'
    data = pf.getdata(fitsP)
    head = pf.getheader(fitsP)
    
    # Create the figure
    fig = pl.figure(figsize=(6.5, 5.0))
    fig = plot_fits_map(data*1000, head, fig=fig, bunit='mJy/beam')
    
    # Write to the pipe
    sio = StringIO.StringIO()
    setattr(sio, "name", "foo.jpg")
    fig.savefig(sio, format='jpg' , dpi=70)
    
    return sio
