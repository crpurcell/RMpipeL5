#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     util_RM.py                                                        #
#                                                                             #
# PURPOSE:  Common procedures used with RM-synthesis scripts.                 #
#                                                                             #
# REQUIRED: Requires the numpy and scipy modules.                             #
#                                                                             #
# MODIFIED: 19-November-2015 by C.Purcell.                                    #
#                                                                             #
# CONTENTS:                                                                   #
#                                                                             #
#  do_rmsynth          ... perform RM-synthesis on Q & U data                 #
#  do_rmclean          ... perform Hogbom RM-clean on a dirty FDF             #
#  fits_make_lin_axis  ... create an array of absica values for a lin axis    #
#  get_RMSF            ... return the RMSF given a weight and sampling        #
#  extrap              ... interpolate and extrapolate an array               #
#  fit_rmsf            ... fit a Gaussian to the main lobe of the RMSF        #
#  detect_peak         ... detect the extent of a peak in a 1D array          #
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


import sys, os
import math as m
import numpy as np
from mpfit import mpfit

# Constants
C = 2.99792458e8


#-----------------------------------------------------------------------------#
def do_rmsynth(dataQ, dataU, lamSqArr, phiArr, weight=None, dType='float32'):
    """
    Perform RM-synthesis on Stokes Q and U cubes
    """

    # Parse the weight argument
    if weight is None:
        wtArr = np.ones(lamSqArr.shape, dtype=dType)
    else:
        wtArr = np.array(weight, dtype=dType)
        if not wtArr.shape  == lamSqArr.shape:
            print "  Err: LambdaSq and Weight arrays must be the same shape."
            sys.exit(1)

    # Sanity check on Q & U data array sizes
    if not dataQ.shape == dataU.shape:
        print "  Err: Stokes Q and U data arrays be the same shape."
        sys.exit(1)
    nDims = len(dataQ.shape)
    if not nDims <= 3:
        print "  Err: data-dimensions must be <= 3."
        sys.exit(1)

    # Check that the first data axis is the same size as the lamda array
    if not dataQ.shape[-1] == lamSqArr.shape[-1]:
        print "  Err: The Stokes Q and U arrays mush be in spectral order."
        print "       # Stokes = %d, # Lamda = %d." % (dataQ.shape[-1],
                                                       lamSqArr.shape[-1])
        sys.exit(1)

    # Create a blanking mask assuming NaNs are blanked in the input data
    # Set the weight = 0 in fully blanked planes
    dataMsk = np.where(np.isnan(dataQ) + np.isnan(dataU), 0, 1)
    dataMsk = np.where(np.sum(np.sum(dataMsk, 0), 0)==0, 1, 0)
    wtArr = np.where(dataMsk==1, 0.0, wtArr)
    del dataMsk

    # Reshape data arrays to allow the same recipies to work on all
    if nDims==1:
        dataQ = np.reshape(dataQ, (1, 1, dataQ.shape[-1]))
        dataU = np.reshape(dataU, (1, 1, dataU.shape[-1]))
    elif nDims==2:
        dataQ = np.reshape(dataQ, (1, dataQ.shape[-2], dataQ.shape[-1]))
        dataU = np.reshape(dataU, (1, dataU.shape[-2], dataU.shape[-1]))
        
    # Calculate the dimensions of the output RM cube
    nY = dataQ.shape[0]
    nX = dataQ.shape[1]
    nPhi = phiArr.shape[0]

    # Initialise the complex Faraday Dispersion Function (FDF) cube
    # Remember, python index order is reversed [2,1,0] = [y,x,phi]
    FDFcube = np.ndarray((nY, nX, nPhi), dtype='complex')
    
    # B&dB equations (24) and (38) give the inverse sum of the weights
    K = 1.0 / np.nansum(wtArr)

    # Get the weighted mean of the LambdaSq distribution (B&dB Eqn. 32)
    lam0Sq = K * np.nansum(wtArr * lamSqArr)
    
    # Mininize the number of inner-loop operations by calculating the
    # argument of the EXP term in B&dB Eqns. (25) and (36) for the FDF
    # Returned array has dimensions (nPhi x nFreq)
    a = (-2.0 * 1j * phiArr)
    b = (lamSqArr - lam0Sq) 
    arg = np.exp( np.outer(a, b) )
        
    # Create a weighted complex polarised surface-brightness cube
    # i.e., observed polarised surface brightness, B&dB Eqns. (8) and (14)
    # Weight-array will broadcast to the spectral dimension
    # Cube has dimensions (nX, nY, nFreq)
    PobsCube = (dataQ + 1j * dataU) * wtArr

    # Do the synthesis at each pixel of the image
    for k in range(nY):
        for i in range(nX):
                
            # Calculate the Faraday Dispersion Function
            # B&dB Eqns. (25) and (36)
            FDFcube[k,i,:] = K * np.nansum(PobsCube[k,i,:] * arg, 1)
            
    # Calculate the complex Rotation Measure Spread Function
    RMSFArr, phiSamp, fwhmRMSF = get_RMSF(lamSqArr, phiArr, wtArr, lam0Sq)

    # Reshape the data and FDF cube back to their original shapes
    if nDims==1:
        dataQ = np.reshape(dataQ, (dataQ.shape[-1]))
        dataU = np.reshape(dataU, (dataU.shape[-1]))
        FDFcube = np.reshape(FDFcube, (FDFcube.shape[-1]))
    elif nDims==2:
        dataQ = np.reshape(dataQ, (dataQ.shape[-2], dataQ.shape[-1]))
        dataU = np.reshape(dataU, (dataU.shape[-2], dataU.shape[-1]))
        FDFcube = np.reshape(FDFcube, (FDFcube.shape[-2], FDFcube.shape[-1]))
    
    return FDFcube, [phiSamp, RMSFArr], lam0Sq, fwhmRMSF


#-----------------------------------------------------------------------------#
def do_rmclean(dirtyFDF, phiArr, lamSqArr, cutoff, maxIter=1000, gain=0.1,
               weight=None, RMSFArr=None, RMSFphiArr=None, fwhmRMSF=None,
               doDerotate=False, dType='float32', doPlots=True):
    """
    Perform Hogbom (Heald) clean on a single RM spectrum.
    """

    # Initial sanity checks --------------------------------------------------#
   
    # Check that dirtyFDF is ID and get its length
    if len(dirtyFDF.shape) != 1:
        print "Err: the dirty FDF is not a 1D array."
        sys.exit(1)
    nFDF = dirtyFDF.shape[0]

    # Check that dirtyFDF is a complex array
    if not np.iscomplexobj(dirtyFDF):
        print "Err: the dirty FDF is not a complex array."
        sys.exit(1)
        
    # Check that phiArr is 1D and get its length
    if len(phiArr.shape) != 1:
        print "Err: the phi array is not a 1D array."
        sys.exit(1)
    nPhi = phiArr.shape[0]

    # Check that the lamSqArr is 1D and get its length
    if len(lamSqArr.shape) != 1:
        print "Err: the lamSqArr array is not a 1D array."
        sys.exit(1)
    nlamSq = lamSqArr.shape[0]

    # Check that phiArr and FDF arrays are the same length
    if nPhi != nFDF:
        print 'Err: the phiArr and dirty FDF are not the same length.'
        sys.exit(1)
    
    # If the RMSF has been passed in then check for correct formatting:
    #  - Twice the number of channels as dirtyFDF
    #  - Must be complex
    if not RMSFArr is None:

        # Check 1D
        if len(RMSFArr.shape) != 1:
            print "Err: input RMSF must be a 1D array."
            sys.exit(1)
        nRMSF = RMSFArr.shape[0]
        
        # Check complex
        if not np.iscomplexobj(RMSFArr):
            print "Err: the RMSF is not a complex array."
            sys.exit(1)
    
        # Check RMSF is at least double the FDF spectrum
        if not (nRMSF >= 2 * nFDF):
            print 'Err: the RMSF must be twice the length of the FDF.'
            sys.exit(1)

        # Check that phiSampArr is also present and the same length
        if RMSFphiArr is None:
            print 'Err: the phi sampling array must be passed with the RMSF.'
            sys.exit(1)
        nRMSFphi = RMSFphiArr.shape[0]
        if not nRMSF==nRMSFphi:
            print 'Err: the RMSF and phi sampling array must be equal length.'
            sys.exit(1)

        # If fwhmRMSF has not been passed in, fit for it now
        parms, status = fit_rmsf(RMSFphiArr, np.abs(RMSFArr))
        if status<1:
            print 'Err: failed to fit the RMSF.'
            sys.exit(1)
        else:
            fwhmRMSF = parms[2]
    
    # If the weight array has been passed in ...
    if not weight is None:
        
        uniformWt = False
        wtArr = np.array(weight, dtype=dType)
    
        # Check weightArr and lamSqArr have the same length
        if not wtArr.shape[0] == lamSqArr.shape[0]:
            print 'Err: the lamSqArr and weightArr are not the same length.'
            sys.exit(1)
            
    # or else use uniform weighting
    else:
        uniformWt = True
        wtArr = np.ones(lamSqArr.shape, dtype=dType)

    if doPlots:
        
        import pylab as pl
        import matplotlib as mpl
        
        # Setup the figure to track the clean
        fig = pl.figure()
        ax = fig.add_subplot(1,1,1)
        yMaxPlot = np.nanmax(np.abs(dirtyFDF))
        yMinPlot = np.nanmin(np.abs(dirtyFDF))
        yRangePlot = yMaxPlot - yMinPlot
        yMaxPlot +=  yRangePlot * 0.05
        yMinPlot -=  yRangePlot * 0.05
        ax.set_ylim([yMinPlot, yMaxPlot])
        fig.show()
        
    # Prerequisite calculations ----------------------------------------------#
    
    # Calculate the normalisation constant.
    # BdB Equations (24) and (38) give the inverse sum of the weights.
    K = 1.0 / np.nansum(wtArr)
    
    # Calculate the default lambda_0^2:
    # the weighted mean of the LambdaSq distribution (B&dB Eqn. 32).
    lam0Sq = K * np.nansum(wtArr * lamSqArr)

    # Calculate the RMSF if it has not been passed in 
    # Equation (26) OF BdB05
    if RMSFArr is None:
        RMSFArr, phiArr2, fwhmRMSF= get_RMSF(lamSqArr, phiArr, wtArr, lam0Sq)
    else:
        phiArr2 = RMSFphiArr
    
    # Find the index of the peak of the RMSF
    indxMaxRMSF = np.nanargmax(RMSFArr)
    
    # Initialise arrays to hold the residual FDF, clean components, clean FDF
    residFDF = dirtyFDF.copy()
    ccArr = np.zeros(phiArr.shape, dtype='complex')
    cleanFDF = np.zeros(phiArr.shape, dtype='complex')

    # HOGBOM CLEAN -----------------------------------------------------------#

    # Calculate the padding in the sampled RMSF
    # Assumes only integer shifts and symmetric
    phiMin = np.min(phiArr)
    nPhiPad = np.argwhere(phiArr2==phiMin)[0][0]

    # Main CLEAN loop
    iterCount = 0
    while ( np.max(np.abs(residFDF)) >= cutoff and iterCount < maxIter ):
        
        # Get the absolute peak channel and values of the residual FDF at peak
        indxPeakFDF = np.nanargmax(np.abs(residFDF))
        peakFDFvals = residFDF[indxPeakFDF]

        # What is the faraday depth at this channel?
        phiPeak = phiArr[indxPeakFDF]
        
        # A clean component (CC) is the max absolute amplitude * loop-gain
        #cc = gain * maxAbsResidFDF
        cc = gain * peakFDFvals
        ccArr[indxPeakFDF] += cc
        
        # At which channel is the CC located at in the RMSF?
        indxPeakRMSF = np.argwhere(phiArr2==phiPeak)[0][0]
        
        # Shift the RMSF in Faraday Depth so that its peak is centred above
        # this CC
        shiftedRMSFArr = np.roll(RMSFArr, indxPeakRMSF-indxMaxRMSF)
        
        # Clip the shifted RMSF to correspond to our FD range
        shiftedRMSFArr = shiftedRMSFArr[nPhiPad:-nPhiPad]
        
        # Subtract the product of the CC shifted RMSF from the residual FDF
        residFDF -= cc * shiftedRMSFArr
 
        # Restore the CC * a Gaussian 'main peak' RMSF to the cleaned FDF
        cleanFDF += cc * np.exp(-2.77258872224 *
                                np.power( (phiArr - phiPeak)/fwhmRMSF, 2.0)) 
        
        # Plot the progress of the clean
        if doPlots:
            ax.cla()
            ax.step(phiArr, np.abs(ccArr), color='g',marker='None',mfc='w',
                    mec='g', ms=10, label='none')
            ax.step(phiArr, np.abs(residFDF), color='r',marker='None',mfc='w',
                    mec='g', ms=10, label='none')
            #ax.step(phiArr, np.abs(shiftedRMSFArr), color='b',marker='None',
            #        mfc='w', mec='g', ms=10, label='none')
            ax.step(phiArr, np.abs(cleanFDF), color='k',marker='None',mfc='w',
                    mec='g', ms=10, label='none')
        
            ax.set_ylim(yMinPlot, yMaxPlot)
            pl.draw()

        # Iterate ...
        iterCount += 1
        print "Iteration %d" % iterCount

    # End clean loop ---------------------------------------------------------#

    # Restore the final residuals to the cleaned FDF
    cleanFDF += residFDF
    
    # Plot the final spectrum
    if doPlots:
        ax.cla()
        ax.step(phiArr, np.abs(dirtyFDF), color='grey',marker='None',mfc='w',
                mec='g', ms=10, label='none')
        ax.step(phiArr, np.abs(cleanFDF), color='k',marker='None',mfc='w',
                mec='g', ms=10, label='none')
        ax.step(phiArr, np.abs(ccArr), color='g',marker='None',mfc='w',
                mec='g', ms=10, label='none')
        #ax.step(phiArr, cleanFDF.real, color='r',marker='None',mfc='w',
        #        mec='g', ms=10, label='none')
        #ax.step(phiArr, cleanFDF.imag, color='b',marker='None',mfc='w',
        #        mec='g', ms=10, label='none')
            
        ax.set_xlabel('phi')
        ax.set_ylabel('Amplitude')
        ax.set_ylim(yMinPlot, yMaxPlot)
        fig.show()
        raw_input()

    return cleanFDF, ccArr, fwhmRMSF, iterCount


#-----------------------------------------------------------------------------#
def fits_make_lin_axis(head, axis=0):
    """
    Create an array containing the axis values, assuming a simple linear
    projection scheme. Axis selection is zero-indexed.
    """
    
    axis = int(axis)
    
    if head['NAXIS'] < axis + 1:
        return []
    
    i = str(int(axis) + 1)
    start = head['CRVAL' + i] + (1 - head['CRPIX' + i]) * head['CDELT' + i]
    stop = (head['CRVAL' + i] + (head['NAXIS' + i] + 1 - head['CRPIX' + i]) * 
            head['CDELT' + i])
    
    return np.arange(start, stop, head['CDELT' + i])
    

#-----------------------------------------------------------------------------#
def get_RMSF(lamSqArr, phiArr, wtArr=None, lam0Sq=None, double=True,
             dType='float32'):
    """
    Calculate the RMSF from inputs
    """
    

    # If the weight array has been passed in ...
    if not wtArr is None:

        uniformWt = False
        wtArr = np.array(wtArr, dtype=dType)

        # Check weightArr and lamSqArr have the same length
        if not wtArr.shape[0] == lamSqArr.shape[0]:
            print 'Err: the lamSqArr and weightArr are not the same length.'
            sys.exit(1)
            
    else:
        uniformWt = True
        wtArr = np.ones(lamSqArr.shape, dtype=dType)
            
    # B&dB equations (24) and (38) give the inverse sum of the weights
    K = 1.0 / np.nansum(wtArr)

    # Get the weighted mean of the LambdaSq distribution (B&dB Eqn. 32)
    if lam0Sq is None:
        lam0Sq = K * np.nansum(wtArr * lamSqArr)

    # For cleaning the RMSF should extend by 1/2 on each side in phi-space
    if double:
        nPhi = phiArr.shape[0]
        nExt = np.ceil(nPhi/2.0)
        resampIndxArr = np.arange(2.0 * nExt + nPhi) - nExt
        phiArr2 = extrap(resampIndxArr, np.arange(nPhi, dtype='int'), phiArr)
    else:
        phiArr2 = phiArr
        
    # Calculate the RM spread function
    a = (-2.0 * 1j * phiArr2)
    b = (lamSqArr - lam0Sq) 
    arg = wtArr * np.exp( np.outer(a, b) )
    RMSFArr = K * np.nansum(arg, 1)

    # Calculate or fit the main-lobe FWHM of the RMSF
    if uniformWt:            
        # B&dB Equation (61)
        fwhmRMSF = 2.0 * m.sqrt(3.0)/(np.nanmax(lamSqArr) -
                                      np.nanmin(lamSqArr))
    else:
        parms, status = fit_rmsf(phiArr2, np.abs(RMSFArr))
        if status<1:
            print 'Err: failed to fit the RMSF.'
            sys.exit(1)
        else:
            fwhmRMSF = parms[2]
            
    return RMSFArr, phiArr2, fwhmRMSF


#-----------------------------------------------------------------------------#
def extrap(x, xp, yp):
    """
    Wrapper to allow np.interp to linearly extrapolate at function ends.
    
    np.interp function with linear extrapolation
    http://stackoverflow.com/questions/2745329/how-to-make-scipy-interpolate
    -give-a-an-extrapolated-result-beyond-the-input-ran

    """
    y = np.interp(x, xp, yp)
    y = np.where(x < xp[0], yp[0]+(x-xp[0])*(yp[0]-yp[1])/(xp[0]-xp[1]), y)
    y = np.where(x > xp[-1], yp[-1]+(x-xp[-1])*(yp[-1]-yp[-2])/(xp[-1]-xp[-2]),
                 y)
    return y


#-----------------------------------------------------------------------------#
def fit_rmsf(xData, yData, thresh=0.3):
    """
    Fit the main lobe of the RMSF with a Gaussian function. Sidelobes beyond
    the first null are masked out.
    """

    # Detect the peak and mask off the sidelobes
    msk = detect_peak(yData, thresh)
    validIndx = np.where(msk==1.0)
    xData = xData[validIndx]
    yData = yData[validIndx]
    
    # Estimate starting parameters
    a = 1.0
    b = xData[np.argmax(yData)]
    w = np.nanmax(xData)-np.nanmin(xData)
    
    # Estimate starting parameters
    inParms=[ {'value': a, 'fixed':False, 'parname': 'amp'},
              {'value': b, 'fixed':False, 'parname': 'offset'},
              {'value': w, 'fixed':False, 'parname': 'width'}]

    # Function which returns another function to evaluate a Gaussian
    def gauss1D(p):
        a, b, w = p
        gfactor = 2.0 * m.sqrt(2.0 * m.log(2.0))
        s = w / gfactor    
        def rfunc(x):
            y = a * np.exp(-(x-b)**2.0 /(2.0 * s**2.0))
            return y
        return rfunc
    
    # Function to evaluate the difference between the model and data.
    # This is minimised in the least-squared sense by the fitter
    def errFn(p, fjac=None):
        status = 0
        return status, gauss1D(p)(xData) - yData
    
    # Use mpfit to perform the fitting
    mp = mpfit(errFn, parinfo=inParms, quiet=True)
    coeffs = mp.params

    return mp.params, mp.status


#-----------------------------------------------------------------------------#
def detect_peak(a, thresh=0.3):
    """
    Detect the extent of the peak in the array by looking for where the slope
    changes to flat. The highest peak is detected and data and followed until
    the slope flattens to a threshold.
    """
    
    iPk = np.argmax(a)
    d = np.diff(a)
    g1 = np.gradient(a)
    g2 = np.gradient(g1)
    
    threshPos = np.nanmax(d) * thresh
    threshNeg = -1 * threshPos

    # Start searching away from the peak zone
    g2L = np.flipud(g2[:iPk])
    g2R = g2[iPk+1:]
    iL = iPk - np.min(np.where(g2L>=0))
    iR = iPk + np.min(np.where(g2R>=0)) + 1
    g1[iL:iR] = np.nan
    
    # Search for the threshold crossing point
    g1L = np.flipud(g1[:iPk])
    g1R = g1[iPk+1:]
    iL = iPk - np.min(np.where(g1L<=threshPos))
    iR = iPk + np.min(np.where(g1R>=threshNeg))
    msk = np.zeros_like(a)
    msk[iL:iR] = 1
    
    return msk
    

