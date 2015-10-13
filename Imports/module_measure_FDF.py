#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     module_measure_FDF.py                                             #
#                                                                             #
# PURPOSE:  Make measurements on a catalogue of Faraday dispersion functions. #
#                                                                             #
# MODIFIED: 13-October-2015 by C. Purcell                                     #
#                                                                             #
# TODO:                                                                       #
#  * Set a flag if the peak is within RMSF_FWHM/2 of the FDF edge.            #
#                                                                             #
#=============================================================================#

import os
import sys
import math as m
import traceback
import numpy as np

import astropy.io.fits as pf

from util_PPC import log_fail
from util_PPC import fail_not_exists
from util_PPC import log_wr
from util_PPC import DataManager

# Constants
C = 2.99792458e8


#-----------------------------------------------------------------------------#
def mod_measure_FDF(catRec, sessionPath, dirty=True, LF=None):
    
    # Default logging to STDOUT
    if LF is None:
        LF = sys.stdout

    # Check required directories exist
    fail_not_exists(sessionPath, 'directory', LF)
    dataPath = sessionPath + "/OUT"
    fail_not_exists(dataPath, 'directory', LF)
    
    # Create a numpy recArray to store RMSF and FDF properties
    dType = [('uniqueName', 'a20'),
             ('phiPeakPIchan_rm2', 'f8'),
             ('dPhiPeakPIchan_rm2', 'f8'),
             ('ampPeakPIchan_Jybm', 'f8'),
             ('ampPeakPIchanEff_Jybm', 'f8'),
             ('dAmpPeakPIchan_Jybm', 'f8'),
             ('snrPIchan', 'f8'),
             ('indxPeakPIchan', 'i8'),
             ('peakFDFimagChan', 'f8'),
             ('peakFDFrealChan', 'f8'),
             ('polAngleChan_deg', 'f8'),
             ('dPolAngleChan_deg', 'f8'),
             ('polAngle0Chan_deg', 'f8'),
             ('dPolAngle0Chan_deg', 'f8'),
             ('phiPeakPIfit_rm2', 'f8'),
             ('dPhiPeakPIfit_rm2', 'f8'),
             ('ampPeakPIfit_Jybm', 'f8'),
             ('ampPeakPIfitEff_Jybm', 'f8'),
             ('dAmpPeakPIfit_Jybm', 'f8'),
             ('snrPIfit', 'f8'),
             ('indxPeakPIfit', 'i8'),
             ('peakFDFimagFit', 'f8'),
             ('peakFDFrealFit', 'f8'),
             ('polAngleFit_deg', 'f8'),
             ('dPolAngleFit_deg', 'f8'),
             ('polAngle0Fit_deg', 'f8'),
             ('dPolAngle0Fit_deg', 'f8'),
             ('thresholdSignalPI', 'f8'),
             ('detectionF', 'i8'),
             ('status', 'i8')]
    fdfRec = np.zeros(len(catRec), dtype=dType)

    # Create a DataManager object to access the stored data products
    dataMan = DataManager(sessionPath, calcParms=False)
    thresholdSignalPI = float(dataMan.pDict["thresholdSignalPI_sigma"])
    
    # Loop through the catalogue entries
    log_wr(LF, '\nMeasuring the properties of the FDF catalogue entries ...')
    for i in range(len(catRec)):
        uniqueName = catRec[i]['uniqueName']
        
        log_wr(LF, "\nProcessing entry %d: '%s'." % (i+1, uniqueName))
        try:
            
            # Read in the complex FDF
            if dirty is False:
                phiArr, FDF = dataMan.get_cleanFDF_byname(uniqueName)
            else:
                phiArr, FDF = dataMan.get_dirtyFDF_byname(uniqueName)

            # Read in the frequency array and calculate the lamSqArr
            freqArr_Hz, weightArr = dataMan.get_freqweight_byname(uniqueName)
            lamArr_m = C / freqArr_Hz
            lamSqArr_m2 = np.power(lamArr_m, 2.0)
        
            # Perform the measurements
            mDict = measure_FDF_parms(FDF,
                                      phiArr,
                                      catRec[i]['fwhmRMSF'],
                                      lamSqArr_m2,
                                      catRec[i]['lam0Sq_m2'],
                                      catRec[i]['rmsMedQUAvg_Jybm'],
                                      snrDoBiasCorrect=5.0)
            # Check for a detection
            if mDict['snrPIchan'] >= thresholdSignalPI:
                mDict['detectionF'] = 1
                log_wr(LF, '> Signal detected at SNR = %.1f' % 
                       mDict['snrPIchan'])
            else:
                mDict['detectionF'] = 0
                log_wr(LF, '> No signal detected (SNR = %.1f).' %
                       mDict['snrPIchan'])

            # TODO: Check if the peak is within half the RMSF of the edge
            #       set a flag saying detection is near edge

        
            # Write the measurements to the catalogue
            for key, val in mDict.iteritems():
                try:
                    fdfRec[i][key] = val
                except Exception:
                    pass
            fdfRec[i]['uniqueName'] = catRec[i]['uniqueName']
            fdfRec[i]['thresholdSignalPI'] = thresholdSignalPI
            fdfRec[i]['status'] = 1
        except Exception:
            fdfRec[i]['uniqueName'] = catRec[i]['uniqueName']
            fdfRec[i]['thresholdSignalPI'] = thresholdSignalPI
            fdfRec[i]['status'] = 0
            log_wr(LF, "Err: failed to process %s." % catRec[i]['uniqueName'])
            log_wr(LF, traceback.format_exc())
            continue

    return fdfRec


#-----------------------------------------------------------------------------#
def measure_FDF_parms(FDF, phiArr, fwhmRMSF, lamSqArr_m2, lam0Sq, dQU,
                      snrDoBiasCorrect=5.0):
    """
    Measure standard parameters of a Faraday Dispersion Function.
    Currently this function assumes that the noise levels in the Stokes Q
    and U spectra are the same.
    """
    
    # Determine the peak channel in the FDF, its amplitude and RM
    absFDF = np.abs(FDF)
    ampPeakPIchan = np.nanmax(absFDF)
    indxPeakPIchan = np.nanargmax(absFDF)
    phiPeakPIchan = phiArr[indxPeakPIchan]
    dPhiPeakPIchan = fwhmRMSF * dQU / (2.0 * ampPeakPIchan)
    snrPIchan = ampPeakPIchan / dQU
        
    # Correct the peak for polarisation bias (POSSUM report 11)
    ampPeakPIchanEff = None
    if snrPIchan >= snrDoBiasCorrect:
        ampPeakPIchanEff = m.sqrt(ampPeakPIchan**2.0 - 2.3 * dQU**2.0)

    # Calculate the polarisation angle from the channel
    peakFDFimagChan = FDF.imag[indxPeakPIchan]
    peakFDFrealChan = FDF.real[indxPeakPIchan]
    polAngleChan_deg = 0.5 * m.degrees(m.atan2(peakFDFimagChan,
                                         peakFDFrealChan))
    dPolAngleChan_deg = m.degrees(dQU**2.0 / (4.0 * ampPeakPIchan**2.0))

    # Calculate the derotated polarisation angle and uncertainty
    polAngle0Chan_deg = m.degrees(m.radians(polAngleChan_deg) -
                                  phiPeakPIchan * lam0Sq)
    nChansGood = np.sum(np.where(lamSqArr_m2==lamSqArr_m2, 1.0, 0.0))
    varLamSqArr_m2 = (np.sum(lamSqArr_m2**2.0) -
                      np.sum(lamSqArr_m2)**2.0/nChansGood) / (nChansGood-1)
    dPolAngle0Chan_rad = \
        m.sqrt( dQU**2.0 / (4.0*(nChansGood-2.0)*ampPeakPIchan**2.0) *
                 ((nChansGood-1)/nChansGood + lam0Sq**2.0/varLamSqArr_m2) )
    dPolAngle0Chan_deg = m.degrees(dPolAngle0Chan_rad)

    # Determine the peak in the FDF, its amplitude and Phi using a
    # 3-point parabolic interpolation
    phiPeakPIfit = None
    dPhiPeakPIfit = None
    ampPeakPIfit = None
    dAmpPeakPIfit = None
    snrPIfit = None
    ampPeakPIfitEff = None
    indxPeakPIfit = None
    peakFDFimagFit = None 
    peakFDFrealFit = None 
    polAngleFit_deg = None
    dPolAngleFit_deg = None
    polAngle0Fit_deg = None
    dPolAngle0Fit_deg = None
    
    if indxPeakPIchan > 0 and indxPeakPIchan < len(FDF)-1:
        phiPeakPIfit, ampPeakPIfit = \
                      calc_parabola_vertex(phiArr[indxPeakPIchan-1],
                                           absFDF[indxPeakPIchan-1],
                                           phiArr[indxPeakPIchan],
                                           absFDF[indxPeakPIchan],
                                           phiArr[indxPeakPIchan+1],
                                           absFDF[indxPeakPIchan+1])
        dPhiPeakPIfit = fwhmRMSF * dQU / (2.0 * ampPeakPIfit)
        dAmpPeakPIfit = m.sqrt(2.0 * ampPeakPIfit**2.0 / dQU**2.0)
        snrPIfit = ampPeakPIfit / dQU

        # Correct the peak for polarisation bias (POSSUM report 11)
        if snrPIfit >= snrDoBiasCorrect:
            ampPeakPIfitEff = m.sqrt(ampPeakPIfit**2.0
                                     - 2.3 * dQU**2.0)
            
        # Calculate the polarisation angle from the fitted peak        
        indxPeakPIfit = np.interp(phiPeakPIfit, phiArr,
                                  np.arange(phiArr.shape[-1], dtype='f4'))
        peakFDFimagFit = np.interp(phiPeakPIfit, phiArr, FDF.imag)
        peakFDFrealFit = np.interp(phiPeakPIfit, phiArr, FDF.real)
        polAngleFit_deg = 0.5 * m.degrees(m.atan2(peakFDFimagFit,
                                                  peakFDFrealFit))
        dPolAngleFit_deg = m.degrees(dQU**2.0 / (4.0 * ampPeakPIfit**2.0))

        # Calculate the derotated polarisation angle and uncertainty
        polAngle0Fit_deg = m.degrees(m.radians(polAngleFit_deg) -
                                     phiPeakPIfit * lam0Sq)
        dPolAngle0Fit_rad = \
            m.sqrt( dQU**2.0 / (4.0*(nChansGood-2.0)*ampPeakPIfit**2.0) *
                    ((nChansGood-1)/nChansGood + lam0Sq**2.0/varLamSqArr_m2) )
        dPolAngle0Fit_deg = m.degrees(dPolAngle0Fit_rad)

    # Store the measurements in a dictionary and return
    mDict = {'phiPeakPIchan_rm2': phiPeakPIchan,
             'dPhiPeakPIchan_rm2': dPhiPeakPIchan,
             'ampPeakPIchan_Jybm': ampPeakPIchan,
             'ampPeakPIchanEff_Jybm': ampPeakPIchanEff,
             'dAmpPeakPIchan_Jybm': dQU,
             'snrPIchan': snrPIchan,
             'indxPeakPIchan': indxPeakPIchan,
             'peakFDFimagChan': peakFDFimagChan,
             'peakFDFrealChan': peakFDFrealChan,
             'polAngleChan_deg': polAngleChan_deg,
             'dPolAngleChan_deg': dPolAngleChan_deg,
             'polAngle0Chan_deg': polAngle0Chan_deg,
             'dPolAngle0Chan_deg': dPolAngle0Chan_deg,
             'phiPeakPIfit_rm2': phiPeakPIfit,
             'dPhiPeakPIfit_rm2': dPhiPeakPIfit,
             'ampPeakPIfit_Jybm': ampPeakPIfit,
             'ampPeakPIfitEff_Jybm': ampPeakPIfitEff,
             'dAmpPeakPIfit_Jybm': dQU,
             'snrPIfit': snrPIfit,
             'indxPeakPIfit': indxPeakPIfit,
             'peakFDFimagFit': peakFDFimagFit,
             'peakFDFrealFit': peakFDFrealFit,
             'polAngleFit_deg': polAngleFit_deg,
             'dPolAngleFit_deg': dPolAngleFit_deg,
             'polAngle0Fit_deg': polAngle0Fit_deg,
             'dPolAngle0Fit_deg': dPolAngle0Fit_deg}

    return mDict


#-----------------------------------------------------------------------------#
def calc_parabola_vertex(x1, y1, x2, y2, x3, y3):
    """
    Calculate the vertex of a parabola given three adjacent points.
    """
    
    D = (x1 - x2) * (x1 - x3) * (x2 - x3)
    A = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / D
    B = (x3*x3 * (y1 - y2) + x2*x2 * (y3 - y1) + x1*x1 * (y2 - y3)) / D
    C = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 *
         (x1 - x2) * y3) / D

    xv = -B / (2.0 * A)
    yv = C - B * B / (4.0 * A)

    return xv, yv
