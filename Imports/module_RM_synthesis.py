#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     module_RM_synthesis.py                                            #
#                                                                             #
# PURPOSE:  Function perform RM-synthesis on spectra in the PPC.              #
#                                                                             #
# MODIFIED: 13-October-2015 by C. Purcell                                     #
#                                                                             #
# TODO:                                                                       #
#   * If doOverwrite is not set, read and return the existing result          #
#                                                                             #
#=============================================================================#
import os
import sys
import math as m
import numpy as np
import traceback

from util_PPC import fail_not_exists
from util_PPC import log_wr
from util_PPC import log_fail
from util_PPC import poly5
from util_PPC import DataManager
from util_RM import do_rmsynth

import astropy.io.fits as pf

# Constants
C = 2.99792458e8


#-----------------------------------------------------------------------------#
def mod_do_RMsynth(specRec, sessionPath, doOverwrite=False, LF=None):
    
    # Default logging to STDOUT
    if LF is None:
        LF = sys.stdout
        
    # Check required directories and files exist
    fail_not_exists(sessionPath, 'directory', LF)
    dataPath = sessionPath + "/OUT"
    fail_not_exists(dataPath, 'directory', LF)
    inParmFile = sessionPath + "/inputs.config"
    fail_not_exists(inParmFile, "file", LF)
    
    # Create a recArray to store RMSF and FDF properties
    dType = [('uniqueName', 'a20'),
             ('weightType', 'a20'),
             ('lam0Sq_m2', 'f8'),
             ('freq0_Hz', 'f8'),
             ('nPhiChan', 'i8'),
             ('fwhmRMSF', 'f8'),
             ('deltaPhiChan_rm2', 'f8'),
             ('phiCentre_rm2', 'f8'),
             ('status', 'i8')]
    rmsfRec = np.zeros(len(specRec), dtype=dType)
        
    # Create a DataManager to access the stored data and inputs
    # Assumption: all data have the same frequency sampling
    dataMan = DataManager(sessionPath, calcParms=True)
    phiArr = dataMan.pDict["phiArr_radm2"]
    weightType = dataMan.pDict["weightType"]
    lamSqArr_m2 = dataMan.pDict["lambdaSqArr_m2"]
    
    # Loop through the catalogue entries
    log_wr(LF, '\nPerforming RM-synthesis on the catalogue entries ...')
    for i in range(len(specRec)):
        uniqueName = specRec[i]['uniqueName']
        log_wr(LF, "\nProcessing entry %d: '%s'." %(i+1, uniqueName))
        
        try:
            
            # Read in the I-model, Q & U spectra & noise spectra
            freqArr_Hz, IArr_Jy, rmsIArr_Jy = \
                        dataMan.get_specI_byname(uniqueName)
            dummy, modIArr_Jy = \
                   dataMan.get_modI_byname(uniqueName, getStored=True)
            dummy, QArr_Jy, rmsQArr_Jy = dataMan.get_specQ_byname(uniqueName)
            dummy, UArr_Jy, rmsUArr_Jy = dataMan.get_specU_byname(uniqueName)
            PIArr_Jy = np.sqrt( np.power(QArr_Jy, 2.0) +
                                np.power(UArr_Jy, 2.0) )
            rmsSpecQUAvg = (rmsQArr_Jy + rmsUArr_Jy)/2.0
                
            # Calculate the weighting as 1/sigma^2 or 1
            if weightType=='Variance':
                weightArr = 1.0 / np.power(rmsSpecQUAvg, 2.0)
            else:
                weightArr = np.ones(rmsSpecQUAvg.shape)
            log_wr(LF,"> Weight type is '%s'." % weightType)
        
            # Form the Q/I & U/I spectra (divide by the model Stokes I)
            specQfrac = QArr_Jy / modIArr_Jy
            specUfrac = UArr_Jy / modIArr_Jy
            specPIfrac = PIArr_Jy / modIArr_Jy
            log_wr(LF, '> Q & U data divided by the Stokes I model.')
        
            # Perform RM-synthesis on the spectrum
            # Returned dirty FDF has amplitude expressed as polarised fraction
            dirtyFDF, [phiSampArr, RMSFArr], lam0Sq, fwhmRMSF = \
               do_rmsynth(specQfrac, specUfrac, lamSqArr_m2, phiArr, weightArr)
            log_wr(LF,'> RM-synthesis completed on Q/I and U/I spectra.')

            # Determine the Stokes I value at lam0Sq from the I polynomial
            # Multiply the dirty FDF by Ifreq0 to recover the PI in Jy
            polyCoeffs = specRec[i]['coeffPolyIspec'].split(',')
            polyCoeffs = [float(x) for x in polyCoeffs]
            freq0_Hz = C / m.sqrt(lam0Sq)
            Ifreq0_Jybm = poly5(polyCoeffs)(freq0_Hz / 1e9)
            dirtyFDF *= Ifreq0_Jybm
            log_wr(LF, '> Dirty FDF recast into input units.')

            # Create a FITS file to store the RM-synthesis results
            dataMan.create_RMSFfits_byname(uniqueName, freqArr_Hz, weightArr,
                                   phiSampArr, RMSFArr, phiArr, dirtyFDF)
            log_wr(LF, '> Dirty FDF saved to FITS file.')
            
        except Exception:
            rmsfRec[i]['uniqueName'] = specRec[i]['uniqueName']
            rmsfRec[i]['status'] = 0
            log_wr(LF, "Err: failed to process %s." % specRec[i]['uniqueName'])
            log_wr(LF, traceback.format_exc())
            continue
        
        # Write flags and metadata to the catalogue recarray
        rmsfRec[i]['uniqueName'] = specRec[i]['uniqueName']
        rmsfRec[i]['weightType'] = weightType
        rmsfRec[i]['lam0Sq_m2'] = lam0Sq
        rmsfRec[i]['freq0_Hz'] = freq0_Hz
        rmsfRec[i]['nPhiChan'] = int(len(phiArr))
        rmsfRec[i]['fwhmRMSF'] = fwhmRMSF
        rmsfRec[i]['deltaPhiChan_rm2'] = np.nanmin(np.diff(phiArr))
        rmsfRec[i]['phiCentre_rm2'] = np.median(phiArr)
        rmsfRec[i]['status'] = 1

    return rmsfRec
