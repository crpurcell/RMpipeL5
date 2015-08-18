#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     module_RM_synthesis.py                                            #
#                                                                             #
# PURPOSE:  Function perform RM-synthesis on spectra in the PPC.              #
#                                                                             #
# MODIFIED: 18-Aug-2015 by C. Purcell                                         #
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
from util_RM import do_rmsynth

import astropy.io.fits as pf

# Constants
C = 2.99792458e8


#-----------------------------------------------------------------------------#
def mod_do_RMsynth(specRec, dataInDir, outDataDir, phiArr,
                   weightType='Variance', doOverwrite=False, LF=None):

    # Default logging to STDOUT
    if LF is None:
        LF = sys.stdout
        
    # Check the data directorys exist
    fail_not_exists(dataInDir, 'directory', LF)
    if not os.path.exists(outDataDir):
        os.mkdir(outDataDir)

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
        
    # Loop through the catalogue entries
    log_wr(LF, '\nPerforming RM-synthesis on the catalogue entries ...')
    for i in range(len(specRec)):
        
        log_wr(LF, "\nProcessing entry %d: '%s'." %
               (i+1, specRec[i]['uniqueName']))
        
        try:
            
            # Read in the I-model, Q & U spectra & noise spectra
            uniqueName = specRec[i]['uniqueName']
            dataFile = dataInDir +  '/' + uniqueName +  '_specI.fits'
            HDULst = pf.open(dataFile, "readonly", memmap=True)
            freqArr_Hz = HDULst[2].data["freq"]
            rmsIArr_Jy = HDULst[2].data["rms"]
            modIArr_Jy = HDULst[2].data["model"]
            HDULst.close()
            dataFile = dataInDir +  '/' + uniqueName +  '_specQ.fits'
            HDULst = pf.open(dataFile, "readonly", memmap=True)
            QArr_Jy = HDULst[2].data["src"]
            rmsQArr_Jy = HDULst[2].data["rms"]
            HDULst.close()
            dataFile = dataInDir +  '/' + uniqueName +  '_specU.fits'
            HDULst = pf.open(dataFile, "readonly", memmap=True)
            UArr_Jy = HDULst[2].data["src"]
            rmsUArr_Jy = HDULst[2].data["rms"]
            HDULst.close()
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
        
            # Calculate the lambda sampling from the frequency array
            lamArr_m = C / freqArr_Hz
            lamSqArr_m2 = np.power(lamArr_m, 2.0)
        
            # Perform RM-synthesis on the spectrum
            # Returned dirty FDF has amplitude expressed as polarised fraction
            tmp = do_rmsynth(specQfrac, specUfrac, lamSqArr_m2, phiArr, 
                             weightArr)
            dirtyFDF, [phiSampArr, RMSFArr], lam0Sq, fwhmRMSF = tmp
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
            outFits = outDataDir + '/' + specRec[i]['uniqueName'] + \
                      "_RMSynth.fits"
            hdu0 = pf.PrimaryHDU(header=pf.Header())
            
            col1 = pf.Column(name="freq", format="f4", array=freqArr_Hz)
            col2 = pf.Column(name="lamsq", format="f4", array=lamSqArr_m2)
            col3 = pf.Column(name="weight", format="f4", array=weightArr)
            hdu1 = pf.BinTableHDU.from_columns([col1, col2, col3])
            
            col1 = pf.Column(name="phiSamp", format="f4", array=phiSampArr)
            col2 = pf.Column(name="RMSFreal", format="f4", array=RMSFArr.real)
            col3 = pf.Column(name="RMSFimag", format="f4", array=RMSFArr.imag)
            hdu2 = pf.BinTableHDU.from_columns([col1, col2, col3])
            
            col1 = pf.Column(name="phi", format="f4", array=phiArr)
            col2 = pf.Column(name="FDFreal", format="f4", array=dirtyFDF.real)
            col3 = pf.Column(name="FDFimag", format="f4", array=dirtyFDF.imag)
            hdu3 = pf.BinTableHDU.from_columns([col1, col2, col3])

            col1 = pf.Column(name="CC", format="f4",
                             array=np.zeros_like(dirtyFDF.real))
            col2 = pf.Column(name="FDFreal", format="f4",
                             array=np.zeros_like(dirtyFDF.real))
            col3 = pf.Column(name="FDFimag", format="f4",
                             array=np.zeros_like(dirtyFDF.real))
            hdu4 = pf.BinTableHDU.from_columns([col1, col2, col3])
            
            hduLst = pf.HDUList([hdu0, hdu1, hdu2, hdu3, hdu4])
            hduLst.writeto(outFits, output_verify="fix", clobber=True)
            hduLst.close()

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
