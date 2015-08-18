#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     module_RM_clean.py                                                #
#                                                                             #
# PURPOSE:  Perform Hogbom RM-clean on dirty FDFs in the PPC.                 #
#                                                                             #
# MODIFIED: 18-August-2015 by C. Purcell                                      #
#                                                                             #
# TODO:                                                                       #
#  * Catch errors in non-fatal way                                            #
#  * If doOverwrite is not set, read and return the existing result           #
#                                                                             #
#=============================================================================#

import os
import sys
import math as m
import numpy as np

import astropy.io.fits as pf

from util_PPC import log_wr
from util_PPC import log_fail
from util_PPC import fail_not_exists
from util_RM import do_rmclean

# Constants
C = 2.99792458e8

#-----------------------------------------------------------------------------#
def mod_do_RMclean(specRec, dataInDir, outDataDir, cleanCutoff_sigma=5,
                   maxIter=1000, gain=0.1, doOverwrite=False, LF=None):

    # Default logging to STDOUT
    if LF is None:
        LF = sys.stdout
        
    # Create the output data directory
    fail_not_exists(dataInDir, 'directory', LF)
    if not os.path.exists(outDataDir):
        os.mkdir(outDataDir)
        
    # Create a recArray to store RMSF and FDF properties
    dType = [('uniqueName', 'a20'),
             ('nIterDone', 'i8'),
             ('cleanCutoff_sigma', 'f8'),
             ('cleanCutoff_Jybm', 'f8')]
    cleanRec = np.zeros(len(specRec), dtype=dType)
    
    # Loop through the catalogue entries
    log_wr(LF, '\nPerforming RM-clean on the catalogue entries ...')
    for i in range(len(specRec)):
        uniqueName = specRec[i]['uniqueName']
        
        log_wr(LF, "\nProcessing entry %d: '%s'." %
               (i+1, uniqueName))

        # Read in the dirty FDF
        dataFile = dataInDir +  '/' + uniqueName +  "_RMSynth.fits"
        HDULst = pf.open(dataFile, "update", memmap=True)
        phiArr = HDULst[3].data["phi"]
        FDFreal = HDULst[3].data["FDFreal"]
        FDFimag = HDULst[3].data["FDFimag"]
        dirtyFDF = (FDFreal + 1j * FDFimag)

        # Read in the RMSF, frequency and weight arrays
        freqArr_Hz = HDULst[1].data["freq"]
        weightArr = HDULst[1].data["weight"]
        RMSFphiArr = HDULst[2].data["phiSamp"]
        RMSFreal = HDULst[2].data["RMSFreal"]
        RMSFimag = HDULst[2].data["RMSFimag"]
        RMSFArr = (RMSFreal + 1j * RMSFimag)
        nFreqChan = len(freqArr_Hz) # REPLACE WITH NVALID?

        # Calculate the lamSqArr
        lamArr_m = C / freqArr_Hz
        lamSqArr_m2 = np.power(lamArr_m, 2.0)

        # Calculate the clean cutoff
        cleanCutoff_Jybm = (cleanCutoff_sigma *
                          specRec[i]['rmsMedQUAvg_Jybm'] / m.sqrt(nFreqChan))
        LF.write('> Using cutoff = %.3f mJy.\n' % (cleanCutoff_Jybm*1000))
                
        # Run the RM-clean procedure
        cleanFDF, ccModel, fwhmRMSF, nIter = \
                  do_rmclean(dirtyFDF,
                             phiArr,
                             lamSqArr_m2,
                             cleanCutoff_Jybm,
                             maxIter,
                             gain,
                             weightArr,
                             RMSFArr,
                             RMSFphiArr,
                             doPlots=False)
        LF.write('> CLEANed dirty FDF in %d iterations.\n' % nIter)

        # Save the clean FDF and CC model to the FITS file
        HDULst[4].data["CC"] = np.abs(ccModel)
        HDULst[4].data["FDFreal"] = cleanFDF.real
        HDULst[4].data["FDFimag"] = cleanFDF.imag
        HDULst.close()
        log_wr(LF, '> Clean FDF and CC model saved to FITS files.')

        # Write flags and metadata to the record array
        cleanRec[i]['uniqueName'] = specRec[i]['uniqueName']
        cleanRec[i]['nIterDone'] = nIter
        cleanRec[i]['cleanCutoff_sigma'] = cleanCutoff_sigma
        cleanRec[i]['cleanCutoff_Jybm'] = cleanCutoff_Jybm

    return cleanRec
