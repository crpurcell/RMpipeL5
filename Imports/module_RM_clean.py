#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     module_RM_clean.py                                                #
#                                                                             #
# PURPOSE:  Perform Hogbom RM-clean on dirty FDFs in the PPC.                 #
#                                                                             #
# MODIFIED: 19-November-2015 by C. Purcell                                    #
#                                                                             #
# TODO:                                                                       #
#  * Catch errors in non-fatal way                                            #
#  * If doOverwrite is not set, read and return the existing result           #
#  * Do we need to account for flagged channels?                              #
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

import astropy.io.fits as pf

from util_PPC import log_wr
from util_PPC import log_fail
from util_PPC import fail_not_exists
from util_PPC import DataManager

from util_RM import do_rmclean

# Constants
C = 2.99792458e8

#-----------------------------------------------------------------------------#
def mod_do_RMclean(specRec, sessionPath, cleanCutoff_sigma=5,
                   maxIter=1000, gain=0.1, doOverwrite=False, LF=None):

    # Default logging to STDOUT
    if LF is None:
        LF = sys.stdout
        
    # Check required directories exist
    fail_not_exists(sessionPath, 'directory', LF)
    dataPath = sessionPath + "/OUT"
    fail_not_exists(dataPath, 'directory', LF)
    inParmFile = sessionPath + "/inputs.config"
    fail_not_exists(inParmFile, "file", LF)
    
    # Create a recArray to store the clean properties
    dType = [('uniqueName', 'a20'),
             ('nIterDone', 'i8'),
             ('cleanCutoff_sigma', 'f8'),
             ('cleanCutoff_Jybm', 'f8')]
    cleanRec = np.zeros(len(specRec), dtype=dType)
    
    # Create a DataManager object to access the stored data products
    dataMan = DataManager(sessionPath, calcParms=False)
    cleanCutoff_sigma = float(dataMan.pDict["cleanCutoff_sigma"])
    maxIter = int(dataMan.pDict["maxCleanIter"])
    gain = float(dataMan.pDict["gain"])
    
    # Loop through the catalogue entries
    log_wr(LF, '\nPerforming RM-clean on the catalogue entries ...')
    for i in range(len(specRec)):
        uniqueName = specRec[i]['uniqueName']        
        log_wr(LF, "\nProcessing entry %d: '%s'." %  (i+1, uniqueName))

        # Read in the dirty FDF, RMSF, frequency and weight arrays
        phiArr, dirtyFDF = dataMan.get_dirtyFDF_byname(uniqueName)
        RMSFphiArr, RMSFArr = dataMan.get_RMSF_byname(uniqueName)
        freqArr_Hz, weightArr = dataMan.get_freqweight_byname(uniqueName)
        nFreqChan = len(freqArr_Hz) 

        # Calculate the lamSqArr
        lamArr_m = C / freqArr_Hz
        lamSqArr_m2 = np.power(lamArr_m, 2.0)

        # Calculate the clean cutoff
        cleanCutoff_Jybm = (cleanCutoff_sigma *
                            specRec[i]['rmsMedQUAvg_Jybm'] / m.sqrt(nFreqChan))
        LF.write('> Using cutoff = %.3f mJy.\n' % (cleanCutoff_Jybm*1e3))
                
        # Run the RM-clean procedure
        cleanFDF, ccModel, fwhmRMSF, nIter = do_rmclean(dirtyFDF,
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
        dataMan.put_cleanFDF_byname(uniqueName, CC=ccModel, cleanFDF=cleanFDF)
        log_wr(LF, '> Clean FDF and CC model saved to FITS files.')
        
        # Write flags and metadata to the record array
        cleanRec[i]['uniqueName'] = specRec[i]['uniqueName']
        cleanRec[i]['nIterDone'] = nIter
        cleanRec[i]['cleanCutoff_sigma'] = cleanCutoff_sigma
        cleanRec[i]['cleanCutoff_Jybm'] = cleanCutoff_Jybm
        
    return cleanRec
