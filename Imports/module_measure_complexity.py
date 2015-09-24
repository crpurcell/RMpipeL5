#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     module_measure_complexity.py                                      #
#                                                                             #
# PURPOSE:  Assess the complexity of a polarised dataset using the Stokes     #
#           parameter spectra and the Faraday dispersion function.            #
#                                                                             #
# MODIFIED: 24-September-2015 by C. Purcell                                   #
#                                                                             #
# TODO:                                                                       #
#                                                                             #
#=============================================================================#

import os
import sys
import math as m
import traceback
import numpy as np

import astropy.io.fits as pf

from util_PPC import DataManager
from util_PPC import poly5
from util_PPC import calc_mom2_FDF
from util_PPC import log_fail
from util_PPC import fail_not_exists
from util_PPC import log_wr
from util_PPC import create_pqu_spectra_RMthin

# Constants
C = 2.99792458e8


#-----------------------------------------------------------------------------#
def mod_assess_complexity(catRec, sessionPath, threshC1=3.0, threshC2=3.0,
                          threshC23=3.0, LF=None):
    
    # Default logging to STDOUT
    if LF is None:
        LF = sys.stdout        
    fail_not_exists(sessionPath, 'directory', LF)

    # Create a numpy recArray to store the complexity measurements
    dType = [('uniqueName', 'a20'),
             ('complexM1', 'f8'),
             ('complexM2', 'f8'),
             ('complexM3', 'f8'),
             ('status', 'f8')]
    compRec = np.zeros(len(catRec), dtype=dType)
    
    # Create a DataManager object to access the stored data products
    dataMan = DataManager(sessionPath, calcParms=False)

    # Loop through the catalogue entries
    log_wr(LF, '\nAssessing complexity in the catalogue entries ...')
    for i in range(len(catRec)):
        uniqueName = catRec[i]['uniqueName']
        log_wr(LF, "\nProcessing entry %d: '%s'." % (i+1, uniqueName))
        try:
            
            # Read the data from disk
            freqArr_Hz, IArr, rmsIArr = dataMan.get_specI_byname(uniqueName)
            dummy, modIArr = dataMan.get_modI_byname(uniqueName,
                                                     oversample=False,
                                                     getStored=True)
            dummy, QArr, rmsQArr = dataMan.get_specQ_byname(uniqueName)
            dummy, UArr, rmsUArr = dataMan.get_specU_byname(uniqueName)
        
            # Divide by the I model to create fractional polarisation spectra
            qArr = QArr / modIArr
            uArr = UArr / modIArr
            dqArr = qArr * np.sqrt( (rmsQArr/QArr)**2.0 +
                                    (rmsIArr/IArr)**2.0 )
            duArr = uArr * np.sqrt( (rmsUArr/UArr)**2.0 +
                                    (rmsIArr/IArr)**2.0 )

            # Query the parameters of the peak in the FDF
            pDict = dataMan.get_FDF_peak_params_byname(uniqueName)
            p = [float(x) for x in pDict["coeffPolyIspec"].split(",")]
            ampPeakIfreq0_Jybm = poly5(p)(pDict["freq0_Hz"]/1e9)
            fracPol = pDict["ampPeakPIfit_Jybm"]/ampPeakIfreq0_Jybm
            phiArr, ccFDF = dataMan.get_ccFDF_byname(uniqueName)
        
            # Run the complexity measurement function
            C1, C2, C3 = assess_complexity(freqArr_Hz, qArr, uArr, dqArr, duArr,
                                           fracPol, pDict["polAngle0Fit_deg"],
                                           pDict["phiPeakPIfit_rm2"])
        
            # Write the measurements to the record array
            compRec[i]["uniqueName"] = uniqueName
            compRec[i]["complexM1"] = C1
            compRec[i]["complexM2"] = C2
            compRec[i]["complexM3"] = C3
            compRec[i]["status"] = 1
        except Exception:
            compRec[i]['uniqueName'] = catRec[i]['uniqueName']
            compRec[i]['status'] = 0
            log_wr(LF, "Err: failed to process %s." % catRec[i]['uniqueName'])
            log_wr(LF, traceback.format_exc())
            continue

    return compRec



#-----------------------------------------------------------------------------#
def assess_complexity(freqArr_Hz, qArr, uArr, dqArr, duArr, fracPol, psi0_deg,
                      RM_radm2, phiArr=None, ccFDF=None):
    
    # Fractional polarised intensity
    pArr = np.sqrt(qArr**2.0 + uArr**2.0 )
    dpArr = np.sqrt(dqArr**2.0 + duArr**2.0 )
    
    # Create a RM-thin model from the best fit peak
    pModArr, qModArr, uModArr = create_pqu_spectra_RMthin(freqArr_Hz, fracPol,
                                                          psi0_deg, RM_radm2)

    # Subtract the a peak RM model to create a residual q & u
    qResidArr = qArr - qModArr
    uResidArr = uArr - uModArr
    pResidArr = np.sqrt(qResidArr**2.0 + uResidArr**2.0)

    # Complexity metric 1
    C1 = (np.sum(np.power((pResidArr), 2.0) / (np.power(dpArr, 2)))
          /(len(pModArr)-1))
    
    # Complexity metric 2
    C2 = (np.sum(np.power((qArr-qModArr), 2.0) + np.power((uArr-uModArr), 2.0))
          /(len(pModArr)-1) )

    # Complexity metric 3 on the clean-component PI FDF (second moment of CC)
    C3 = None
    if phiArr is not None and ccFDF is not None:
        print "MEEP"
        if not np.all(ccFDF==0):
            C3 = calc_mom2_FDF(ccFDF, phiArr)

    return C1, C2, C3
    
    
    
    
    
    
    
    
