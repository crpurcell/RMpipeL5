#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     module_RM_synthesis.py                                            #
#                                                                             #
# PURPOSE:  Function perform RM-synthesis on spectra in the PPC.              #
#                                                                             #
# MODIFIED: 28-Nov-2014 by C. Purcell                                         #
#                                                                             #
# TODO:     * Catch errors in non-fatal way                                   #
#           * If doOverwrite is not set, read and return the existing result  #
#                                                                             #
#=============================================================================#
import os
import sys
import math as m
import numpy as np

from util_RM import *
from util_PPC import log_fail
from util_PPC import log_wr
from util_PPC import fail_not_exists
from util_PPC import poly5


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
             ('phiCentre_rm2', 'f8')]
    rmsfRec = np.zeros(len(specRec), dtype=dType)
        
    # Loop through the catalogue entries
    log_wr(LF, '\nPerforming RM-synthesis on the catalogue entries ...')
    for i in range(len(specRec)):
        
        log_wr(LF, "\nProcessing entry %d: '%s'." %
               (i+1, specRec[i]['uniqueName']))
            
        # Read in the I-model, Q & U spectra
        inFile = dataInDir + '/' + specRec[i]['uniqueName'] + \
                 '_specImodel.dat'
        freqArr_Hz, specImodel = np.loadtxt(inFile, unpack=True)
        inFile = dataInDir + '/' + specRec[i]['uniqueName'] +  '_specQ.dat'
        freqArr_Hz, specQ = np.loadtxt(inFile, unpack=True)
        inFile = dataInDir + '/' + specRec[i]['uniqueName'] +  '_specU.dat'
        freqArr_Hz, specU = np.loadtxt(inFile, unpack=True)
        specPI = np.sqrt( np.power(specQ, 2.0) + np.power(specU, 2.0) )
        
        # Read in the noise spectrum
        inFile = dataInDir + '/' + specRec[i]['uniqueName'] + \
                 '_rmsSpecQUavg.dat'
        freqArr_Hz, rmsSpecQUAvg = np.loadtxt(inFile, unpack=True)

        # Calculate the weighting as 1/sigma^2 or 1
        if weightType=='Variance':
            weightArr = 1.0 / np.power(rmsSpecQUAvg, 2.0)
        else:
            weightArr = np.ones(rmsSpecQUAvg.shape)
        log_wr(LF,"> Weight type is '%s'." % weightType)
        
        # Form the Q/I & U/I spectra (divide by the model Stokes I)
        specQfrac = specQ / specImodel
        specUfrac = specU / specImodel
        specPIfrac = specPI / specImodel
        log_wr(LF, '> Q & U data divided by the Stokes I model.')
        
        # Calculate the lambda sampling from the frequency array
        lamArr_m = C / freqArr_Hz
        lamSqArr_m2 = np.power(lamArr_m, 2.0)
        
        # Perform RM-synthesis on the spectrum
        # Returned dirty FDF has amplitude expressed as polarised fraction
        tmp = do_rmsynth(specQfrac, specUfrac, lamSqArr_m2, phiArr, weightArr)
        dirtyFDF, [phiSampArr, RMSFArr], lam0Sq, fwhmRMSF = tmp
        log_wr(LF,'> RM-synthesis completed on Q/I and U/I spectra.')

        # Determine the Stokes I value at lam0Sq from the Stokes I polynomial
        # Multiply the dirty FDF by Ifreq0 to recover the PI in Jy
        polyCoeffs = specRec[i]['coeffPolyIspec'].split(',')
        polyCoeffs = [float(x) for x in polyCoeffs]
        freq0_Hz = C / m.sqrt(lam0Sq)
        Ifreq0_Jybm = poly5(polyCoeffs)(freq0_Hz / 1e9)
        dirtyFDF *= Ifreq0_Jybm
        log_wr(LF, '> Dirty FDF recast into input units.')

        # Save the RMSF to simple text files
        np.savetxt(outDataDir + '/' + specRec[i]['uniqueName'] +
                   '_RMSF.dat', zip(phiSampArr, RMSFArr.real, RMSFArr.imag))

        # Save the weight array to a simple text file
        np.savetxt(outDataDir + '/' + specRec[i]['uniqueName'] +
                   '_Weight.dat', zip(freqArr_Hz, weightArr))
        
        # Save the dirty FDF to a simple text file
        np.savetxt(outDataDir + '/' + specRec[i]['uniqueName'] +
                   '_dirtyFDF.dat', zip(phiArr, dirtyFDF.real, dirtyFDF.imag))
        log_wr(LF, '> Dirty FDF saved to ASCII file.')

        # Write flags and metadata to the catalogue file
        rmsfRec[i]['uniqueName'] = specRec[i]['uniqueName']
        rmsfRec[i]['weightType'] = weightType
        rmsfRec[i]['lam0Sq_m2'] = lam0Sq
        rmsfRec[i]['freq0_Hz'] = freq0_Hz
        rmsfRec[i]['nPhiChan'] = int(len(phiArr))
        rmsfRec[i]['fwhmRMSF'] = fwhmRMSF
        rmsfRec[i]['deltaPhiChan_rm2'] = np.nanmin(np.diff(phiArr))
        rmsfRec[i]['phiCentre_rm2'] = np.median(phiArr)

    return rmsfRec
