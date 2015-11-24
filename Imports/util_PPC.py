#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     util_PPC.py                                                       #
#                                                                             #
# PURPOSE:  Helper functions for the POSSUM pipeline.                         #
#                                                                             #
# REQUIRED: Requires numpy and astropy.                                       #
#                                                                             #
# MODIFIED: 20-November-2015 by C. Purcell                                    #
#                                                                             #
# CONTENTS:                                                                   #
#                                                                             #
#  PipelineInputs       ... class to contain the pipeline inputs and file     #
#  DataManager          ... class to interface with the database and data     #
#  PlotParms            ... class containing query and plotting parameters    #
#  cleanup_str_input    ... condense multiple newlines and spaces in a string #
#  config_read          ... read a key=value format text file                 #
#  deg2dms              ... convert decimal degrees to dms string             #
#  fail_not_exists      ... check for dir/file existence, exit if missing     #
#  set_statusfile       ... set text in an ASCII status file                  #
#  read_statusfile      ... read text in an ASCII status file                 #
#  read_dictfile        ... save a dictionary to a text file                  #
#  write_dictfile       ... write a dictionary to a text file                 #
#  load_vector_fail     ... load a single-column vector from a text file      #
#  cat_to_recarray      ... convert ASCII catalogue file to a record array    #
#  extract_spec_noise   ... extract source and noise spectra from FITS file(s)#
#  calc_sumbox_norm     ... calculate the beam normalisation for a box MxM    #
#  calc_mom2_FDF        ... calculate the 2nd moment of the CC                #
#  poly5                ... function to evaluate a 5th order polynomial       #
#  log_fail             ... record a fatal error in the log file and exit     #
#  log_wr               ... record a message in the log file                  #
#  nanmedian            ... np.median ignoring NaNs                           #
#  MAD                  ... calculate the madfm                               #
#  calc_stats           ... calculate the statistics of an array              #
#  sort_nicely          ... sort a list in the order a human would            #
#  twodgaussian         ... return an array containing a 2D Gaussian          #
#  create_pqu_spectra_RMthin ... return fractional spectra for a thin source  #
#  create_IQU_spectra_RMthin ... return IQU spectra for a thin source         #
#  create_pqu_resid_RMthin ... return fractional spectra - a thin component   #
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

import os
import sys
import copy
import re
import time
import math as m
import numpy as np
import numpy.ma as ma
import ConfigParser
import sqlite3
import csv
import json

import astropy.io.fits as pf

from util_DB import select_into_arr
from util_DB import get_tables_description

C = 2.997924538e8 # Speed of light [m/s]


#-----------------------------------------------------------------------------#
class PipelineInputs:
    """
    Class to hold the inputs and defaults for the POSSUM pipeline.
    Methods assume unique keys across all configuration sections. There are
    two levels of parameters: 1) input parameters linked to a file and stored
    in a ConfigParser object, 2) derived parameters stored in a dictionary.
    Derived parameters pertain to the frequency and Faraday depth sampling and
    are set to sensible defaults tuned to the data.
    """

    def __init__(self, configFile="", calcParms=False, resetPhiSamp=False):
        """
        Initialise the config from an existing file or defaults defined here.
        """
        self.configFile = configFile
        self.sessionDir = os.path.split(configFile)[0]
        self.config = ConfigParser.ConfigParser()
        self.config.optionxform = str
        self.derivedParmDict = {}
        if os.path.exists(self.configFile):
            self.read_file()
        else:
            self._setup_defaults()
        if calcParms:
            self.calculate_derived_parms(resetPhiSamp=resetPhiSamp)
        
    def _setup_defaults(self):
        """
        Hard-coded initial defaults for the RM pipeline.
        """
        self.config.add_section("Dataset")
        self.config.set("Dataset", "dataPath", ".")
        self.config.set("Dataset", "sumBoxPix", "3")
        self.config.set("Dataset", "polyOrd", "2")
        
        self.config.add_section("RMsynthesis")
        self.config.set("RMsynthesis", "dPhi_radm2", "10.0")
        self.config.set("RMsynthesis", "phiMax_radm2", "0.0")
        self.config.set("RMsynthesis", "oversampling", "10.0")
        self.config.set("RMsynthesis", "weightType", "Variance")
        
        self.config.add_section("RMclean")
        self.config.set("RMclean", "cleanAlgorithm", "Heald")
        self.config.set("RMclean", "gain", "0.1")
        self.config.set("RMclean", "cleanCutoff_sigma", "5.0")
        self.config.set("RMclean", "maxCleanIter", "1000")
        
        self.config.add_section("Thresholds")
        self.config.set("Thresholds", "thresholdPolBias_sigma", "5.0")
        self.config.set("Thresholds", "thresholdSignalPI_sigma", "8.0")
        self.config.set("Thresholds", "thresholdDoClean_sigma", "3.0")
        
        self.config.add_section("Complexity")
        self.config.set("Complexity", "threshC1", "3.0")
        self.config.set("Complexity", "threshC2", "3.0")
        self.config.set("Complexity", "threshC3", "3.0")
        
    def read_file(self):
        """
        Read the inputs from an existing configuration file.
        """
        self.config.read(self.configFile)
    
    def write_file(self):
        """
        Save a configuration file (overwrite if it exists).
        """
        if os.path.exists(self.configFile):
            os.remove(self.configFile)
        CFG = open(self.configFile, "w")
        self.config.write(CFG)
        CFG.close()

    def get_flat_dict(self, includeDerived=False):
        """
        Return a dictionary of key=val assuming unique keys across all
        sections.
        """
        pDict = {}
        for section in self.config.sections():
            for k, v in self.config.items(section):
                pDict[k] = v
        if includeDerived:
            for k, v in self.derivedParmDict.iteritems():
                pDict[k] = v            
        return pDict
    
    def set_flat_dict(self, pDict):
        """
        Set values in the ConfigParser object using a dictionary. Values
        will only be updated if the keywords already exist in the object.
        A dictionary of non-matching keys are returned. To add new
        key/value pairs to the ConfigParser, use the 'config.set' method.
        """
        pDict1 = copy.copy(pDict)
        sectionLst = self.config.sections()
        for section in sectionLst:
            for k in self.config.options(section):
                if k in pDict1:
                    self.config.set(section, k, str(pDict1[k]))
                    pDict1.pop(k)
        return pDict1

    def calculate_derived_parms(self, resetPhiSamp=False):
        """
        Calculate the secondary parameters needed to run the RM-pipeline.
        """
        
        # Read the dataType.txt file
        dataPath = self.config.get("Dataset", "dataPath")
        dataTypeFile = dataPath + '/dataType.txt'
        if not os.path.exists(dataTypeFile):
            dataTypeFile = self.sessionDir + '/dataType.txt'
        dataType = read_statusfile(dataTypeFile)

        # Frequency and wavelength sampling
        # These are set by the dataset and cannot be altered by the user
        inFreqFile = dataPath + '/freqs_Hz.txt'
        if not os.path.exists(inFreqFile):
            inFreqFile =  self.sessionDir + '/freqs_Hz.txt'
        freqArr_Hz = load_vector_fail(inFreqFile, "float32", False)
        dFreq_Hz = np.nanmin(np.abs(np.diff(freqArr_Hz)))
        lambdaArr_m = C / freqArr_Hz
        lambdaSqArr_m2 = np.power(lambdaArr_m, 2.0)
        lambdaSqRange_m2 = ( np.nanmax(lambdaSqArr_m2) -
                             np.nanmin(lambdaSqArr_m2) )        
        dLambdaSqMin_m2 = np.nanmin(np.abs(np.diff(lambdaSqArr_m2)))
        dLambdaSqMax_m2 = np.nanmax(np.abs(np.diff(lambdaSqArr_m2)))

        # Default Faraday depth limits
        # The Faraday sampling can be forced by existing config options
        fwhmRMSF_radm2 = 2.0 * m.sqrt(3.0) / lambdaSqRange_m2
        if resetPhiSamp:
            oversampling = float(self.config.get("RMsynthesis",
                                                 "oversampling"))
            dPhi_radm2 = fwhmRMSF_radm2 / oversampling
            phiMax_radm2 = m.sqrt(3.0) / dLambdaSqMax_m2
        else:
            dPhi_radm2 = float(self.config.get("RMsynthesis", "dPhi_radm2"))
            phiMax_radm2 = float(self.config.get("RMsynthesis", 
                                                 "phiMax_radm2"))

        # Force the minimum phiMax
        phiMax_radm2 = max(phiMax_radm2, 600.0)
            
        # Faraday depth sampling. Zero always centred on middle channel
        nChanRM = round(abs((phiMax_radm2 - 0.0) / dPhi_radm2)) * 2.0 + 1.0
        startPhi_radm2 = - (nChanRM-1.0) * dPhi_radm2 / 2.0
        stopPhi_radm2 = + (nChanRM-1.0) * dPhi_radm2 / 2.0
        phiArr_radm2 = np.arange(startPhi_radm2, stopPhi_radm2 + dPhi_radm2,
                                 dPhi_radm2)

        # Save the derived parameters
        self.derivedParmDict["dataType"] = dataType
        self.derivedParmDict["freqArr_Hz"] = freqArr_Hz
        self.derivedParmDict["dFreq_Hz"] = dFreq_Hz
        self.derivedParmDict["nChanFreq"] = len(freqArr_Hz)
        self.derivedParmDict["lambdaSqArr_m2"] = lambdaSqArr_m2
        self.derivedParmDict["lambdaSqRange_m2"] = lambdaSqRange_m2
        self.derivedParmDict["dLambdaSqMin_m2"] = dLambdaSqMin_m2
        self.derivedParmDict["dLambdaSqMax_m2"] = dLambdaSqMax_m2
        self.derivedParmDict["fwhmRMSF_radm2"] = fwhmRMSF_radm2        
        self.config.set("RMsynthesis", "dPhi_radm2", str(dPhi_radm2))
        self.config.set("RMsynthesis", "phiMax_radm2", str(phiMax_radm2))
        self.derivedParmDict["phiArr_radm2"] = phiArr_radm2
        self.derivedParmDict["nChanRM"] = int(nChanRM)
        self.derivedParmDict["startPhi_radm2"] = startPhi_radm2
        self.derivedParmDict["stopPhi_radm2"] = stopPhi_radm2
        self.derivedParmDict["phiCentre_radm2"] = 0.0
        
    def inparm_verify(self):
        """
        Verify that required input parameters are all present.
        """
        requiredKeyLst = ["dataPath",
                          "sumBoxPix",
                          "polyOrd",
                          "dPhi_radm2",
                          "phiMax_radm2",
                          "oversampling",
                          "weightType",
                          "cleanAlgorithm",
                          "gain",
                          "cleanCutoff_sigma",
                          "maxCleanIter",
                          "thresholdPolBias_sigma",
                          "thresholdSignalPI_sigma",
                          "thresholdDoClean_sigma",
                          "threshC1",
                          "threshC2",
                          "threshC3"]
        pDict = self.get_flat_dict(includeDerived=False)
        keys = pDict.keys()
        missingLst = [x for x in requiredKeyLst if x not in keys]

        if len(missingLst) > 0:
            return missingLst
        else:
            return []
        

#-----------------------------------------------------------------------------#
class DataManager:
    """
    Class to interface with the database, datasets and results.

    Public methods:
        indx2name                  ... translate unique name to table row
        name2indx                  ... translate table row to unique name
        get_specI_byname           ... return the I spectrum array
        get_specI_byindx
        get_modI_byname            ... return the I model array
        get_modI_byindx
        get_specQ_byname           ... return the Q spectrum array
        get_specQ_byindx
        get_specU_byname           ... return the U spectrum array
        get_specU_byindx
        get_freqweight_byname      ... return the frequency and weight arrays
        get_freqweight_byindx
        get_RMSF_byname            ... return the complex RMSF array
        get_RMSF_byindx
        get_dirtyFDF_byname        ... return the dirty FDF array (complex)
        get_dirtyFDF_byindx
        create_RMSFfits_byname     ... create FITS file to store RM-synth
        create_RMSFfits_byindx
        get_cleanFDF_byname        ... return the clean FDF array (complex)
        get_cleanFDF_byindx
        get_ccFDF_byname           ... return the clean-component PI FDF array
        get_ccFDF_byindx
        put_cleanFDF_byname        ... insert a clean FDF into FITS file
        put_cleanFDF_byindx
        get_Imodel_coeffs_byname   ... return the polynomial coefficient
        get_Imodel_coeffs_byindx
        get_RMSF_params_byname     ... return the parameters of the RMSF
        get_RMSF_params_byindx
        get_FDF_peak_params_byname ... query the DB for params of peak
        get_FDF_peak_params_byindx
        get_thin_qumodel_byname    ... return the thin fractional QU model
        get_thin_qumodel_byindx
        get_stampI_byname          ... return the I postage stamp data & header
        get_stampI_byindx
        get_stampP_byname          ... return the P postage stamp data & header
        get_stampP_byindx
        get_database_schema        ... return a dict containing the DB schema
        query_database             ... run a SQL query on the DB
        close                      ... close a conection to the DB
        export_table               ... save a table to disk
    """

    def __init__(self, sessionPath, calcParms=True):
        self.sessionPath = sessionPath
        self.dbFile =  sessionPath + "/session.sqlite"
        self.pipeInpObj = None
        self.pDict = None
        self.summaryRec = None
        self.tempRec = None
        
        # Load the session inputs
        self._load_pipeInputs(calcParms=calcParms)
        
        # Connect to the database and load a summary table
        conn = sqlite3.connect(self.dbFile)
        cursor = conn.cursor()
        self._load_summaryTab(cursor)
        cursor.close()
        conn.close()

    def _load_pipeInputs(self, calcParms):

        # Create a PipelineInput class, in the process reading in the 
        # pipeline input file and calculating the derived parameters
        inputFile = self.sessionPath + "/inputs.config"
        self.pipeInpObj = PipelineInputs(inputFile, calcParms,
                                         resetPhiSamp=False)
        missingLst = self.pipeInpObj.inparm_verify()
        if len(missingLst)>0:
            print "Err: Required input parameters missing. - %s" % missingLst
            sys.exit(0)
        self.pDict = self.pipeInpObj.get_flat_dict(includeDerived=calcParms)
        
    def _load_summaryTab(self, cursor):
        sql = """
        SELECT
        sourceCat.uniqueName,
        spectraParms.fluxMedI_Jybm as fluxI_mJy,
        dirtyFDFparms.ampPeakPIfit_Jybm as peakPI_mJybm,
        dirtyFDFparms.phiPeakPIfit_rm2 as RM_radm2,
        dirtyFDFparms.snrPIfit as SNR,
        polAngleFit_deg as PA_deg
        FROM sourceCat
        LEFT JOIN spectraParms
        ON sourceCat.uniqueName = spectraParms.uniqueName
        LEFT JOIN dirtyFDFparms
        ON sourceCat.uniqueName = dirtyFDFparms.uniqueName
        """
        self.summaryRec = select_into_arr(cursor, sql)

    def indx2name(self, indx):
        return self.summaryRec["uniqueName"][indx]
    
    def name2indx(self, name):
        return np.argwhere(self.summaryRec["uniqueName"]==name)[0][0]

    def get_specI_byname(self, uniqueName):
        specDir = self.sessionPath + '/OUT'
        dataFile = specDir +  '/' + uniqueName +  '_specI.fits'
        HDULst = pf.open(dataFile, "readonly", memmap=True)
        freqArr_Hz = HDULst[2].data["freq"]
        IArr_Jy = HDULst[2].data["src"]
        rmsIArr_Jy = HDULst[2].data["rms"]
        HDULst.close()
        return freqArr_Hz, IArr_Jy, rmsIArr_Jy
    
    def get_specI_byindx(self, indx):
        uniqueName = self.indx2name(indx)
        return self.get_specI_byname(uniqueName)

    def get_modI_byname(self, uniqueName, oversample=False, getStored=False):
        specDir = self.sessionPath + '/OUT'
        dataFile = specDir +  '/' + uniqueName +  '_specI.fits'
        HDULst = pf.open(dataFile, "readonly", memmap=True)
        freqArr_Hz = HDULst[2].data["freq"]
        if oversample:
            freqArr_Hz = np.linspace(freqArr_Hz[0], freqArr_Hz[-1], 10000)
            p = self.get_Imodel_coeffs_byname(uniqueName)
            modIArr_Jy = poly5(p)(freqArr_Hz/1e9)
        else:
            modIArr_Jy = HDULst[2].data["model"]
        HDULst.close()
        return freqArr_Hz, modIArr_Jy
    
    def get_modI_byindx(self, indx, oversample=False, getStored=False):
        uniqueName = self.indx2name(indx)
        return self.get_modI_byname(uniqueName, oversample, getStored)
    
    def get_specQ_byname(self, uniqueName):
        specDir = self.sessionPath + '/OUT'
        dataFile = specDir +  '/' + uniqueName +  '_specQ.fits'
        HDULst = pf.open(dataFile, "readonly", memmap=True)
        freqArr_Hz = HDULst[2].data["freq"]
        QArr_Jy = HDULst[2].data["src"]
        rmsQArr_Jy = HDULst[2].data["rms"]
        HDULst.close()
        return freqArr_Hz, QArr_Jy, rmsQArr_Jy
        
    def get_specQ_byindx(self, indx):
        uniqueName = self.indx2name(indx)
        return self.get_specQ_byname(uniqueName)

    def get_specU_byname(self, uniqueName):
        specDir = self.sessionPath + '/OUT'
        dataFile = specDir +  '/' + uniqueName +  '_specU.fits'
        HDULst = pf.open(dataFile, "readonly", memmap=True)
        freqArr_Hz = HDULst[2].data["freq"]
        UArr_Jy = HDULst[2].data["src"]
        rmsUArr_Jy = HDULst[2].data["rms"]
        HDULst.close()
        return freqArr_Hz, UArr_Jy, rmsUArr_Jy
        
    def get_specU_byindx(self, indx):
        uniqueName = self.indx2name(indx)
        return self.get_specU_byname(uniqueName)
    
    def get_freqweight_byname(self, uniqueName):
        specDir = self.sessionPath + '/OUT'
        dataFile = specDir +  '/' + uniqueName +  '_RMSynth.fits'
        HDULst = pf.open(dataFile, "readonly", memmap=True)
        freqArr_Hz = HDULst[1].data["freq"]
        weightArr = HDULst[1].data["weight"]
        HDULst.close()
        return freqArr_Hz, weightArr

    def get_freqweight_byindx(self, indx):
        uniqueName = self.indx2name(indx)
        return self.get_freqweight_byname(uniqueName)
    
    def get_RMSF_byname(self, uniqueName):
        specDir = self.sessionPath + '/OUT'
        dataFile = specDir +  '/' + uniqueName +  '_RMSynth.fits'
        HDULst = pf.open(dataFile, "readonly", memmap=True)        
        phiArr = HDULst[2].data["phiSamp"]
        RMSFreal = HDULst[2].data["RMSFreal"]
        RMSFimag = HDULst[2].data["RMSFimag"]
        HDULst.close()
        RMSFArr =  (RMSFreal + 1j * RMSFimag)
        return phiArr, RMSFArr
        
    def get_RMSF_byindx(self, indx):
        uniqueName = self.indx2name(indx)
        return self.get_RMSF_byname(uniqueName)

    def get_dirtyFDF_byname(self, uniqueName):
        specDir = self.sessionPath + '/OUT'
        dataFile = specDir +  '/' + uniqueName +  '_RMSynth.fits'
        HDULst = pf.open(dataFile, "readonly", memmap=True)        
        phiArr = HDULst[3].data["phi"]
        FDFreal_Jy = HDULst[3].data["FDFreal"]
        FDFimag_Jy = HDULst[3].data["FDFimag"]
        HDULst.close()
        FDFArr_Jy =  (FDFreal_Jy + 1j * FDFimag_Jy)
        return phiArr, FDFArr_Jy

    def get_dirtyFDF_byindx(self, indx):
        uniqueName = self.indx2name(indx)
        return self.get_dirtyFDF_byname(uniqueName)

    def create_RMSFfits_byname(self, uniqueName, freqArr_Hz, weightArr,
                               phiSampArr, RMSFArr, phiArr, dirtyFDF):
        specDir = self.sessionPath + '/OUT'
        dataFile = specDir +  '/' + uniqueName +  '_RMSynth.fits'
        hdu0 = pf.PrimaryHDU(header=pf.Header())

        lamArr_m = C / freqArr_Hz
        lamSqArr_m2 = np.power(lamArr_m, 2.0)        
        
        col1 = pf.Column(name="freq", format="f4", array=freqArr_Hz)
        col2 = pf.Column(name="lamsq", format="f4", array=lamSqArr_m2)
        col3 = pf.Column(name="weight", format="f4", array=weightArr)
        hdu1 = pf.new_table([col1, col2, col3])
        #hdu1 = pf.BinTableHDU.from_columns([col1, col2, col3])
            
        col1 = pf.Column(name="phiSamp", format="f4", array=phiSampArr)
        col2 = pf.Column(name="RMSFreal", format="f4", array=RMSFArr.real)
        col3 = pf.Column(name="RMSFimag", format="f4", array=RMSFArr.imag)
        hdu2 = pf.new_table([col1, col2, col3])
        #hdu2 = pf.BinTableHDU.from_columns([col1, col2, col3])
            
        col1 = pf.Column(name="phi", format="f4", array=phiArr)
        col2 = pf.Column(name="FDFreal", format="f4", array=dirtyFDF.real)
        col3 = pf.Column(name="FDFimag", format="f4", array=dirtyFDF.imag)
        hdu3 = pf.new_table([col1, col2, col3])
        #hdu3 = pf.BinTableHDU.from_columns([col1, col2, col3])

        col1 = pf.Column(name="CC", format="f4",
                         array=np.zeros_like(dirtyFDF.real))
        col2 = pf.Column(name="FDFreal", format="f4",
                         array=np.zeros_like(dirtyFDF.real))
        col3 = pf.Column(name="FDFimag", format="f4",
                         array=np.zeros_like(dirtyFDF.real))
        hdu4 = pf.new_table([col1, col2, col3])
        #hdu4 = pf.BinTableHDU.from_columns([col1, col2, col3])
            
        hduLst = pf.HDUList([hdu0, hdu1, hdu2, hdu3, hdu4])
        hduLst.writeto(dataFile, output_verify="fix", clobber=True)
        hduLst.close()

    def create_RMSFfits_byindx(self, indx, **kwargs):
        uniqueName = self.indx2name(indx)
        return self.create_RMSFfits_byname(self, uniqueName, **kwargs)
        
    def get_cleanFDF_byname(self, uniqueName):
        specDir = self.sessionPath + '/OUT'
        dataFile = specDir +  '/' + uniqueName +  '_RMSynth.fits'
        HDULst = pf.open(dataFile, "readonly", memmap=True)        
        phiArr = HDULst[3].data["phi"]
        FDFreal_Jy = HDULst[4].data["FDFreal"]
        FDFimag_Jy = HDULst[4].data["FDFimag"]
        HDULst.close()
        FDFArr_Jy =  (FDFreal_Jy + 1j * FDFimag_Jy)
        return phiArr, FDFArr_Jy
    
    def get_cleanFDF_byindx(self, indx):
        uniqueName = self.indx2name(indx)
        return self.get_cleanFDF_byname(uniqueName)
    
    def get_ccFDF_byname(self, uniqueName):
        specDir = self.sessionPath + '/OUT'
        dataFile = specDir +  '/' + uniqueName +  '_RMSynth.fits'
        HDULst = pf.open(dataFile, "readonly", memmap=True)        
        phiArr = HDULst[3].data["phi"]
        ccFDF_Jy = HDULst[4].data["CC"]        
        return phiArr, ccFDF_Jy

    def get_ccFDF_byindx(self, indx):
        uniqueName = self.indx2name(indx)
        return self.get_ccFDF_byname(uniqueName)

    def put_cleanFDF_byname(self, uniqueName, CC=None, cleanFDF=None):
        specDir = self.sessionPath + '/OUT'
        dataFile = specDir +  '/' + uniqueName +  '_RMSynth.fits'
        HDULst = pf.open(dataFile, "update", memmap=True)
        if CC is not None:
            HDULst[4].data["CC"] = np.abs(CC)
        if cleanFDF is not None:
            HDULst[4].data["FDFreal"] = cleanFDF.real
            HDULst[4].data["FDFimag"] = cleanFDF.imag
        HDULst.close()
        
    def put_cleanFDF_byindx(self, uniqueName, **kwargs):
        uniqueName = self.indx2name(indx)
        put_cleanFDF_byname(self, uniqueName, **kwargs)
        
    def get_Imodel_coeffs_byname(self, uniqueName):
        # Connect to the database and fetch the parameters
        conn = sqlite3.connect(self.dbFile)
        cursor = conn.cursor()
        sql = """
        SELECT coeffPolyIspec
        FROM spectraParms
        WHERE uniqueName='%s'
        """ % (uniqueName)
        resultArr = self.query_database(sql, buffer=False)
        cursor.close()
        coeffLst = []
        if len(resultArr)>0:
            coeffLst = [float(x) for x in resultArr[0][0].split(",")]
        return coeffLst

    def get_Imodel_coeffs_byindx(self, indx):
        uniqueName = self.indx2name(indx)
        return self.get_Imodel_params_byname(uniqueName)
        
    def get_RMSF_params_byname(self, uniqueName):
        
        # Connect to the database and fetch the parameters
        conn = sqlite3.connect(self.dbFile)
        cursor = conn.cursor()
        sql = """
        SELECT
        lam0Sq_m2,
        freq0_Hz,
        nPhiChan,
        deltaPhiChan_rm2,
        phiCentre_rm2,
        weightType,
        fwhmRMSF
        FROM dirtyFDFparms
        WHERE uniqueName='%s'
        """ % (uniqueName)
        resultArr = self.query_database(sql, buffer=False)
        cursor.close()
        conn.close()
        pDict = {}
        if len(resultArr)>0:
            for key,val in zip(resultArr.dtype.names, resultArr[0]):
                pDict[key] = val
        return pDict        

    def get_RMSF_params_byindx(self, indx):
        uniqueName = self.indx2name(indx)
        return self.get_RMSF_params_byname(uniqueName)
    
    def get_FDF_peak_params_byname(self, uniqueName, doClean=False):

        if doClean:
            FDFtable = "cleanFDFparms"            
        else:
            FDFtable = "dirtyFDFparms"
    
        # Connect to the database and fetch the parameters of the peak
        conn = sqlite3.connect(self.dbFile)
        cursor = conn.cursor()
        sql = "SELECT "
        sql += "dirtyFDFparms.lam0Sq_m2, "
        sql += "dirtyFDFparms.freq0_Hz, "
        sql += "%s.phiPeakPIfit_rm2, " % (FDFtable)
        sql += "%s.dPhiPeakPIfit_rm2, " % (FDFtable)
        sql += "%s.ampPeakPIfit_Jybm, " % (FDFtable)
        sql += "%s.dAmpPeakPIfit_Jybm, " % (FDFtable)
        sql += "%s.polAngleFit_deg, " % (FDFtable)
        sql += "%s.dPolAngleFit_deg, " % (FDFtable)
        sql += "%s.polAngle0Fit_deg, " % (FDFtable)
        sql += "%s.dPolAngle0Fit_deg, " % (FDFtable)
        sql += "%s.detectionF, " % (FDFtable)
        sql += "cleanFDFparms.cleanCutoff_sigma, "
        sql += "cleanFDFparms.cleanCutoff_Jybm, "        
        sql += "spectraParms.coeffPolyIspec "
        sql += "FROM spectraParms INNER JOIN dirtyFDFparms "
        sql += "ON dirtyFDFparms.uniqueName = spectraParms.uniqueName "
        sql += "LEFT JOIN cleanFDFparms "
        sql += "ON cleanFDFparms.uniqueName = spectraParms.uniqueName "
        sql += "WHERE spectraParms.uniqueName='%s'" % (uniqueName)
        resultArr = self.query_database(sql, buffer=False)
        cursor.close()
        conn.close()
        pDict = {}
        if len(resultArr)>0:
            for key,val in zip(resultArr.dtype.names, resultArr[0]):
                pDict[key] = val
        return pDict

    def get_FDF_peak_params_byindx(self, indx, doClean=False):
        uniqueName = self.indx2name(indx)
        return self.get_FDF_peak_params_byname(uniqueName, doClean)
        
    def get_thin_qumodel_byname(self, uniqueName, oversample=False, 
                                doClean=False):
        
        # Determine the frequency sampling
        freqArr_Hz, modIArr_Jy = self.get_modI_byname(uniqueName)
        if oversample:
            freqArr_Hz = np.linspace(freqArr_Hz[0], freqArr_Hz[-1], 10000)
            
        # Get the parameters of the peak in the FDF
        pDict = self.get_FDF_peak_params_byname(uniqueName, doClean)
        p = [float(x) for x in pDict["coeffPolyIspec"].split(",")]

        # Calculate the fractional polarisation at lambda0/freq0
        ampPeakIfreq0_Jybm = poly5(p)(pDict["freq0_Hz"]/1e9)
        fracPol = pDict["ampPeakPIfit_Jybm"]/ampPeakIfreq0_Jybm

        # Create a model spectrum
        pArr, qArr, uArr = create_pqu_spectra_RMthin(freqArr_Hz, fracPol,
                                                   pDict["polAngle0Fit_deg"], 
                                                   pDict["phiPeakPIfit_rm2"])
        return freqArr_Hz, qArr, uArr

    def get_thin_qumodel_byindx(self, indx, oversample=False, 
                                doClean=False):
        uniqueName = self.indx2name(indx)
        return self.get_thin_qumodel_byname(uniqueName, oversample, doClean)

    def get_stampI_byname(self, uniqueName):
        dataDir = self.sessionPath + '/OUT'
        dataFile = dataDir +  '/' + uniqueName +  '_specI.fits'
        head = pf.getheader(dataFile)
        data = pf.getdata(dataFile)[0]
        return data, head

    def get_stampI_byindx(self, indx):
        uniqueName = self.indx2name(indx)
        return self.get_stampI_byname(uniqueName)
    
    def get_stampP_byname(self, uniqueName):
        dataDir = self.sessionPath + '/OUT'
        dataQFile = dataDir +  '/' + uniqueName +  '_specQ.fits'
        head = pf.getheader(dataQFile)
        dataQ = pf.getdata(dataQFile)[0]
        dataUFile = dataDir +  '/' + uniqueName +  '_specU.fits'
        dataU = pf.getdata(dataUFile)[0]
        data = np.sqrt(dataQ**2.0+dataU**2.0)        
        return data, head
    
    def get_stampP_byindx(self, indx):
        uniqueName = self.indx2name(indx)
        return self.get_stampP_byname(uniqueName)
    
    def get_database_schema(self):
        conn = sqlite3.connect(self.dbFile)
        cursor = conn.cursor()
        descDict = get_tables_description(cursor)
        cursor.close()
        conn.close()
        return descDict

    def query_database(self, sql, buffer=True, doCommit=False):
        conn = sqlite3.connect(self.dbFile)
        cursor = conn.cursor()
        resultArr = None
        resultArr = select_into_arr(cursor, sql)
        if buffer:
            self.tempRec = resultArr
        if doCommit:
            conn.commit()
        cursor.close()
        conn.close()
        return resultArr    

    def close(self):
        cursor.close()
        conn.close()

    def export_table(self, tableName, format="CSV", outDir=".", 
                     saveBuffer=False):
        if saveBuffer is False:
            sql = "SELECT * FROM %s" % tableName
            self.query_database(sql, buffer=True)
        outFileName = outDir + "/" +tableName + ".csv"
        with open(outFileName, "wb") as FH:
            writer = csv.writer(FH)
            colNames = self.tempRec.dtype.names
            writer.writerow(colNames)
            writer.writerows(self.tempRec)
        FH.close()
        print "Table '%s' written to file:\n'%s'" % (tableName, outFileName)


#-----------------------------------------------------------------------------#
class PlotParms:
    """
    Class to store plotting parameters and queries.
    """
    
    def __init__(self, queryFile=None):
        self.configDict = {}
        self.queryLst = []
        self.queryLabLst = []
        loadDefaults = False
        if queryFile is not None:
            try:
                self.configDict = parse_queryfile(queryFile)
                self.queryLst = self.configDict.pop('QUERY')
                self.queryLabLst = self.configDict.pop('QLABEL')
            except Exception:
                loadDefaults = True
        if loadDefaults:
            self._setup_defaults()
                
    def _setup_defaults(self):
        self.configDict["DBFILE"] = ""
        self.configDict["TYPE"] = "Histogram"
        self.configDict["DOLOGX"] = "0"
        self.configDict["DOLOGY"] = "0"
        self.configDict["NBINS"] = "10"
        self.configDict["ZPOWER"] = "1.0"
        self.configDict["XDATAMIN"] = ""
        self.configDict["XDATAMAX"] = ""
        self.configDict["YDATAMIN"] = ""
        self.configDict["YDATAMAX"] = ""
        self.configDict["ZDATAMIN"] = ""
        self.configDict["ZDATAMAX"] = ""
        self.configDict["XLABEL"] = ""
        self.configDict["YLABEL"] = ""
        self.configDict["ZLABEL"] = ""
        self.configDict["TITLE"] = ""
        self.queryLst = []
        self.queryLabLst = ["Query 1", "Query 2", "Query 3", "Query 4"]

        
#-----------------------------------------------------------------------------#
def parse_queryfile(filename):
    """Parse a file containing the SQL queries and key=value pairs."""
    
    configDict = dict()              # Dictionary to hold keyword-value pairs
    queryLst=[]                      # List to hold queries
    queryLabLst = []                 # List to hold legend labels
    CONFIGFILE = open(filename, "r")

    # Compile a few useful regular expressions
    spaces = re.compile('\s+')
    comma_and_spaces = re.compile(',\s+')
    comma_or_space = re.compile('[\s|,]')
    brackets = re.compile('[\[|\]\(|\)|\{|\}]')
    comment = re.compile('#.*')
    quotes = re.compile('\'[^\']*\'')
    key_val = re.compile('^.+=.+')

    # Read in the input file, line by line
    for line in CONFIGFILE:
        line = line.rstrip("\n\r")

        # Filter for comments and blank lines
        if not comment.match(line) and key_val.match(line):

            # Weed out internal comments & split on 1st '='
            line = comment.sub('',line)
            (keyword, value) = line.split('=',1)

            # If the line contains a value            
            keyword = keyword.strip()          # kill external whitespace
            keyword = spaces.sub('', keyword)  # kill internal whitespaces
            value = value.strip()              # kill external whitespace
            value = spaces.sub(' ', value)     # shrink internal whitespace
            value = comma_and_spaces.sub(',', value) # kill ambiguous spaces

            # Separate out the queries 
            if value:
                if (keyword=='QUERY'):
                    queryLst.append(value)
                elif (keyword=='QLABEL'):
                    queryLabLst.append(value)
                else:                
                    configDict[keyword] = value
            configDict['QUERY'] = queryLst
            configDict['QLABEL'] = queryLabLst

    CONFIGFILE.close()
    
    return configDict
        
        


#-----------------------------------------------------------------------------#
def cleanup_str_input(textBlock):
    
    # Compile a few useful regular expressions
    spaces = re.compile(r"[^\S\r\n]+")
    newlines = re.compile(r"\n+")
    rets = re.compile(r"\r+")

    # Strip multiple spaces etc
    textBlock = textBlock.strip()
    textBlock = rets.sub('\n', textBlock)
    textBlock = newlines.sub('\n', textBlock)
    textBlock = spaces.sub(' ', textBlock)

    return textBlock

    
#-----------------------------------------------------------------------------#
def config_read(filename, delim='=', doValueSplit=True):
    """
    Read a configuration file and output a 'KEY=VALUE' dictionary.
    """

    configTable = {}
    CONFIGFILE = open(filename, "r")
    
    # Compile a few useful regular expressions
    spaces = re.compile('\s+')
    commaAndSpaces = re.compile(',\s+')
    commaOrSpace = re.compile('[\s|,]')
    brackets = re.compile('[\[|\]\(|\)|\{|\}]')
    comment = re.compile('#.*')
    quotes = re.compile('\'[^\']*\'')
    keyVal = re.compile('^.+' + delim + '.+')

    # Read in the input file, line by line
    for line in CONFIGFILE:

        valueLst=[]
        line = line.rstrip("\n\r")

        # Filter for comments and blank lines
        if not comment.match(line) and keyVal.match(line):

            # Weed out internal comments & split on 1st space
            line = comment.sub('',line)
            (keyword, value) = line.split(delim,1)

            # If the line contains a value
            keyword = keyword.strip()              # kill external whitespace
            keyword = spaces.sub('', keyword)      # kill internal whitespaces
            value = value.strip()                  # kill external whitespace
            value = spaces.sub(' ', value)         # shrink internal whitespace
            value = value.replace("'", '')         # kill quotes
            value = commaAndSpaces.sub(',', value) # kill ambiguous spaces

            # Split comma/space delimited value strings
            if doValueSplit:
                valueLst = commaOrSpace.split(value)
                if len(valueLst)<=1:
                    valueLst = valueLst[0]
                configTable[keyword] = valueLst
            else:
                configTable[keyword] = value

    return configTable


#-----------------------------------------------------------------------------#
def deg2dms(deg, delim=':', doSign=False, nPlaces=2):
    """
    Convert a float in degrees to 'dd mm ss' format.
    """

    try:
        angle = abs(deg)
        sign=1
        if angle!=0: sign = angle/deg
        
        # Calcuate the degrees, min and sec
        dd = int(angle)
        rmndr = 60.0*(angle - dd)
        mm = int(rmndr)
        ss = 60.0*(rmndr-mm)

        # If rounding up to 60, carry to the next term
        if float("%05.2f" % ss) >=60.0:
            mm+=1.0
            ss = ss - 60.0
        if float("%02d" % mm) >=60.0:
            dd+=1.0
            mm = mm -60.0
        if nPlaces> 0:
            formatCode = "%0" + "%s.%sf" % (str(2 + nPlaces + 1), str(nPlaces))
        else:
            formatCode = "%02.0f"
        if sign>0:
            if doSign:
                formatCode = "+%02d%s%02d%s" + formatCode
            else:
                formatCode = "%02d%s%02d%s" + formatCode
        else:
            formatCode = "-%02d%s%02d%s" + formatCode
        return formatCode % (dd, delim, mm, delim, ss)
        
    except Exception:
        return None

    
#-----------------------------------------------------------------------------#
def fail_not_exists(item, type='file', LF=None):
    """
    Check for a file or directory, exit if not found. Messages to STDOUT or a
    log file via a file handle.
    """
    if LF is None:
        LF = sys.stdout
    if not os.path.exists(item):
        log_fail(LF, "Err: The %s '%s' does not exist." % (type, item))
    else:
        log_wr(LF, "Found the %s '%s'." % (type, item))


#-----------------------------------------------------------------------------#
def set_statusfile(statusFile, status=0):
    """
    Write a text file containing an integer representing a status.
    """
    
    if os.path.exists(statusFile):
        os.remove(statusFile)
    SF = open(statusFile, 'w')
    SF.write("%d\n" % status)
    SF.close()


#-----------------------------------------------------------------------------#
def read_statusfile(statusFile):
    """
    Read a single line text file.
    """
    
    SF = open(statusFile, 'r')
    statStr = SF.readline().rstrip("\n\r")
    SF.close()
    
    return statStr


#-----------------------------------------------------------------------------#
def read_dictfile(dictFile):

    return json.load(open(dictFile, "r"))


#-----------------------------------------------------------------------------#
def write_dictfile(data, dictFile):
    
    json.dump(data, open(dictFile, "w"))



#-----------------------------------------------------------------------------#
def load_vector_fail(inFile, dtype='float32', tolist=False, LF=None):
    """
    Use the numpy loadtxt function to read a vector from a single-column
    text file. Can process multiple datatypes and return a numpy array or list. 
    """
    
    if LF is None:
        LF = sys.stdout

    try:
        a = np.loadtxt(inFile, dtype=dtype, ndmin=1)
        log_wr(LF, "Successfully loaded '%s' vector from '%s'." % (dtype,
                                                                   inFile))
    except Exception:
        log_fail(LF, "Err: Failed to load '%s' vector from '%s'." % (dtype,
                                                                     inFile))
    if tolist:
        return a.tolist()
    else:
        return a


#-----------------------------------------------------------------------------#
def cat_to_recarray(inCatFile, dtype, delim=" ", addUniqueName=False, 
                    doWarn=True, LF=None):
    """
    Read and parse the catalogue file. The 'dtype' argument should be a list
    of tuples, one tuple per column. Each tuple should have two entries: the
    column name and the dtype code (e.g., 'f8' for a float64, 'a10' for a
    10 character string). For example:    
    dtype=[('inName', 'a20'), ('rms', 'f8'), ('x_deg', 'f8'), ('y_deg', '<f8')]
    Lines which have the wrong number of entries or where entries fail to
    convert are skipped with a warning
    """
    if LF is None:
        LF = sys.stdout
        
    # Unpack the column names and formats
    try:
        colNameLst, dtypeLst = zip(*dtype)
        colNameLst = list(colNameLst)
        dtypeLst = list(dtypeLst)
    except Exception:
        log_fail(LF, "Err: Catalogue dtype formatter is inconsistant.")
        log_fail(LF, "dtype=%s" % dtype)

    # Check that the uniqueName field has been defined
    if not 'uniqueName' in colNameLst:
        log_fail(LF,
              "Err: Field 'uniqueName' missing from 'sourceCat' table schema.")
    else:
        indxUname = colNameLst.index('uniqueName')
                 
    # Compile a few useful regular expressions
    spaces = re.compile('\s+')
    comment = re.compile('#.*')
    
    DATFILE = open(inCatFile, "r")
    row = 0
    catLst = []
    for line in DATFILE:
        row += 1
        entryDict = {}
        line = line.rstrip("\n\r")
        if comment.match(line):            
            continue
        line = comment.sub('', line)     # remove internal comments
        line = line.strip()              # kill external whitespace
        line = spaces.sub(' ', line)     # shrink internal whitespace
        line = line.split(delim)
        
        # Add the uniqueName column
        if addUniqueName:
            line.insert(indxUname, '')
        
        nCols = len(line)
        if len(line)!=len(colNameLst):
            if doWarn:
                log_wr(LF,'Warn: Missmatch in number of catalogue entries.')
                log_wr(LF,'Row %d: #entries=%d, #columns=%d.' % (row, nCols,
                                                           len(colNameLst)))
            continue
        else:
            catLst.append(line)
    DATFILE.close()
    
    # Convert the 2D list to a ndarray and then to a record array
    if len(catLst)==0:
        log_fail(LF, 'Err: Zero valid entries in the catalogue file!')
    log_wr(LF,'Read %d entries from the catalogue file.' % len(catLst))
    catArr = np.array(catLst)
    catArr = np.core.records.fromarrays(catArr.transpose(), 
                                        names = ','.join(colNameLst),
                                        formats = ','.join(dtypeLst))
    return catArr


#-----------------------------------------------------------------------------#
def calc_sumbox_norm(beamFWHM_pix, side_pix):
    """
    Calculate the normalisation required for a spectrum summed over MxM pixels.
    """
    
    fNorm = 0.0
    start = int( -1 * m.floor(side_pix/2) )
    stop = int( m.floor(side_pix/2) + 1 )
    for i in range(start, stop):
        for j in range(start, stop):
            fNorm += m.exp( -1 * (i**2.0 + j**2.0) * 4.0 * m.log(2.0) / 
                            beamFWHM_pix**2.0)
            
    return 1.0 / fNorm


#-----------------------------------------------------------------------------#
def calc_mom2_FDF(ccFDF, phiArr):
    """
    Calculate the 2nd moment of the clean component spectrum.
    """
    
    K = np.sum( np.abs(ccFDF) )
    phiMean = np.sum( phiArr * np.abs(ccFDF) ) / K
    phiMom2 = np.sqrt( np.sum( np.power((phiArr - phiMean), 2.0) *
                                np.abs(ccFDF) ) / K )
    
    return phiMom2

    
#-----------------------------------------------------------------------------
def poly5(p):
    """
    Function which returns another function to evaluate a polynomial.
    """

    C5, C4, C3, C2, C1, C0 = p
    def rfunc(x):
        y = C5*x**5.0 + C4*x**4.0 + C3*x**3.0 + C2*x**2.0 + C1*x + C0
        return y
    return rfunc


#-----------------------------------------------------------------------------#
def log_fail(LF, errStr):
    """
    Record a fatal error in the log file, echo to STDOUT and exit.
    """
    
    print errStr
    if not LF==sys.stdout:
        LF.write(errStr + '\n')
        LF.close()
    sys.exit(1)


#-----------------------------------------------------------------------------#
def log_wr(LF, message):
    """
    Record a message in a log file and echo to STDOUT.
    """
    
    print message
    if not LF==sys.stdout:
        LF.write(message + '\n')


#-----------------------------------------------------------------------------#
def nanmedian(arr, **kwargs):
    """
    Returns median ignoring NANs.
    """
    
    return ma.median( ma.masked_where(arr!=arr, arr), **kwargs )


#-----------------------------------------------------------------------------#
def MAD(a, c=0.6745, axis=None):
    """
    Median Absolute Deviation along given axis of an array:
    median(abs(a - median(a))) / c
    c = 0.6745 is the constant to convert from MAD to std
    """
    
    a = ma.masked_where(a!=a, a)
    if a.ndim == 1:
        d = ma.median(a)
        m = ma.median(ma.fabs(a - d) / c)
    else:
        d = ma.median(a, axis=axis)
        if axis > 0:
            aswp = ma.swapaxes(a,0,axis)
        else:
            aswp = a
        m = ma.median(ma.fabs(aswp - d) / c, axis=0)

    return m


#-----------------------------------------------------------------------------#
def calc_stats(a, maskzero=False):
    """
    Calculate the statistics of an array.
    """
    
    statsDict = {}
    a = np.array(a)

    # Mask off bad values and count valid pixels
    if maskzero:
        a = np.where( np.equal(a, 0.0), np.nan, a)
    am = ma.masked_invalid(a)
    statsDict['npix'] = np.sum(~am.mask)
    
    if statsDict['npix']>=2:
        statsDict['stdev'] = float(np.std(am))
        statsDict['mean'] = float(np.mean(am))
        statsDict['median'] = float(nanmedian(am))
        statsDict['max'] = float(np.max(am))
        statsDict['min'] = float(np.min(am))
        statsDict['centmax'] = list(np.unravel_index(np.argmax(am),
                                                     a.shape))
        statsDict['madfm'] = float(MAD(am.flatten()))
        statsDict['success'] = True
        
    else:
        statsDict['npix'] == 0
        statsDict['stdev']   = 0.0
        statsDict['mean']    = 0.0
        statsDict['median']  = 0.0
        statsDict['max']     = 0.0
        statsDict['min']     = 0.0
        statsDict['centmax'] = (0.0, 0.0)
        statsDict['madfm']   = 0.0
        statsDict['success'] = False
        
    return statsDict


#-----------------------------------------------------------------------------#
def sort_nicely(l):
    """
    Sort a list in the order a human would.
    """
    
    convert = lambda text: int(text) if text.isdigit() else text 
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    l.sort( key=alphanum_key ) 


#-----------------------------------------------------------------------------#
def twodgaussian(params, shape):
    """
    Build a 2D Gaussian ellipse as parameterised by 'params' for a region with
    'shape'
        params - [amp, xo, yo, cx, cy, pa] where:
                amp - amplitude
                xo  - centre of Gaussian in X
                yo  - centre of Gaussian in Y
                cx  - width of Gaussian in X (sigma or c, not FWHM)
                cy  - width of Gaussian in Y (sigma or c, not FWHM)
                pa  - position angle of Gaussian, aka theta (radians)
        shape - (y, x) dimensions of region
    Returns a 2D numpy array with shape="shape" 
    """
    
    assert(len(shape) == 2)
    amp, xo, yo, cx, cy, pa = params
    y, x = np.indices(shape)
    st = m.sin(pa)**2
    ct = m.cos(pa)**2
    s2t = m.sin(2*pa)
    a = (ct/cx**2 + st/cy**2)/2
    b = s2t/4 *(1/cy**2-1/cx**2)
    c = (st/cx**2 + ct/cy**2)/2
    v = amp*np.exp(-1*(a*(x-xo)**2 + 2*b*(x-xo)*(y-yo) + c*(y-yo)**2))
    
    return v


#-----------------------------------------------------------------------------#
def create_pqu_spectra_RMthin(freqArr_Hz, fracPol, psi0_deg, RM_radm2):
    """Return fractional P/I, Q/I & U/I spectra for a Faraday thin source"""
    
    # Calculate the p, q and u Spectra
    pArr = fracPol * np.ones_like(freqArr_Hz)
    lamSqArr_m2 = np.power(C/freqArr_Hz, 2.0)
    pArr = fracPol * np.ones_like(lamSqArr_m2)
    quArr = pArr * np.exp( 2j * (np.radians(psi0_deg) +
                                 RM_radm2 * lamSqArr_m2 ) )
    qArr = quArr.real
    uArr = quArr.imag
    
    return pArr, qArr, uArr


#-----------------------------------------------------------------------------#
def create_IQU_spectra_RMthin(freqArr_Hz, fluxI, SI, fracPol, psi0_deg, 
                              RM_radm2, freq0_Hz=None):
    """Return Stokes I, Q & U spectra for a Faraday thin source"""

    pArr, qArr, uArr = create_pqu_spectra_RMthin(freqArr_Hz,
                                                 fracPol,
                                                 psi0_deg, 
                                                 RM_radm2)
    if freq0_Hz is None:
        freq0_Hz = freqArr_Hz[0]
    IArr = fluxI * np.power(freqArr_Hz/freq0_Hz, SI)
    PArr = IArr * pArr
    QArr = IArr * qArr
    UArr = IArr * uArr

    return IArr, QArr, UArr


#-----------------------------------------------------------------------------#
def create_pqu_resid_RMthin(qArr, uArr, freqArr_Hz, fracPol, psi0_deg,
                            RM_radm2):
    """Subtract a RM-thin component from the fractional q and u data."""

    pModArr, qModArr, uModArr = create_pqu_spectra_RMthin(freqArr_Hz,
                                                          fracPol,
                                                          psi0_deg,
                                                          RM_radm2)
    qResidArr = qArr - qModArr
    uResidArr = uArr - uModArr
    pResidArr = np.sqrt(qResidArr**2.0 + uResidArr**2.0)

    return pResidArr, qResidArr, uResidArr


#-----------------------------------------------------------------------------#
def xfloat(x, default=None):

    if x is None or x is "":
        return default
    try:
        return float(x)
    except Exception:
        return default
