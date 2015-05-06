#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     util_PPC.py                                                       #
#                                                                             #
# PURPOSE:  Helper functions for the POSSUM pipeline.                         #
#                                                                             #
# REQUIRED: Requires the numpy and astropy.                                   #
#                                                                             #
# MODIFIED: 30-Apr-2015 by C. Purcell                                         #
#                                                                             #
# CONTENTS:                                                                   #
#                                                                             #
#  PipelineInputs       ... class to contain the pipeline inputs and file     #
#  config_read          ... read a key=value format text file                 #
#  deg2dms              ... convert decimal degrees to dms string             #
#  fail_not_exists      ... check for dir/file existence, exit if missing     #
#  set_statusfile       ... set text in an ASCII status file                  #
#  read_statusfile      ... read text in an ASCII status file                 #
#  read_paths           ... read the paths from an external file              #
#  load_vector_fail     ... load a single-column vector from a text file      #
#  cat_to_recarray      ... convert ASCII catalogue file to a record array    #
#  extract_spec_noise   ... extract source and noise spectra from FITS file(s)#
#  calc_sumbox_norm     ... calculate the beam normalisation for a box MxM    #
#  calc_mom2_FDF        ... calculate the 2nd moment of the CC                #
#  poly5                ... function to evaluate a 5th order polynomial       #
#  log_fail             ... record a fatal error in the log file and exit     #
#  log_wr               ... record a message in the log file                  #
#  inparm_verify        ... verify presence of required pipeline inputs       #
#  spawnDaemon          ... call an executable file as a detached procress    #
#  nanmedian            ... np.median ignoring NaNs                           #
#  MAD                  ... calculate the madfm                               #
#  calc_stats           ... calculate the statistics of an array              #
#  calc_clipped_stats   ... calculate the stats after sigma-clipping          #
#  sort_nicely          ... sort a list in the order a human would            #
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

C = 2.997924538e8 # Speed of light [m/s]

#-----------------------------------------------------------------------------#
class PipelineInputs:
    """
    Class to hold the driving inputs and defaults for the POSSUM pipeline.
    Methods assume unique keys across all configuration sections. There are
    two levels of parameters: 1) input parameters linked to a file and stored
    in a ConfigParser object, 2) derived parameters stored in a dictionary.
    Derived parameters pertain to the frequency and Faraday depth sampling and
    are set to sensible

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
        Hard-coded defaults for the RM pipeline.
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
        self.config.set("Thresholds", "phiMom2Thresh_sigma", "1.5")
        
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
        Return a dictionary of key=val assuming unique keys across all sections.
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
        #dataTypeFile = dataPath + '/dataType.txt'
        dataTypeFile = self.sessionDir + '/dataType.txt'
        dataType = read_statusfile(dataTypeFile)

        # Frequency and wavelength sampling
        # These are set by the dataset and cannot be altered by the user
        #inFreqFile = dataPath + '/freqs_Hz.txt'
        inFreqFile =  self.sessionDir + '/freqs_Hz.txt'
        freqArr_Hz = load_vector_fail(inFreqFile, "float32", False)
        dFreq_Hz = np.nanmin(np.abs(np.diff(freqArr_Hz)))
        lambdaArr_m = C / freqArr_Hz
        lambdaSqArr_m2 = np.power(lambdaArr_m, 2.0)
        lambdaSqRange_m2 = ( np.nanmax(lambdaSqArr_m2) -
                             np.nanmin(lambdaSqArr_m2) )        
        dLambdaSqMin_m2 = np.nanmin(np.abs(np.diff(lambdaSqArr_m2)))

        # Default Faraday depth limits
        # The Faraday sampling can be forced by existing config options
        fwhmRMSF_radm2 = 2.0 * m.sqrt(3.0) / lambdaSqRange_m2
        if resetPhiSamp:
            oversampling = float(self.config.get("RMsynthesis", "oversampling"))
            dPhi_radm2 = fwhmRMSF_radm2 / oversampling
            phiMax_radm2 = m.sqrt(3.0) / dLambdaSqMin_m2
        else:
            dPhi_radm2 = float(self.config.get("RMsynthesis", "dPhi_radm2"))
            phiMax_radm2 = float(self.config.get("RMsynthesis", "phiMax_radm2"))
            
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
                          "phiMom2Thresh_sigma"]
        pDict = self.get_flat_dict(includeDerived=False)
        keys = pDict.keys()
        missingLst = [x for x in requiredKeyLst if x not in keys]

        if len(missingLst) > 0:
            return missingLst
        else:
            return []
        

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
def read_paths(pathFile):
    """
    Read the paths in a text file and return a dictionary of key=val pairs.
    """

    # Default required paths to current directory
    varDict = {'pipeDir':'.', 'cgiURL':'.'}
    
    if os.path.exists(pathFile):
        F = open(pathFile)
        for line in F:
            line = line.rstrip("\n\r")
            (keyword, value) = line.split('=', 1)
            varDict[keyword] = value
        F.close()

    return varDict


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
def cat_to_recarray(inCatFile, dtype, delim=" ", LF=None):
    """
    Read and parse the catalogue file. The 'dtype' argument should be a list
    of tuples, one tuple per column. Each tuple should have two entries: the
    column name and the dtype code (e.g., 'f8' for a float64, 'a10' for a
    10 character string). For example:    
    dtype=[('inName', 'a20'), ('rms', 'f8'), ('x_deg', 'f8'), ('y_deg', '<f8')]
    Lines which have the wrong number of entries or where entries fail to
    convert are skipped.
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
        line.insert(indxUname, '')
        
        nCols = len(line)
        if len(line)!=len(colNameLst):
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
    Calculate the 2nd moment of the clean component distribution.
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
def inparm_verify(parmDict):
    """
    Verify that the input parameters to the pipeline are all present
    """
    
    requiredKeyLst = ['dataPath',
                      'inFreqFile',
                      'gain',
                      'sumBoxPix',
                      'dLambdaSqMin_m2',
                      'phiCentre',
                      'lambdaSqRange',
                      'dPhi',
                      'nChanRM',
                      'weightType',
                      'thresholdSignalPI',
                      'thresholdPolBias',
                      'thresholdDoClean',
                      'phiMom2Thresh',
                      'cleanAlgorithm',
                      'maxCleanIter',
                      'cleanCutoff_sigma']

    keys = parmDict.keys()
    missingLst = [x for x in requiredKeyLst if x not in keys]

    if len(missingLst) > 0:
        return missingLst
    else:
        return []

    
#-----------------------------------------------------------------------------#
def spawnDaemon(path_to_executable, *args):
    """Spawn a completely detached subprocess (i.e., a daemon).
    E.g.:
    spawnDaemon('../bin/producenotify.py', 'producenotify.py', 'xx')
    """

    # The standard I/O file descriptors are redirected to /dev/null.
    if (hasattr(os, "devnull")):
       REDIRECT_TO = os.devnull
    else:
       REDIRECT_TO = "/dev/null"
    
    # fork the first time (to make a non-session-leader child process)
    try:
        pid = os.fork()
    except OSError, e:
        raise RuntimeError("1st fork failed: %s [%d]" % (e.strerror, e.errno))
    if pid != 0:
        # parent (calling) process is all done
        return

    # detach from controlling terminal (to make child a session-leader)
    os.setsid()
    try:
        pid = os.fork()
    except OSError, e:
        raise RuntimeError("2nd fork failed: %s [%d]" % (e.strerror, e.errno))
        raise Exception, "%s [%d]" % (e.strerror, e.errno)
    if pid != 0:
        # child process is all done
        os._exit(0)
        
    # grandchild process now non-session-leader, detached from parent
    # grandchild process must now close all open files
    try:
        maxfd = os.sysconf("SC_OPEN_MAX")
    except (AttributeError, ValueError):
        maxfd = 1024

    for fd in range(0,maxfd):
        try:
           os.close(fd)
        except OSError: # ERROR, fd wasn't open to begin with (ignored)
           pass

    # redirect stdin, stdout and stderr to /dev/null
    os.open(REDIRECT_TO, os.O_RDWR) # standard input (0)
    os.dup2(0, 1)
    os.dup2(0, 2)

    # and finally let's execute the executable for the daemon!
    try:
       
       os.execv(path_to_executable, (path_to_executable,) +  tuple(args))
    except Exception, e:
       # oops, we're cut off from the world, let's just give up
       os._exit(255)



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
        # I don't want the array to change so I have to copy it?
        if axis > 0:
            aswp = ma.swapaxes(a,0,axis)
        else:
            aswp = a
        m = ma.median(ma.fabs(aswp - d) / c, axis=0)

    return m


#-----------------------------------------------------------------------------#
def calc_stats(a, maskzero=False):
    """
    Calculate the statistics of an array (masked version).
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
def calc_clipped_stats(data, clip=3.0, nIter=10, maskzero=False):
    """
    Calculate the statistics of an array given a sigma clip.
    """
    
    if maskzero:
        data = np.where( np.equal(data, 0.0), np.nan, data)
        
    ms = calc_stats_ma(data)

    if clip>0 and nIter>0:
        convergeFlg = 0
        itCnt = 0
        while convergeFlg==0 and itCnt<nIter:
            meanOld, stdOld, madOld = ms['mean'], ms['stdev'], ms['madfm']
            minVal = ms['mean'] - (clip * ms['madfm'])
            maxVal = ms['mean'] + (clip * ms['madfm'])
            
            # Blank values outside the 3-sigma range
            dataMsk = np.where(np.greater(data, maxVal), np.nan, data)
            dataMsk = np.where(np.less(data, minVal), np.nan, dataMsk)
            
            # Measure the statistics
            ms = calc_stats(dataMsk)
            dataMsk = []
    
            if ms['mean'] == meanOld and ms['madfm'] == madOld:
                convergFlg = 1
            itCnt += 1

    return ms


#-----------------------------------------------------------------------------#
def sort_nicely(l):
    """
    Sort a list in the order a human would.
    """
    
    convert = lambda text: int(text) if text.isdigit() else text 
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    l.sort( key=alphanum_key ) 

