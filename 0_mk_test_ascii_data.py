#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     mk_test_ascii_data.py                                             #
#                                                                             #
# USAGE:    ./mk_test_ascii_data.py                                           #
#                                                                             #
# PURPOSE:  Create a small ASCII dataset for the purposes of testing the      #
#           RM-pipeline. The script outputs a set ASCII files containing      #
#           frequency and Stokes vectors based on source parameters read from #
#           an input catalogue file. Properties of the output data are set    #
#           in variables at the top of this script. A template describing the #
#           shape of the rms noise curve may also be provided in an external  #
#           ASCII file [freq_Hz, amp].                                        #
#                                                                             #
# MODIFIED: 29-January-2016 by C. Purcell                                     #
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

# Frequency parameters
startFreq_Hz = 0.7e9
endFreq_Hz = 1.8e9
nChans = 1100

# Noise level
rmsNoise_mJy = 0.1

# Variations in the noise vs frequency can be specified in an external file.

# Properties of injected sources are given by an external CSV catalogue file.
# Two types of model may be specified, assuming a common flux & spectral index.

# Frequency at which the flux is specified in the catalogue
freq0_Hz = startFreq_Hz

# END USER EDITS -------------------------------------------------------------#

import os
import sys
import argparse
import shutil
import math as m
import numpy as np
import astropy.io.fits as pf
import astropy.wcs.wcs as pw

from Imports.util_PPC import create_IQU_spectra_burn
from Imports.util_PPC import create_IQU_spectra_diff
from Imports.util_PPC import csv_read_to_list
from Imports.util_PPC import split_repeat_lst 
from Imports.util_PPC import calc_stats
from Imports.util_RM import extrap

C = 2.99792458e8

#-----------------------------------------------------------------------------#
def main():
    """
    Start the create_IQU_ascii_data function if called from the command line.
    """
    
    # Help string to be shown using the -h option
    descStr = """
    Create a new dataset directory and populate it with ASCII files containing
    Stokes I, Q and U spectra. Each output file contains four columns
    corresponding to [freq_Hz, StokesI_Jy, StokesQ_Jy, StokesU_Jy] vectors for
    one source.

    The spectra are populated with polarised sources whose properties are given
    in an external CSV-format catalogue file. Two types of model may be 
    specified, assuming a common flux & spectral index:

        # MODEL TYPE 1: One or more components affected by Burn depolarisation.
        #
        # Column  |  Description
        #---------------------------------------------------
        # [0]     |  Model type (1)
        # [1]     |  X coordinate (deg)
        # [2]     |  Y coordinate (deg)
        # [3]     |  Major axis (arcsec)
        # [4]     |  Minor axis (arcsec)
        # [5]     |  Position angle (deg)
        # [6]     |  Total flux (mJy)
        # [7]     |  Spectral index
        # Component 1:
        # [8]     |  Intrinsic polarisation angle (deg)
        # [9]     |  Fractional polarisation
        # [10]    |  Faraday depth (radians m^-2)
        # [11]    |  Farday dispersion (radians m^-2)
        # Component 2:
        # [12]    |  Intrinsic polarisation angle (deg)
        # [13]    |  Fractional polarisation
        # [14]    |  Faraday depth (radians m^-2)
        # [15]    |  Farday dispersion (radians m^-2)
        # Component 3:
        # [16]    |  ...
        #---------------------------------------------------

        # MODEL TYPE 2: One or more stacked layers with differential Faraday
        # rotation (Sokoloff 1998, Eqn. 9).
        #
        # Column  |  Description
        #---------------------------------------------------
        # [0]     |  Model type (2)
        # [1]     |  X coordinate (deg)
        # [2]     |  Y coordinate (deg)
        # [3]     |  Major axis (arcsec)
        # [4]     |  Minor axis (arcsec)
        # [5]     |  Position angle (deg)
        # [6]     |  Total flux (mJy)
        # [7]     |  Spectral index
        # Component 1:
        # [8]     |  Intrinsic polarisation angle (deg)
        # [9]     |  Fractional polarisation
        # [10]    |  Faraday depth (radians m^-2)
        # Component 2:
        # [11]    |  Intrinsic polarisation angle (deg)
        # [12]    |  Fractional polarisation
        # [13]    |  Faraday depth (radians m^-2)
        # Component 3:
        # [14]    |  ...
        #---------------------------------------------------

    Properties of the data (frequency sampling, noise level) are given at the
    top of this script, including an optional template for the shape of the
    noise curve. 

    In addition to the ASCII files, the script outputs a simple ASCII catalogue
    and a SQL description of that catalogue. The catalogue file is used to
    drive the pipeline and the SQL descripton file tells the pipeline the
    format of the catalogue. This allows the user to define custom columns in
    the input catalogue, which are then incorporated into the results database.

    Example:

    ./0_mk_test_ascii_data.py catalogue.csv testASCIIData/
    """

    # Parse the command line options
    parser = argparse.ArgumentParser(description=descStr,
                                 formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("inCatFile", metavar="catalogue.csv", nargs=1,
                        help="Input catalogue file in CSV format")
    parser.add_argument("dataPath", metavar="PATH/TO/DATA",
                        default="testASCIIData/", nargs="?",
                        help="Path to new data directory [testASCIIData/]")
    parser.add_argument('-n', dest='noiseTmpFile', metavar="NOISE.TXT",
                        help="File providing a template noise curve")
    parser.add_argument('-f', dest='flagFreqStr', metavar='f1,f2,f1,f2,...',
                        default="", help="Frequency ranges to flag out")
    args = parser.parse_args()
    inCatFile = args.inCatFile[0]
    dataPath = args.dataPath
    noiseTmpFile = args.noiseTmpFile
    flagRanges_Hz = []
    if len(args.flagFreqStr)>0:
        try:
            flagFreqLst = args.flagFreqStr.split(",")
            flagFreqLst = [float(x) for x in flagFreqLst]
            flagRanges_Hz = zip(*[iter(flagFreqLst)]*2)
        except Exception:
            "Warn: Failed to parse frequency flagging string!"

    # Read the RMS noise template
    try:
        noiseTmpArr = np.loadtxt(noiseTmpFile, unpack=True)
    except Exception:
        noiseTmpArr = None
        print "Failed to load noise template '%s'." % noiseTmpFile
        print "Assuming flat noise profile."
        
    
    # Call the function to create the ASCII data files on disk
    nSrc = create_IQU_ascii_data(dataPath, inCatFile, startFreq_Hz,
                                 endFreq_Hz, nChans, rmsNoise_mJy,
                                 noiseTmpArr, flagRanges_Hz)

    # Print summary to user 
    sessionPath = "testSessASCII/"
    outCatFile = dataPath.rstrip("/") + "/testCat.txt"
    sqlFile = dataPath.rstrip("/") + "/testCatDesc.sql"
    print
    print "-" * 80
    print ">>> How to run the RM pipeline:"
    print "-" * 80
    print "A test dataset had been created in the directory '%s/':" \
          % dataPath.rstrip("/")
    print "> %d ASCII files [freq, I, Q, U]" \
          % nSrc
    print "> A simple catalogue in the file '%s' " % outCatFile
    print "> A SQL catalogue description in the file '%s'" % sqlFile
    print
    print "To run the RM-pipeline execute the following commands in order:"
    print
    print "./1_verify_ascii_data.py %s/" % dataPath.rstrip("/")
    print "./2_create_ascii_session.py %s/ %s/ %s %s" % \
          (sessionPath.rstrip("/"), dataPath.rstrip("/"),  outCatFile, sqlFile)
    print "# Edit the file '%s/inputs.config' (optional)" \
          % sessionPath.rstrip("/")
    print "./3_extract_spectra.py %s/"  % sessionPath.rstrip("/")
    print "./4_do_RM-synthesis.py %s/" % sessionPath.rstrip("/")
    print "./5_do_RM-clean.py %s/" % sessionPath.rstrip("/")
    print "./6_measure_complexity.py %s/" % sessionPath.rstrip("/")
    print
    print "NOTE: information and help on each script can be viewed by ",
    print "executing each\ncommand followed by a '-h' flag, e.g.: \n"
    print "./0_mk_test_ascii_data.py -h"
    print 


#-----------------------------------------------------------------------------#
def create_IQU_ascii_data(dataPath, inCatFile, startFreq_Hz, endFreq_Hz, 
                          nChans, rmsNoise_mJy, noiseTmpArr=None,
                          flagRanges_Hz=[]):
    """Create a set of ASCII files containing Stokes I Q & U spectra."""

    # Sample frequency space and flag the bad frequency ranges   
    freqArr_Hz = np.linspace(startFreq_Hz, endFreq_Hz, nChans)
    for i in range(len(freqArr_Hz)):
        for fRng in flagRanges_Hz:            
            if freqArr_Hz[i]>=fRng[0] and freqArr_Hz[i]<=fRng[1]:
                freqArr_Hz[i]=np.nan
    freqArr_Hz = freqArr_Hz[np.where(freqArr_Hz==freqArr_Hz)]

    # Create normalised noise array from a template or assume all ones.
    if noiseTmpArr is None:
        noiseArr = np.ones(freqArr_Hz.shape, dtype="f8")
    else:
        xp = noiseTmpArr[0]
        yp = noiseTmpArr[1]
        mDict = calc_stats(yp)
        yp /= mDict["median"]
        noiseArr = extrap(freqArr_Hz, xp, yp)
        
    # Check the catalogue file exists
    if not os.path.exists(inCatFile):
        print "Err: File does not exist '%s'." % inCatFile
        sys.exit()
    catInLst = csv_read_to_list(inCatFile, doFloat=True)

    # Create the output directory path
    dataPath = dataPath.rstrip("/")
    print "Creating test dataset in '%s/'" % dataPath
    dirs = dataPath.split("/")
    for i in range(1, len(dirs)):
        dirStr = "/".join(dirs[:i])
        if not os.path.exists(dirStr):
            os.mkdir(dirStr)
    if os.path.exists(dataPath):
        shutil.rmtree(dataPath, True)
    os.mkdir(dataPath)

    # Create a catalogue description file
    sqlFile = dataPath + "/testCatDesc.sql"
    sqlFH = open(sqlFile, "w")
    descStr = """
CREATE TABLE sourceCat (
uniqueName varchar(50),
fileName varchar(50),
x_deg double,
y_deg double);
    """
    sqlFH.write("%s\n" % descStr)
    sqlFH.close()

    # Open an output catalogue file
    outCatFile = dataPath + "/testCat.txt"
    outCatFH = open(outCatFile, "w")
    outCatFH.write("#Name fileName x_deg  y_deg\n")
    
    # Loop through the sources, calculate the spectra and save to disk
    successCount = 0
    for i in range(len(catInLst)):
        e = catInLst[i]
        modelType = int(e[0])

        # Type 1 = multiple Burn depolarisation affected components
        if modelType==1:
            
            # Parse the parameters of multiple components
            preLst, parmArr = split_repeat_lst(e[1:],7,4)
            
            # Create the model spectra from multiple thin components
            # modified by external depolarisation
            IArr_Jy, QArr_Jy, UArr_Jy = \
                create_IQU_spectra_burn(freqArr_Hz = freqArr_Hz,
                                        fluxI = preLst[5]/1e3, # mJy->Jy
                                        SI = preLst[6],
                                        fracPolArr = parmArr[0],
                                        psi0Arr_deg = parmArr[1],
                                        RMArr_radm2 = parmArr[2],
                                        sigmaRMArr_radm2 = parmArr[3],
                                        freq0_Hz = freq0_Hz)
                
        # Type 2 = multiple internal depolarisation affected components
        elif modelType==2:
            
            # Parse the parameters of multiple components
            preLst, parmArr = split_repeat_lst(e[1:],7,3)
            
            # Create the model spectra from multiple components
            # modified by internal Faraday depolarisation
            IArr_Jy, QArr_Jy, UArr_Jy = \
                create_IQU_spectra_diff(freqArr_Hz = freqArr_Hz,
                                        fluxI = preLst[5]/1e3, # mJy->Jy
                                        SI = preLst[6],
                                        fracPolArr = parmArr[0],
                                        psi0Arr_deg = parmArr[1],
                                        RMArr_radm2 = parmArr[2],
                                        freq0_Hz = freq0_Hz)
        else:
            continue
        
        # Add scatter to the data to simulate noise
        rmsNoise_Jy = rmsNoise_mJy/1e3
        IArr_Jy += (np.random.normal(scale=rmsNoise_Jy, size=IArr_Jy.shape)
                    * noiseArr)
        QArr_Jy += (np.random.normal(scale=rmsNoise_Jy, size=QArr_Jy.shape)
                    * noiseArr)
        UArr_Jy += (np.random.normal(scale=rmsNoise_Jy, size=UArr_Jy.shape)
                    * noiseArr)
        dIArr_Jy = noiseArr * rmsNoise_Jy
        dIArr_Jy *= np.random.normal(loc=1.0, scale=0.05, size=noiseArr.shape)
        dQArr_Jy = noiseArr * rmsNoise_Jy 
        dQArr_Jy *= np.random.normal(loc=1.0, scale=0.05, size=noiseArr.shape)
        dUArr_Jy = noiseArr * rmsNoise_Jy 
        dUArr_Jy *= np.random.normal(loc=1.0, scale=0.05, size=noiseArr.shape)
        
        # Save spectra to disk
        outFileName = "Source%d.dat" % (i+1)
        outFilePath = dataPath + "/" + outFileName
        print "Writing ASCII file '%s' ..." % outFileName,
        np.savetxt(outFilePath,
                   np.column_stack((freqArr_Hz, IArr_Jy, QArr_Jy, UArr_Jy,
                                    dIArr_Jy, dQArr_Jy, dUArr_Jy)))
        print "done."
        
        # Add to the catalogue file
        outCatFH.write("Source%d %s %f %f\n" % ((i+1),
                                                outFileName,
                                                preLst[0],
                                                preLst[1]))
        successCount += 1

    # Clean up
    outCatFH.close()
    
    return successCount



#-----------------------------------------------------------------------------#
if __name__=="__main__":
    main()
