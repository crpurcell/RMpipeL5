#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     mk_test_ascii_data.py                                             #
#                                                                             #
# USAGE:    ./mk_test_ascii_data.py                                           #
#                                                                             #
# PURPOSE:  Create a small ASCII dataset for the purposes of testing the      #
#           RM-pipeline. The script outputs a set ASCII files containing      #
#           frequency and Stokes vesctors. Edit the values at the top of the  #
#           script and run.                                                   #
#                                                                             #
# MODIFIED: 24-November-2015 by C. Purcell                                    #
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

# Example session path
sessionPath = "testSessASCII/"

# Data parameters
startFreq_Hz =  1.0e9
endFreq_Hz = 3.0e9
nChans = 30
rmsNoise_mJy = 0.1

# Properties of injected Faraday thin point sources
#-------------------------------------------------------------#
# [0]        [1] [2]      [3]       [4]      
# fluxI_mJy, SI, fracPol, psi0_deg, RM_radm2
#-------------------------------------------------------------#
srcInLst = [ [5.0, -0.7, 0.2, 30.0, 19.0],
             [10.0, -0.1, 0.6, 80.0, 10.0],
             [15.0, -0.5, 0.7, 10.0, -50.0],
             [9.0,  0.0, 0.1, 0.0, 60.0],
             [7.0, -0.2, 0.5, 45.0, -32.0],
             [2.3, +0.5, 0.7, 120.0, -90.0],
             [0.3, +0.0, 0.1, 120.0, -90.0]]
#-------------------------------------------------------------#
freq0_Hz = startFreq_Hz  # Frequency at which the flux is specified

# END USER EDITS -------------------------------------------------------------#

import os
import sys
import argparse
import shutil
import math as m
import numpy as np
import astropy.io.fits as pf
import astropy.wcs.wcs as pw

from Imports.util_PPC import twodgaussian
from Imports.util_PPC import create_IQU_spectra_RMthin
from Imports.util_FITS import strip_fits_dims
from Imports.util_FITS import create_simple_fits_hdu

C = 2.99792458e8

#-----------------------------------------------------------------------------#
def main():
    """
    Start the create_IQU_ascii_data function if called from the command line.
    """
    
    # Help string to be shown using the -h option
    descStr = """
    Create a new dataset directory and populate it with ASCII files containing
    Stokes I, Q and U spectra. Each file contains four columns corresponding to
    [freq_Hz, StokesI_Jy, StokesQ_Jy, StokesU_Jy] vectors for one source.
    The spectra are populated with Faraday thin sources defined in a table at
    the the top of this script:

      [ [fluxI_mJy, SI, fracPol, psi0_deg, RM_radm2], ... ]
    
      fluxI_mJy ... flux of source in first channel (mJy)
      SI        ... frequency spectral index
      fracPol   ... fractional polarisation (constant across frequency range)
      psi0_deg  ... intrinsic polarisation angle )
      RM_radm2  ... rotation measure (rad/m^2)

    Please edit the variables at the top of the script to change the properties
    of the output data.

    In addition to the ASCII files, the script outputs a simple ASCII catalogue
    and a SQL description of that catalogue. The catalogue file is used to
    drive the pipeline and the SQL descripton file tells the pipeline the
    format of the catalogue. This allows the user to define custom columns in
    the input catalogue, which are then incorporated into the results database.

    Example:

    ./0_mk_test_ascii_data.py testASCIIData/
    """

    # Parse the command line options
    parser = argparse.ArgumentParser(description=descStr,
                                 formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("dataPath", metavar="PATH/TO/DATA",
                        default="testASCIIData/", nargs="?",
                        help="Path to new data directory [testASCIIData/]")
    args = parser.parse_args()
    dataPath = args.dataPath

    # Call the function to create FITS data
    create_IQU_ascii_data(dataPath, srcInLst, startFreq_Hz, endFreq_Hz, nChans,
                          rmsNoise_mJy)

    # Print summary to user
    nSrc = len(srcInLst)
    catFile = dataPath.rstrip("/") + "/testCat.txt"
    sqlFile = dataPath.rstrip("/") + "/testCatDesc.sql"
    print
    print "-" * 80
    print ">>> How to run the RM pipeline:"
    print "-" * 80
    print "A test dataset had been created in the directory '%s/':" \
          % dataPath.rstrip("/")
    print "> %d ASCII files [freq, I, Q, U]" \
          % nSrc
    print "> A simple catalogue in the file '%s' " % catFile
    print "> A SQL catalogue description in the file '%s'" % sqlFile
    print
    print "To run the RM-pipeline execute the following commands in order:"
    print
    print "./1_verify_ascii_data.py %s/" % dataPath.rstrip("/")
    print "./2_create_ascii_session.py %s/ %s/ %s %s" % \
          (sessionPath.rstrip("/"), dataPath.rstrip("/"),  catFile, sqlFile)
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
def create_IQU_ascii_data(dataPath, srcLst, startFreq_Hz, endFreq_Hz, nChans,
                          rmsNoise_mJy, x_deg=0.0, y_deg=0.0):
    """
    Create a set of ASCII files containing Stokes I Q & U spectra.
    """

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

    # Open a catalogue file
    catFile = dataPath + "/testCat.txt"
    catFH = open(catFile, "w")
    catFH.write("#Name fileName x_deg  y_deg\n")
    
    # Loop through the sources, calculate the spectra and save to disk
    freqArr_Hz = np.linspace(startFreq_Hz, endFreq_Hz, nChans)
    for i in range(len(srcLst)):
        IArr_Jy, QArr_Jy, UArr_Jy = \
                 create_IQU_spectra_RMthin(freqArr_Hz,
                                           srcLst[i][0]/1e3,  # fluxI_mJy -> Jy
                                           srcLst[i][1],      # SI
                                           srcLst[i][2],      # fracPol
                                           srcLst[i][3],      # psi0_deg
                                           srcLst[i][4],      # RM_radm2
                                           freq0_Hz)
        
        # Add the Gaussian noise 
        IArr_Jy += np.random.normal(scale=rmsNoise_mJy/1e3,  # mJy -> Jy
                                    size=IArr_Jy.shape)
        QArr_Jy += np.random.normal(scale=rmsNoise_mJy/1e3,  # mJy -> Jy
                                    size=QArr_Jy.shape)
        UArr_Jy += np.random.normal(scale=rmsNoise_mJy/1e3,  # mJy -> Jy
                                    size=UArr_Jy.shape)
        dIArr_Jy = np.ones_like(IArr_Jy) * rmsNoise_mJy/1e3  
        dQArr_Jy = np.ones_like(QArr_Jy) * rmsNoise_mJy/1e3
        dUArr_Jy = np.ones_like(UArr_Jy) * rmsNoise_mJy/1e3
        
        # Save spectra to disk
        outFileName = "Source%d.dat" % (i+1)
        outFilePath = dataPath + "/" + outFileName
        print "Writing ASCII file '%s' ..." % outFileName,
        np.savetxt(outFilePath,
                   np.column_stack((freqArr_Hz, IArr_Jy, QArr_Jy, UArr_Jy,
                                    dIArr_Jy, dQArr_Jy, dUArr_Jy)))
        print "done."
        
        # Add to the catalogue file
        catFH.write("Source%d %s %f %f\n" % ((i+1), outFileName, x_deg, y_deg))

    # Clean up
    catFH.close()


#-----------------------------------------------------------------------------#
if __name__=="__main__":
    main()
