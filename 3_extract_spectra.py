#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     3_extract_spectra.py                                              #
#                                                                             #
# USAGE:    ./3_extract_spectra.py PATH/TO/SESSION                            #
#                                                                             #
# PURPOSE:  Extract spectra from the I, Q & U FITS files linked to a          #
#           RM pipeline session.                                              #
#                                                                             #
# MODIFIED: 19-November-2015 by C. Purcell                                    #
#                                                                             #
# TODO:                                                                       #
#                                                                             #
#  * Check that blanked (NaN) planes/pixels are dealt with correctly.         #
#    - certainly not dealt with in mpfit.                                     #
#  * Implement an LPF alternative, as a fallback mode if the I fit fails?     #
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
import shutil
import time
import argparse
import math as m
import numpy as np
import sqlite3

from Imports.util_PPC import PipelineInputs
from Imports.util_PPC import fail_not_exists
from Imports.util_PPC import log_wr
from Imports.util_PPC import log_fail
from Imports.util_PPC import load_vector_fail
from Imports.util_PPC import read_dictfile
from Imports.util_PPC import write_dictfile

from Imports.util_DB import register_sqlite3_numpy_dtypes
from Imports.util_DB import select_into_arr
from Imports.util_DB import insert_arr_db

from Imports.module_spec_extract_area import mod_spec_extract

# Constants
C = 2.99792458e8

# Turn off print statements buffering
sys.stdout = os.fdopen(sys.stdout.fileno(), "w", 0)

# Map numpy data-types to the limited sqlite3 data-types
register_sqlite3_numpy_dtypes()


#-----------------------------------------------------------------------------#
def main():
    """
    Start the run_spectral_extraction procedure if called from the command line.
    """
    
    # Help string to be shown using the -h option
    descStr = """
    Loop through the FITS files linked to the current session and extract I, Q
    and U data at the positions of sources. The catalogue should already
    exist in the database, loaded by the 'create_image_session.py' script.
    That script has also created a 'PATH/TO/SESSION/inputs.config' file, used
    to set the pipeline input parameters. Data extracted for each source
    will be saved to a directory called 'PATH/TO/SESSION/OUT'. For each source
    in the catalogue the following data are saved to a FITS format file:
      * A cube centred on each source, or offset if abutting an edge.
      * A single-plane mask image showing the extraction aperture.
      * A one dimentional spectrum from the source, RMS and frequency axis.
    If the output files already exist, default behaviour is to redo the
    measurements of the spectra. 

    Note: The measurements on the spectra are saved to the SQLite database in
    the file 'PATH/TO/SESSION/session.sqlite'.
    
    Example:
    
    ./3_extract_spectra.py testSession/
    """

    # Parse the command line options
    parser = argparse.ArgumentParser(description=descStr,
                                 formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("sessionPath", metavar="PATH/TO/SESSION", nargs=1,
                        help="Path to the new session directory [no default]")
    parser.add_argument("-o", dest="doOverwrite", action="store_true",
                        help="Overwrite previously extracted files.")
    parser.add_argument("-r", dest="doReset", action="store_true",
                        help="Completely reset the OUT/ directory.")
    args = parser.parse_args()
    sessionPath = args.sessionPath[0]
    doOverwrite = args.doOverwrite
    doReset = args.doReset

    # Call the spectral extraction function
    run_spectral_extraction(sessionPath, doOverwrite, doReset)


#-----------------------------------------------------------------------------#
def run_spectral_extraction(sessionPath, doOverwrite=False, doReset=False):

    sessionPath = sessionPath.rstrip("/")

    # Check the required directory structure exists or exit
    fail_not_exists(sessionPath, "directory")
    
    # Open the existing logfile to append to
    logFile = sessionPath + "/pipeline.log"
    LF = open(logFile, "a", 0)
    log_wr(LF, "\n>>> Beginning spectral extraction process.")
    log_wr(LF, "Time: %s" % time.ctime() )

    # Check that the session has been created successfully
    statusFile = sessionPath + "/status.json"
    fail_not_exists(statusFile, "file", LF)
    statusDict = read_dictfile(statusFile)
    if int(statusDict["session"])<1:
        log_fail(LF, "Err: Session status file reports session was not " + \
                     "created successfully.")

    # Read and parse the pipeline input file
    inParmFile = sessionPath + "/inputs.config"
    fail_not_exists(inParmFile, "file", LF)
    try:
        pipeInpObj = PipelineInputs(inParmFile)
        log_wr(LF, "Successfully parsed the input parameter file.")
    except Exception:
        log_fail(LF, "Err: Failed to parse the input parameter file.")

    # Verify that the parameter file has all the correct entries
    missingLst = pipeInpObj.inparm_verify()
    if len(missingLst)>0:
        log_fail(LF, "Err: Required input parameters missing. - %s" %
                 missingLst)
    else:
        log_wr(LF, "All required input parameters present.")
    pDict = pipeInpObj.get_flat_dict()
        
    # Check that the specified dataset exists
    dataPath = pDict["dataPath"].rstrip("/")
    fail_not_exists(dataPath, "directory", LF)

    # Read in the lists of FITS files from the datasets list file
    dataListFileI = dataPath + "/fileLstI.txt"
    fail_not_exists(dataListFileI, "file", LF)
    fitsLstI = load_vector_fail(dataListFileI, "str", True, LF)
    fitsLstI = [dataPath + "/" + x for x in fitsLstI]
    dataListFileQ = dataPath + "/fileLstQ.txt"
    fail_not_exists(dataListFileQ, "file", LF)
    fitsLstQ = load_vector_fail(dataListFileQ, "str", True, LF)
    fitsLstQ = [dataPath + "/" + x for x in fitsLstQ]
    dataListFileU = dataPath + "/fileLstU.txt"
    fail_not_exists(dataListFileU, "file", LF)
    fitsLstU = load_vector_fail(dataListFileU, "str", True, LF)
    fitsLstU = [dataPath + "/" + x for x in fitsLstU]

    # Check that the correct numbers of planes are reported
    if len(fitsLstI)<2 or len(fitsLstQ)<2 or len(fitsLstU)<2:
        log_wr(LF, "Err: Dataset reports too few FITS files:")
        log_fail(LF, "I(%d), Q(%s), U(%s)" % (len(fitsLstI), len(fitsLstQ),
                                              len(fitsLstU)))

    # Get the frequency vector from the text file
    inFreqFile = dataPath + "/freqs_Hz.txt"
    fail_not_exists(inFreqFile, "file", LF)
    freqArr_Hz = load_vector_fail(inFreqFile, "float32", False, LF)
    
    # (Re)create the output directory if necessary
    outDataDir = sessionPath + "/OUT"
    if os.path.exists(outDataDir):
        log_wr(LF, "Directory '%s' exists." % outDataDir)
        if doReset:
            log_wr(LF, "Resetting existing OUT directory.")
            shutil.rmtree(outDataDir, True)
            os.mkdir(outDataDir)
    else:
        log_wr(LF, "Creating new OUT directory.")
        os.mkdir(outDataDir)

    # Connect to the database
    dbFile = sessionPath + "/session.sqlite"
    try:
        log_wr(LF, "> Connecting to existing DB file '%s' ..." % dbFile)
        conn = sqlite3.connect(dbFile)
        cursor = conn.cursor()
    except Exception:
        log_fail(LF, "Err: Failed to connect to '%s'." % dbFile)
    
    # Load the coordinates from the catalogue into memory
    sql = """
    SELECT uniqueName, x_deg, y_deg
    FROM sourceCat
    """
    log_wr(LF, "Loading the catalogue from the database.")
    catRec = select_into_arr(cursor, sql)
    nRows = len(catRec)
    if nRows==0:
        log_fail(LF, "Err: query returned zero entries. Exiting.")
    log_wr(LF, "%s rows returned." % nRows)
    
    # RUN SPECTRAL EXTRACTION MODULE -----------------------------------------#
    specRec = mod_spec_extract(catRec,
                               fitsLstI,
                               fitsLstQ,
                               fitsLstU,
                               freqArr_Hz=freqArr_Hz,
                               sumBox_pix=int(pDict["sumBoxPix"]),
                               polyOrd=int(pDict["polyOrd"]),
                               outDataDir=outDataDir,
                               doOverwrite=doOverwrite,
                               LF=LF)
        
    # Write the relevant tables to the database
    insert_arr_db(cursor, specRec, "spectraParms")
    conn.commit()
    log_wr(LF, "Database updated with spectral parameters.")

    # Close the connection to the database
    cursor.close()
    conn.close()

    # Update the status file to reflect successful extraction
    statusDict["extract"] = 1
    write_dictfile(statusDict, sessionPath + "/status.json")
    

#-----------------------------------------------------------------------------#
if __name__=="__main__":
    main()

