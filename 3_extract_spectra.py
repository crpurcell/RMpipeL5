#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     3_extract_spectra.py                                              #
#                                                                             #
# USAGE:    ./3_extract_spectra.py PATH/TO/SESSION                            #
#                                                                             #
# PURPOSE:  Extract spectra from the I, Q & U FITS files linked to a          #
#           POSSUM pipeline session.                                          #
#                                                                             #
# MODIFIED: 19-May-2015 by C. Purcell                                         #
#                                                                             #
# TODO:                                                                       #
#                                                                             #
#  * Implement an LPF alternative, as a fallback mode if the I fit fails.     #
#  * Check that blanked (NaN) planes/pixels are dealt with correctly.         #
#    - certainly not dealt with in mpfit.                                     #
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

from Imports.util_DB import register_sqlite3_numpy_dtypes
from Imports.util_DB import select_into_arr
from Imports.util_DB import insert_arr_db

from Imports.module_spec_extract import mod_spec_extract

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
    and U spectra at the positions of sources. The catalogue should already
    exist  in the database, loaded by the 'create_image_session.py' script.
    That script has also created a 'PATH/TO/SESSION/inputs.config' file, used
    to set the pipeline input parameters. Spectra extracted for each source
    will be saved to a directory called PATH/TO/SESSION/OUT, which is
    optionally overwritten if it already exists. Each spectra is saved as an
    ASCII vector file named for the position on the sky.
    """

    # Parse the command line options
    parser = argparse.ArgumentParser(description=descStr)
    parser.add_argument("sessionPath", metavar="PATH/TO/SESSION", nargs=1,
                        help="Path to the new session directory [no default]")
    parser.add_argument("-o", dest="doOverwrite", action="store_true",
                        help="Overwrite old extraction files")
    args = parser.parse_args()
    sessionPath = args.sessionPath[0]
    doOverwrite = args.doOverwrite

    # Call the spectral extraction function
    run_spectral_extraction(sessionPath, doOverwrite)


#-----------------------------------------------------------------------------#
def run_spectral_extraction(sessionPath, doOverwrite=False):
    
    sessionPath = sessionPath.rstrip("/")

    # Check the required directory structure exists or exit
    fail_not_exists(sessionPath, "directory")
    
    # Open the existing logfile to append to
    logFile = sessionPath + "/pipeline.log"
    LF = open(logFile, "a", 0)
    log_wr(LF, "\n>>> Beginning spectral extraction process.")
    log_wr(LF, "Time: %s" % time.ctime() )

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
        if doOverwrite:
            log_wr(LF, "Removing existing OUT directory.")
            shutil.rmtree(outDataDir, True)
        else:
            log_wr(LF, "Err: directory '%s' already exists." % outDataDir)
            log_fail(LF, "Use -o option to overwrite.")
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
                               freqArr_Hz,
                               sumBox_pix=int(pDict["sumBoxPix"]),
                               polyOrd=int(pDict["polyOrd"]),
                               outDataDir=outDataDir,
                               LF=LF)
    
    # Write the relevant tables to the database
    insert_arr_db(cursor, specRec, "spectraParms")
    conn.commit()

    # Close the connection to the database
    cursor.close()
    conn.close()
    log_wr(LF, "Database updated with spectral parameters.")


#-----------------------------------------------------------------------------#
if __name__=="__main__":
    main()

