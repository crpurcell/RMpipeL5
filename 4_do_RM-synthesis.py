#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     4_do_RM-synthesis.py                                              #
#                                                                             #
# USAGE:    ./4_do_RM-synthesis.py PATH/TO/SESSION                            #
#                                                                             #
# PURPOSE:  Perform RM-synthesis on the extracted spectra in the current      #
#           pipeline session.                                                 #
#                                                                             #
# MODIFIED: 29-September-2015 by C. Purcell                                   #
#                                                                             #
#=============================================================================#

import os
import sys
import shutil
import argparse
import time
import traceback
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
from Imports.util_DB import update_arr_db

from Imports.module_RM_synthesis import mod_do_RMsynth
from Imports.module_measure_FDF import mod_measure_FDF

# Constants
C = 2.99792458e8

# Turn off print statements buffering
sys.stdout = os.fdopen(sys.stdout.fileno(), "w", 0)

# Map numpy data-types to the limited sqlite3 data-types
register_sqlite3_numpy_dtypes()


#-----------------------------------------------------------------------------#
def main():
    """
    Start the run_RM_synthesis function if called from the command line.
    """
    
    # Help string to be shown using the -h option
    descStr = """
    Loop through the catalogue in the database and perform RM-synthesis on the
    Stokes I, Q and U spectra that have previously been extracted. Measure the
    properties of the resultant Faraday dispersion function and save to a new
    table in the database.

    Each complex Faraday dispersion function is saved in FITS format to
    a directory called 'PATH/TO/SESSION/OUT'. The FITS files are named for
    the source in the catalogue.

    Note: All measurements on the dirty Faraday dispersion function are saved
    to the SQLite database in the file 'PATH/TO/SESSION/session.sqlite'.

    Example:

    ./4_do_RM-synthesis.py testSession/
    """

    # Parse the command line options
    parser = argparse.ArgumentParser(description=descStr,
                                 formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("sessionPath", metavar="PATH/TO/SESSION", nargs=1,
                        help="Path to the new session directory [no default]")
    parser.add_argument("-o", dest="doOverwrite", action="store_true",
                        help="Overwrite old RM-synthesis results")
    args = parser.parse_args()
    sessionPath = args.sessionPath[0]
    doOverwrite = args.doOverwrite

    # Call the RM-synthesis function
    run_RM_synthesis(sessionPath, doOverwrite)


#-----------------------------------------------------------------------------#
def run_RM_synthesis(sessionPath, doOverwrite=False):
    
    sessionPath = sessionPath.rstrip("/")
    
    # Check the required directory structure exists or exit
    fail_not_exists(sessionPath, "directory")
    
    # Open the existing logfile to append to
    logFile = sessionPath + "/pipeline.log"
    LF = open(logFile, "a", 0)
    log_wr(LF, "\n>>> Beginning RM-synthesis.")
    log_wr(LF, "Time: %s" % time.ctime() )

    # Check that preceeding steps have been done
    statusFile = sessionPath + "/status.json"
    fail_not_exists(statusFile, "file", LF)
    statusDict = read_dictfile(statusFile)
    if int(statusDict["session"])<1:
        log_fail(LF, "Err: Session status file reports session was not " + \
                     "created successfully.")
    if int(statusDict["extract"])<1:
        log_fail(LF, "Err: Session status file reports spectral extraction " + \
                     "was not done.")

    # Read and parse the pipeline input file
    inParmFile = sessionPath + "/inputs.config"
    fail_not_exists(inParmFile, "file", LF)
    try:
        pipeInpObj = PipelineInputs(inParmFile)
        log_wr(LF, "Successfully parsed the input parameter file.")
    except Exception:
        log_wr(LF, "Err: Failed to parse the input parameter file.")
        log_fail(LF, traceback.format_exc())

    # Verify that the parameter file has all the correct entries
    missingLst = pipeInpObj.inparm_verify()
    if len(missingLst)>0:
        log_fail(LF, "Err: Required input parameters missing. - %s" %
                 missingLst)
    else:
        log_wr(LF, "All required input parameters present.")

    # Calculate the Faraday depth sampling etc.
    log_wr(LF, "Calculating the Faraday depth sampling")
    pipeInpObj.calculate_derived_parms(resetPhiSamp=False)
    pDict = pipeInpObj.get_flat_dict(includeDerived=True)
    phiArr = pDict["phiArr_radm2"]
    log_wr(LF, "> PhiArr = %.2f to %.2f by %.2f" % (phiArr[0], phiArr[-1],
                                                   float(pDict["dPhi_radm2"])))
        
    # Check that the input/output data directory exists
    specPath = sessionPath + "/OUT"
    fail_not_exists(specPath, "directory", LF)
    
    # Check that the specified dataset exists
    dataPath = pDict["dataPath"].rstrip("/")
    fail_not_exists(dataPath, "directory", LF)
    
    # Get the frequency vector from the text file
    #inFreqFile = dataPath + "/freqs_Hz.txt"
    #fail_not_exists(inFreqFile, "file", LF)
    #freqArr_Hz = load_vector_fail(inFreqFile, "float32", False, LF)
    #lamSqArr_m2 = np.power(C/freqArr_Hz, 2.0)
    lamSqArr_m2 = pDict["lambdaSqArr_m2"]
    
    # Connect to the database
    dbFile = sessionPath + "/session.sqlite"
    try:
        log_wr(LF, "> Connecting to existing DB file '%s' ..." % dbFile)
        conn = sqlite3.connect(dbFile)
        cursor = conn.cursor()
    except Exception:
        log_wr(LF, "Err: Failed to connect to '%s'." % dbFile)
        log_fail(LF, traceback.format_exc())
        
    # Load the spectral parameter table into memory
    sql = """
    SELECT uniqueName, coeffPolyIspec, rmsMedQUAvg_Jybm, extractStatus
    FROM spectraParms
    """
    log_wr(LF, "Loading the spectral parameters from the database.")
    specRec = select_into_arr(cursor, sql)
    nRows = int(len(specRec))
    if nRows==0:
        log_fail(LF, "Err: query returned zero entries. Exiting.")
    log_wr(LF, "%s rows returned." % nRows)
    
    # RUN THE RM-SYNTHESIS MODULE --------------------------------------------#
    rmsfRec = mod_do_RMsynth(specRec,
                             sessionPath,
                             phiArr,
                             weightType=pDict["weightType"],
                             doOverwrite=doOverwrite,
                             LF=LF)

    # Write the RMSF parameters to the database
    cursor.execute("DELETE FROM dirtyFDFparms")
    insert_arr_db(cursor, rmsfRec, "dirtyFDFparms")
    conn.commit()
    log_wr(LF, "Database updated with RMSF parameters.")
    
    msg = "\n" + "-"*80 + "\n"
    msg += "Proceeding to measure the properties of the FDFs."
    msg = "\n" + "-"*80 + "\n"
    log_wr(LF, msg)
    
    # Load the spectral & RMSF parameters back into memory
    sql = """
    SELECT spectraParms.uniqueName,
    spectraParms.coeffPolyIspec,
    spectraParms.rmsMedQUAvg_Jybm,
    dirtyFDFparms.fwhmRMSF,
    dirtyFDFparms.lam0Sq_m2
    FROM spectraParms LEFT JOIN dirtyFDFparms
    ON spectraParms.uniqueName = dirtyFDFparms.uniqueName
    """
    log_wr(LF, "Reloading the spectral & RMSF parameters from the database.")
    catRec = select_into_arr(cursor, sql)
    nRows = int(len(catRec))
    if nRows==0:
        log_fail(LF, "Err: query returned zero entries. Exiting.")
    log_wr(LF, "%s rows returned." % nRows)
    
    # RUN THE FDF MEASUREMENT MODULE------------------------------------------#
    fdfRec = mod_measure_FDF(catRec,
                             sessionPath,
                             float(pDict["thresholdSignalPI_sigma"]),
                             dirty=True,
                             LF=LF)
    
    # Update the table with the FDF measurements
    update_arr_db(cursor, fdfRec, "dirtyFDFparms", "uniqueName")
    conn.commit()
    log_wr(LF, "Database updated with FDF parameters.")
    
    # Close the connection to the database
    cursor.close()
    conn.close()

    # Update the status file to reflect successful extraction
    statusDict["rmsynth"] = 1
    write_dictfile(statusDict, sessionPath + "/status.json")


#-----------------------------------------------------------------------------#
if __name__=="__main__":
    main()

    
    
