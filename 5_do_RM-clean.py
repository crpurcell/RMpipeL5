#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     5_do_RM-clean.py                                                  #
#                                                                             #
# USAGE:    ./5_do_RM-clean.py PATH/TO/SESSION                                #
#                                                                             #
# PURPOSE:  Perform RM-clean on the Faraday dispersion Functions in a         #
#           POSSUM pipeline session.                                          #
#                                                                             #
# MODIFIED: 21-Apr-2015 by C. Purcell                                         #
#                                                                             #
#=============================================================================#

import os
import sys
import shutil
import argparse
import math as m
import numpy as np
import sqlite3
from scipy.stats.stats import nanmedian

from Imports.util_PPC import *
from Imports.util_DB import *
from Imports.module_RM_clean import *
from Imports.module_measure_fdf import *

# Constants
C = 2.99792458e8

# Turn off print statements buffering
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

# Map numpy data-types to the limited sqlite3 data-types
register_sqlite3_numpy_dtypes()


#-----------------------------------------------------------------------------#
def main():
    """
    Start the run_RM_clean if called from the command line.
    """
    
    # Help string to be shown using the -h option
    descStr = """
    Loop through the catalogue and perform RM-synthesis on the Stokes I, Q
    and U spectra that have previously been extracted.
    """

    # Parse the command line options
    parser = argparse.ArgumentParser(description=descStr)
    parser.add_argument('sessionPath', metavar='PATH/TO/SESSION', nargs=1,
                        help='Path to the new session directory [no default]')
    parser.add_argument('-o', dest='doOverwrite', action='store_true',
                        help='Overwrite old RM-clean results')
    args = parser.parse_args()
    sessionPath = args.sessionPath[0]
    doOverwrite = args.doOverwrite

    # Call the RM-synthesis function
    run_RM_clean(sessionPath, doOverwrite)


#-----------------------------------------------------------------------------#
def run_RM_clean(sessionPath, doOverwrite=False):

    sessionPath = sessionPath.rstrip('/')
    
    # Check the required directory structure exists or exit
    sessionRootDir, sessionName = os.path.split(sessionPath)
    sessionRootDir = '.' if sessionRootDir=='' else sessionRootDir
    fail_not_exists(sessionPath, 'directory')
    
    # Set the session status file to '1'
    # TODO: set this in the database
    statusFile = sessionPath + '/status.txt'
    set_statusfile(statusFile, 1)
    
    # Open the existing logfile to append to
    logFile = sessionPath + '/pipeline.log'
    LF = open(logFile, 'a', 0)
    log_wr(LF, ">>> Beginning RM-clean.")

    # Read and parse the pipeline input file
    inParmFile = sessionPath + '/inputs.config'
    fail_not_exists(inParmFile, 'file', LF)
    try:
        pipeInpObj = PipelineInputs(inParmFile)
        pDict = pipeInpObj.get_flat_dict()
        log_wr(LF, 'Successfully parsed the input parameter file.')
    except Exception:
        log_fail(LF, 'Err: Failed to parse the input parameter file.')

    # Verify that the parameter file has all the correct entries
    missingLst = pipeInpObj.inparm_verify()
    if len(missingLst)>0:
        log_fail(LF, 'Err: Required input parameters missing. - %s' %
                 missingLst)
    else:
        log_wr(LF, 'All required input parameters present.')

    # Check that the input/output data directory exists
    specPath = sessionPath + '/OUT'
    fail_not_exists(specPath, 'directory', LF)

    # Check that the specified dataset exists
    dataPath = pDict['dataPath'].rstrip('/')
    fail_not_exists(dataPath, 'directory', LF)

    # Get the frequency vector from the text file
    inFreqFile = dataPath + '/freqs_Hz.txt'
    fail_not_exists(inFreqFile, 'file', LF)
    freqArr_Hz = load_vector_fail(inFreqFile, "float32", False, LF)
    lamSqArr_m2 = np.power(C/freqArr_Hz, 2.0)
    
    # Connect to the database
    dbFile = sessionPath + '/session.sqlite'
    try:
        log_wr(LF, "> Connecting to existing DB file '%s' ..." % dbFile)
        conn = sqlite3.connect(dbFile)
        cursor = conn.cursor()
    except Exception:
        log_fail(LF, "Err: Failed to connect to '%s'." % dbFile)
        
    # Load the parameters into memory
    sql = """
    SELECT uniqueName, rmsMedQUAvg_Jybm
    FROM spectraParms
    """
    log_wr(LF, "Loading the spectral parameters from the database.")
    specRec  = select_into_arr(cursor, sql)
    nRows = int(len(specRec))
    if nRows==0:
        log_fail(LF, "Err: query returned zero entries. Exiting.")
    log_wr(LF, "%s rows returned." % nRows)

    # RUN THE RM-CLEAN MODULE ------------------------------------------------#
    cleanRec = mod_do_RMclean(specRec,
                              specPath,
                              specPath,
                              float(pDict['cleanCutoff_sigma']),
                              int(pDict['maxCleanIter']),
                              float(pDict['gain']),
                              doOverwrite=doOverwrite,
                              LF=LF)
    
    # Write the RMSF parameters to the database
    cursor.execute("DELETE FROM cleanFDFparms")
    insert_arr_db(cursor, cleanRec, 'cleanFDFparms')
    conn.commit()
    log_wr(LF, 'Database updated with CLEANING parameters.')
    
    # Load the spectral & RMSF parameters into memory
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
                             specPath,
                             lamSqArr_m2,
                             float(pDict['thresholdSignalPI_sigma']),
                             fileSuffix='_cleanFDF.dat',
                             LF=LF)

    # Update the table with the FDF measurements
    update_arr_db(cursor, fdfRec, 'cleanFDFparms', 'uniqueName')
    conn.commit()
    log_wr(LF, 'Database updated with FDF parameters.')
    
    # Close the connection to the database
    cursor.close()
    conn.close()
    

#-----------------------------------------------------------------------------#
if __name__=="__main__":
    main()
