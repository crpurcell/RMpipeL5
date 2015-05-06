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
# MODIFIED: 20-Apr-2015 by C. Purcell                                         #
#                                                                             #
#=============================================================================#

import os
import sys
import shutil
import argparse
import time
import math as m
import numpy as np
import sqlite3
#from scipy.stats.stats import nanmedian

from Imports.util_PPC import *
from Imports.util_DB import *
from Imports.module_RM_synthesis import *
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
    Start the run_RM_synthesis function if called from the command line.
    """
    
    # Help string to be shown using the -h option
    descStr = """
    Loop through the catalogue in the database and perform RM-synthesis on the
    Stokes I, Q and U spectra that have previously been extracted. Measure the
    properties of the resultant Faraday dispersion function.
    """

    # Parse the command line options
    parser = argparse.ArgumentParser(description=descStr)
    parser.add_argument('sessionPath', metavar='PATH/TO/SESSION', nargs=1,
                        help='Path to the new session directory [no default]')
    parser.add_argument('-o', dest='doOverwrite', action='store_true',
                        help='Overwrite old RM-synthesis results')
    args = parser.parse_args()
    sessionPath = args.sessionPath[0]
    doOverwrite = args.doOverwrite

    # Call the RM-synthesis function
    run_RM_synthesis(sessionPath, doOverwrite)


#-----------------------------------------------------------------------------#
def run_RM_synthesis(sessionPath, doOverwrite=False):
    
    sessionPath = sessionPath.rstrip('/')
    
    # Check the required directory structure exists or exit
    fail_not_exists(sessionPath, 'directory')
    
    # Open the existing logfile to append to
    logFile = sessionPath + '/pipeline.log'
    LF = open(logFile, 'a', 0)
    log_wr(LF, "\n>>> Beginning RM-synthesis.")
    log_wr(LF, "Time: %s" % time.ctime() )

    # Read and parse the pipeline input file
    inParmFile = sessionPath + '/inputs.config'
    fail_not_exists(inParmFile, 'file', LF)
    try:
        pipeInpObj = PipelineInputs(inParmFile)
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

    # Calculate the Faraday depth sampling  etc
    log_wr(LF, 'Calculating the Faraday depth sampling')
    pipeInpObj.calculate_derived_parms(resetPhiSamp=False)
    pDict = pipeInpObj.get_flat_dict(includeDerived=True)
    phiArr = pDict['phiArr_radm2']
    log_wr(LF, '> PhiArr = %.2f to %.2f by %.2f' % (phiArr[0], phiArr[-1],
                                                   float(pDict['dPhi_radm2'])))
    
        
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
    
    # Calculate the RM sampling from the input parameters
    # Connect to the database
    dbFile = sessionPath + '/session.sqlite'
    try:
        log_wr(LF, "> Connecting to existing DB file '%s' ..." % dbFile)
        conn = sqlite3.connect(dbFile)
        cursor = conn.cursor()
    except Exception:
        log_fail(LF, "Err: Failed to connect to '%s'." % dbFile)
        
    # Load the spectral parameters into memory    
    sql = """
    SELECT uniqueName, coeffPolyIspec, rmsMedQUAvg_Jybm
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
                             specPath,
                             specPath,
                             phiArr,
                             weightType=pDict['weightType'],
                             doOverwrite=doOverwrite,
                             LF=LF)
    
    # Write the RMSF parameters to the database
    cursor.execute("DELETE FROM dirtyFDFparms")
    insert_arr_db(cursor, rmsfRec, 'dirtyFDFparms')
    conn.commit()
    log_wr(LF, 'Database updated with RMSF parameters.')

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
                             fileSuffix='_dirtyFDF.dat',
                             LF=LF)

    # Update the table with the FDF measurements
    update_arr_db(cursor, fdfRec, 'dirtyFDFparms', 'uniqueName')
    conn.commit()
    log_wr(LF, 'Database updated with FDF parameters.')
    
    # Close the connection to the database
    cursor.close()
    conn.close()
    


#-----------------------------------------------------------------------------#
if __name__=="__main__":
    main()

    
    
