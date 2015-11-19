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
# MODIFIED: 19-November-2015 by C. Purcell                                    #
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
import argparse
import time
import traceback
import math as m
import numpy as np
import sqlite3

from Imports.util_PPC import fail_not_exists
from Imports.util_PPC import log_wr
from Imports.util_PPC import log_fail
from Imports.util_PPC import read_dictfile
from Imports.util_PPC import write_dictfile

from Imports.util_DB import register_sqlite3_numpy_dtypes
from Imports.util_DB import select_into_arr
from Imports.util_DB import insert_arr_db
from Imports.util_DB import update_arr_db

from Imports.module_RM_clean import mod_do_RMclean
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
    Start the run_RM_clean if called from the command line.
    """
    
    # Help string to be shown using the -h option
    descStr = """
    Loop through the catalogue in the database and perform RM-CLEAN on the
    Faraday dispersion functions that have previously been calculated. Measure
    the properties of the clean Faraday dispersion function and save to a new
    table in the database.

    Each clean Faraday dispersion function and clean-component spectrum is
    saved in ASCII format to a directory called 'PATH/TO/SESSION/OUT'. The
    ASCII files are named for the position of the source on the sky.

    Note: All measurements on the clean Faraday dispersion function are saved
    to the SQLite database in the file 'PATH/TO/SESSION/session.sqlite'.

    Example:

    ./5_do_RM-clean.py testSession/
    """

    # Parse the command line options
    parser = argparse.ArgumentParser(description=descStr,
                                 formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("sessionPath", metavar="PATH/TO/SESSION", nargs=1,
                        help="Path to the new session directory [no default]")
    parser.add_argument("-o", dest="doOverwrite", action="store_true",
                        help="Overwrite old RM-clean results")
    args = parser.parse_args()
    sessionPath = args.sessionPath[0]
    sessionPath = sessionPath.rstrip("/")
    doOverwrite = args.doOverwrite
    
    # Check the required directory structure exists or exit
    fail_not_exists(sessionPath, "directory")
    
    # Open the existing logfile to append to
    logFile = sessionPath + "/pipeline.log"
    LF = open(logFile, "a", 0)
    log_wr(LF, ">>> Beginning RM-clean.")
    log_wr(LF, "Time: %s" % time.ctime() )

    # Check that preceeding steps have been done
    statusFile = sessionPath + "/status.json"
    fail_not_exists(statusFile, "file", LF)
    statusDict = read_dictfile(statusFile)
    if int(statusDict["session"])<1:
        log_fail(LF, "Err: Session status file reports session was not " + \
                 "created successfully.")
    if int(statusDict["extract"])<1:
        log_fail(LF, "Err: Session status file reports spectral " + \
                 "extraction was not done.")
    if int(statusDict["rmsynth"])<1:
        log_fail(LF, "Err: Session status file reports RM-synthesis " + \
                 "was not done.")

    # Connect to the database
    dbFile = sessionPath + "/session.sqlite"
    try:
        log_wr(LF, "> Connecting to existing DB file '%s' ..." % dbFile)
        conn = sqlite3.connect(dbFile)
        cursor = conn.cursor()
    except Exception:
        log_wr(LF, "Err: Failed to connect to '%s'." % dbFile)
        log_fail(LF, traceback.format_exc())
        
    # Load necessary columns from the spectral parameter table into a recarray
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
                              sessionPath,
                              doOverwrite=doOverwrite,
                              LF=LF)
    
    # Write the RM-clean parameters to the database
    # [uniqueName, nIterDone, cleanCutoff_sigma, cleanCutoff_Jybm]
    cursor.execute("DELETE FROM cleanFDFparms")
    insert_arr_db(cursor, cleanRec, "cleanFDFparms")
    conn.commit()
    log_wr(LF, "Database updated with CLEANING parameters.")
    
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
                             sessionPath,
                             dirty=False,
                             LF=LF)

    # Update the table with the FDF measurements
    update_arr_db(cursor, fdfRec, "cleanFDFparms", "uniqueName")
    conn.commit()
    log_wr(LF, "Database updated with FDF parameters.")
    
    # Close the connection to the database
    cursor.close()
    conn.close()
    
    # Update the status file to reflect successful extraction
    statusDict["rmclean"] = 1
    write_dictfile(statusDict, sessionPath + "/status.json")

#-----------------------------------------------------------------------------#
if __name__=="__main__":
    main()
