#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     6_measure_complexity.py                                           #
#                                                                             #
# USAGE:    ./6_measure_complexity.py PATH/TO/SESSION                         #
#                                                                             #
# PURPOSE:  Create residual q-u spectra after subtracting a thin component    #
#           and measure the Faraday complexity.                               #
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

from Imports.module_measure_complexity import mod_measure_complexity

# Constants
C = 2.99792458e8

# Turn off print statements buffering
sys.stdout = os.fdopen(sys.stdout.fileno(), "w", 0)

# Map numpy data-types to the limited sqlite3 data-types
register_sqlite3_numpy_dtypes()

#-----------------------------------------------------------------------------#
def main():
    """
    Start the measure_complexity function from the command line.
    """
    
    # Help string to be shown using the -h option
    descStr = """
    This script loops through the catalogue in the database and subtract a
    Faraday thin model from the Q and U data if a polarised peak was
    detected in the Faraday dispersion function (FDF). The residual q=Q/I and
    u=U/I data are examined for evidence of additional components, indicating
    Faraday complexity. If a clean-component FDF has been saved on disk then a
    moment analysis is performed as an alternative complexity measurement.

    The model Stokes Q and U data are saved in the 'PATH/TO/SESSION/OUT'
    directory in FITS files named for the source in the catalogue. Measurements
    are saved to the SQLite database in the file 'session.sqlite'.

    Example:

    ./6_measure_complexity.py testSession/
    """
    
    # Parse the command line options
    parser = argparse.ArgumentParser(description=descStr,
                                 formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("sessionPath", metavar="PATH/TO/SESSION", nargs=1,
                        help="Path to the new session directory [no default]")
    args = parser.parse_args()
    sessionPath = args.sessionPath[0]
    sessionPath = sessionPath.rstrip("/")
    
    # Check the required directory structure exists or exit
    sessionRootDir, sessionName = os.path.split(sessionPath)
    sessionRootDir = "." if sessionRootDir=="" else sessionRootDir
    fail_not_exists(sessionPath, "directory")
    
    # Open the existing logfile to append to
    logFile = sessionPath + "/pipeline.log"
    LF = open(logFile, "a", 0)

    # Check that preceeding steps have been done
    statusFile = sessionPath + "/status.json"
    fail_not_exists(statusFile, "file", LF)
    statusDict = read_dictfile(statusFile)
    if int(statusDict["session"])<1:
        log_fail(LF, "Err: Session status file reports session was not " + \
                 "created successfully.")
    if int(statusDict["extract"])<1:
        log_fail(LF, "Err: Session status file reports spectral " + \
                 " extraction was not done.")
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

    # Load the catalogue from the database
    sql = """
    SELECT spectraParms.uniqueName
    FROM spectraParms
    """
    log_wr(LF, "Loading catalogue from the database.")
    catRec = select_into_arr(cursor, sql)
    nRows = int(len(catRec))
    if nRows==0:
        log_fail(LF, "Err: query returned zero entries. Exiting.")
        log_wr(LF, "%s rows returned." % nRows)
        
    
    # RUN THE COMPLEXITY ASSESSMENT MODULE -----------------------------------#
    compRec = mod_measure_complexity(catRec,
                                     sessionPath,
                                     LF=LF)
    
    # Write the complexity measurements to the database
    cursor.execute("DELETE FROM complexMeasures")
    insert_arr_db(cursor, compRec, "complexMeasures")
    conn.commit()
    log_wr(LF, "Database updated with complexity measurements.")
    
    # Close the connection to the database
    cursor.close()
    conn.close()
    
    # Update the status file to reflect successful measurement
    statusDict["complexity"] = 1
    write_dictfile(statusDict, sessionPath + "/status.json")

    
#-----------------------------------------------------------------------------#
if __name__=="__main__":
    main()

    
