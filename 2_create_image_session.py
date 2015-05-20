#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     create_image_session.py                                           #
#                                                                             #
# USAGE:    create_image_session.py PATH/TO/DATA  PATH/TO/SESSION CATFILE     #
#                                                        [CATFORMATFILE]      #
#                                                                             #
# PURPOSE:  Create a new pipeline processing session in a new directory.      #
#           A dataset and source catalogue must be chosen at this time and    #
#           linked to the session. The data must have been processed using    #
#           the 'verify_image_data.py' script. The ASCII catalogue file is    #
#           described by a SQL statement in the optional CATFORMATFILE and    #
#           defaults to the format produced by the Aegean source-finder.      #
#           This script performs some sanity-checks and creates a pipeline    #
#           input file with sensible input parameters determined from the     #
#           properties of the data.                                           #
#                                                                             #
# MODIFIED: 20-May-2015 by C. Purcell                                         #
#                                                                             #
#=============================================================================#

import os
import sys
import shutil
import traceback
import argparse
import time
import sqlite3
import ConfigParser
import string
import math as m
import numpy as np

from Imports.util_PPC import PipelineInputs
from Imports.util_PPC import fail_not_exists
from Imports.util_PPC import log_wr
from Imports.util_PPC import log_fail
from Imports.util_PPC import load_vector_fail
from Imports.util_PPC import deg2dms
from Imports.util_PPC import cat_to_recarray

from Imports.util_DB import register_sqlite3_numpy_dtypes
from Imports.util_DB import create_db
from Imports.util_DB import insert_arr_db
from Imports.util_DB import schema_to_tabledef

# Constants
C = 2.99792458e8

# Turn off print statements buffering
sys.stdout = os.fdopen(sys.stdout.fileno(), "w", 0)

# Map numpy data-types to the limited sqlite3 data-types
register_sqlite3_numpy_dtypes()


#-----------------------------------------------------------------------------#
def main():
    """
    Start the create_image_session function if called from the command line.
    """
    
    # Help string to be shown using the -h option
    descStr = """
    Create a new pipeline processing session under the PATH/TO/SESSION
    directory and link to the data in the PATH/TO/DATA directory. 'SESSION' is
    the session name and this directory is created if necessary or overwritten
    if it already exists (using the -o flag). The data (in FITS format) must
    exist and have been verified using the 'verify_image_data.py' script. An
    ASCII source catalogue file must be provided. The catalogue format defaults
    to that produced by the Aegean source-finder. Alternative formats may be
    specified by a SQL statement in the CATFORMATFILE, given as a command line
    argument (see the example file in Imports/templates/catDescDefault.sql).
    This script creates a pipeline input file 'PATH/TO/SESSION/inputs.config',
    populated with sensible defaults based on the parameters of the data.

    Note: The input catalogue and all results are saved to a SQLite database
    in the file 'PATH/TO/SESSION/session.sqlite'.
    
    Example:

    ./2_create_image_session.py testData/ testSession/ testData/testCat.dat
                                testData/testCatDesc.sql
    """
    
    # Parse the command line options
    parser = argparse.ArgumentParser(description=descStr,
                                 formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-o", dest="doOverwrite", action="store_true",
                        help="Overwrite (delete) existing session")
    parser.add_argument("dataPath", metavar="PATH/TO/DATA", nargs=1,
                        help="Path to data directory [no default]")
    parser.add_argument("sessionPath", metavar="PATH/TO/SESSION", nargs=1,
                        help="Path to the new session directory [no default]")
    parser.add_argument("catPath", metavar="CATFILE", nargs=1,
                        help="Name of the ASCII catalogue file [no default]")
    parser.add_argument("catFormatPath", metavar="CATFORMATFILE", nargs="?",
                        default="Imports/templates/catDescDefault.sql",
                  help="Name of the file describing the catalogue (optional)")
    args = parser.parse_args()
    doOverwrite = args.doOverwrite
    dataPath = args.dataPath[0]
    sessionPath = args.sessionPath[0]
    catPath = args.catPath[0]
    catFormatPath = args.catFormatPath

    # Call the create_session function
    create_image_session(dataPath, sessionPath, catPath, catFormatPath,
                         doOverwrite)


#-----------------------------------------------------------------------------#
def create_image_session(dataPath, sessionPath, catPath, catFormatPath,
                         doOverwrite=False):
    """
    Create an image session directory, link to a dataset and catalogue and
    populate a default pipeline driving file with sensible values.
    """
    
    dataPath = dataPath.rstrip("/")
    sessionPath = sessionPath.rstrip("/")
    catPath = catPath.rstrip("/")
    pDict = {}
        
    # Create a new session directory. Erase, or exit if it already exists
    sessionRootDir, sessionName = os.path.split(sessionPath)
    sessionRootDir = "." if sessionRootDir=="" else sessionRootDir
    if os.path.exists(sessionPath):
        if doOverwrite:
            print "Deleting exiting session '%s'."  % sessionPath
            shutil.rmtree(sessionPath)
        else:
            print "Err: directory already exists '%s'." % sessionPath
            print "Use -o option to overwrite."
            sys.exit()
    print "Creating new session directory '%s'." % sessionPath
    dirs = sessionPath.split("/")
    for i in range(1, len(dirs)+1):
        dirStr = "/".join(dirs[:i])
        if not os.path.exists(dirStr):
            os.mkdir(dirStr)

    # Remove the old logfile and open a new one
    logFile = sessionPath + "/pipeline.log"
    LF = open(logFile, "w", 0)
    try:
        log_wr(LF, "Logging messages for session to '%s'." % logFile)
        log_wr(LF, "Start time: %s" % time.ctime() )
    except Exception:
        print "Err: Failed to open the log file '%s'." % logFile
        sys.exit()
    
    # Check other inputs exist
    dataPath = "." if dataPath=="" else dataPath
    fail_not_exists(dataPath, "directory")
    pDict["dataPath"] = dataPath
    catRootDir, catFile = os.path.split(catPath)
    fail_not_exists(catPath)
    catFormatRootDir, catFormatFile = os.path.split(catFormatPath)
    fail_not_exists(catFormatPath)

    # Check that the catalogue description file has X and Y columns
    log_wr(LF, "Parsing catalogue description file '%s'" % catFormatPath)
    catDefDict, catSQLdict = schema_to_tabledef(catFormatPath,
              addColDict={"sourceCat":"uniqueName varchar(20) PRIMARY KEY"})
    catHasX = False
    catHasY = False
    countUniqueName = 0
    for e in catDefDict["sourceCat"]:
        if e[0]=="x_deg":
            catHasX = True
        if e[0]=="y_deg":
            catHasY = True
        if e[0]=="uniqueName":
            countUniqueName += 1
    if catHasX is False:
        log_fail(LF, "Catalogue description missing 'x_deg' column.")
    if catHasY is False:
        log_fail(LF, "Catalogue description missing 'y_deg' column.")
    if countUniqueName>1:
        log_fail(LF, "Catalogue description contains reserved" + \
                     "'uniqueName' column. Please rename in description.")
    log_wr(LF, "Catalogue description contains 'x_deg' and 'y_deg' columns.")

    # Read the input catalogue to a record array
    log_wr(LF, "Reading the catalogue into memory ...")
    try:
        catRec = cat_to_recarray(catPath, catDefDict["sourceCat"],
                                 delim=" ", LF=LF)
    except Exception:
        log_wr(LF, "Failed to parse or read catalogue.")
        log_fail(LF, traceback.format_exc())

    # Populate the unique-name field
    log_wr(LF, "Assigning a unique name to each entry (position based).")
    for i in range(len(catRec)):
        catRec[i]["uniqueName"] = \
             deg2dms(catRec[i]["x_deg"]/15.0, delim="", nPlaces=1) + \
             deg2dms(catRec[i]["y_deg"], delim="", doSign=True, nPlaces=2)

    # Read the SQL schema for the internal database
    schemaFile = "Imports/templates/DBSchema.sql"
    log_wr(LF, "Parsing the database schema from '%s'" % schemaFile)
    tableDefDict, tableSQLdict = schema_to_tabledef(schemaFile)
    
    # Create a new persistent database using the table schema
    dbFile = sessionPath + "/session.sqlite"
    success = create_db(dbFile, dict(catSQLdict.items()+tableSQLdict.items()))
    if success:
        log_wr(LF, "Successfully created the database '%s'." % dbFile)
    else:
        log_fail(LF, "Err: failed to create the database '%s'." % dbFile)

    # Insert the catalogue into the database file
    log_wr(LF, "Inserting the catalogue into the database ...")
    conn = sqlite3.connect(dbFile)
    cursor = conn.cursor()
    try:
        insert_arr_db(cursor, catRec, "sourceCat")
        conn.commit()
    except Exception:
        log_wr(LF, "Err: Failed to insert the catalogue into the database.")
        log_fail(LF, traceback.format_exc())    
    cursor.close()
    conn.close()
    
    # Check for the frequency sampling file in the data directory
    inFreqFile = dataPath + "/freqs_Hz.txt"
    if os.path.exists(inFreqFile):
        pDict["inFreqFile"] = "freqs_Hz.txt"
    else:
        log_fail(LF, "File containing the frequency vector is missing: '%s'." \
                 % inFreqFile)
        
    # Check for the datatype file in the data directory
    inTypeFile = dataPath + "/dataType.txt"
    if os.path.exists(inTypeFile):
        pDict["inTypeFile "] = "dataType.txt"
    else:
        log_fail(LF, "File containing the data type is missing: '%s'." \
                 % inTypeFile)

    # Create a new pipeline input object. Default parameters are taken from
    # the file 'Imports/templates/defaultInputs.config'
    pipeInpObj = PipelineInputs("Imports/templates/defaultInputs.config")

    # Set the dataset in the input object and calculate the defaults
    log_wr(LF, "Calculating pipeline inputs for linked dataset.")
    try:        
        pipeInpObj.config.set("Dataset", "dataPath", dataPath)
        pipeInpObj.calculate_derived_parms(resetPhiSamp=True)
        pD = pipeInpObj.get_flat_dict(includeDerived=True)    
    except Exception:
        log_wr(LF, "Err: Failed to calculate the Faraday depth sampling.")
        log_fail(LF, traceback.format_exc())

    # Feedback
    log_wr(LF, "> min(dLambda-squared) = %.3g m^2" % \
           float(pD["dLambdaSqMin_m2"]))
    log_wr(LF, "> max(phi) = %.1f rad m^-2" % float(pD["stopPhi_radm2"]))
    log_wr(LF, "> dPhi = %.1f" % float(pD["dPhi_radm2"]))
    log_wr(LF, "> nChan(phi) = %d" % int(pD["nChanRM"]))
        
    # Point to the input file for the current session and write to disk
    pipeInFile = sessionPath + "/inputs.config"
    log_wr(LF, "Writing pipeline driving file to '%s'" % pipeInFile)
    pipeInpObj.configFile = pipeInFile
    pipeInpObj.write_file()
    log_wr(LF, "> Edit this file to make changes to input parameters.")
    
    # Copy the catalogue and format file to the session directory
    log_wr(LF, "Copying dataset metadata files to the session directory:")
    log_wr(LF, "> '%s' -> '%s'" % (catPath, sessionPath + "/" + catFile))
    shutil.copyfile(catPath, sessionPath + "/" + catFile)
    log_wr(LF, "> '%s' -> '%s'" % (catFormatPath,
                                   sessionPath + "/" + catFormatFile)) 
    shutil.copyfile(catFormatPath, sessionPath + "/" + catFormatFile)
    log_wr(LF, "> '%s' -> '%s'" % (inFreqFile, sessionPath + "/freqs_Hz.txt"))
    shutil.copyfile(inFreqFile, sessionPath + "/freqs_Hz.txt")
    log_wr(LF, "> '%s' -> '%s'" % (inTypeFile, sessionPath + "/dataType.txt"))
    shutil.copyfile(inTypeFile, sessionPath + "/dataType.txt")


#-----------------------------------------------------------------------------#
if __name__=="__main__":
    main()

