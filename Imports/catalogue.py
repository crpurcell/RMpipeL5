#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     catalogue.py                                                      #
#                                                                             #
# PURPOSE:  Definition of a class and methods to store the input catalogue    #
#           and tables of results for each stage of the POSSUM pipeline.      #
#                                                                             #
# REQUIRED:                                                                   #
#                                                                             #
# MODIFIED: 05-Feb-2014 by C. Purcell                                         #
#                                                                             #
# CONTENTS:                                                                   #
#                                                                             #
#=============================================================================#
import os
import sys
import shutil
import numpy as np
import sqlite3

from util_PPC import schema_to_tabledef
from util_PPC import cat_to_recarray
from util_PPC import deg2dms
from util_PPC import log_fail
from util_PPC import log_wr

from rec import *

#-----------------------------------------------------------------------------#
class tabStore:
    """
    Class to hold the input catalogues and results produced by the POSSUM
    pipeline.
    """
    
    def __init__(self, schemaFile, inCatFile, delim=" ", LF=None):
        """
        Create record arrays to store the results of the POSSUM pipeline.
        The record arrays are formatted by parsing the SQL CREATE statements
        in the database schema file. The input catalogue is parsed and inserted
        into the 'sourceCat' table, also using the CREATE statement as a guide.
        """
        
        # Tables are record arrays stored in a dictionary
        self.tbl = {}
        self.createSQL = {}
        self.nRows = 0
        
        # Parse the record array dtypes from the SQL schema
        # tableDefDict = {'tableName': [('colName1', 'dtype1'),  ... ]}
        tableDefDict, self.createSQL = schema_to_tabledef(schemaFile)

        # Read the input catalogue to a record array
        catRec = cat_to_recarray(inCatFile, tableDefDict['sourceCat'],
                                 delim=" ", LF=LF)
        self.tbl['sourceCat'] = catRec
        self.nRows = len(catRec)

        # Populate the unique-name field
        for i in range(self.nRows): 
            catRec[i]['uniqueName'] = \
                deg2dms(catRec[i]['x_deg']/15.0, delim='', nPlaces=1) + \
                deg2dms(catRec[i]['y_deg'], delim='', doSign=True, nPlaces=2)

        # Create empty recarray tables defined by the SQL schema
        for tabName, tabDef in tableDefDict.iteritems():
            if tabName=='sourceCat':
                continue
            self.tbl[tabName] = np.zeros(catRec.shape,
                                         dtype=tableDefDict[tabName])
            self.tbl[tabName]['uniqueName'] = \
                                            self.tbl['sourceCat']['uniqueName']
            
    def __call__(self, tabName):
        """
        Return the table if the class is called with a tableName argument.
        """
        if tabName in self.tbl:
            return self.tbl[tabName]
        else:
            return None

    def create_DB(self, dbFile, LF=None):
        """
        Create an empty SQLite3 database on disk.
        """
        if os.path.exists(dbFile):
            os.remove(dbFile)
        try:
            conn = sqlite3.connect(dbFile)
            cursor = conn.cursor()
            for sql in self.createSQL.values():
                cursor.execute(sql)    
            conn.commit()
            cursor.close()
            conn.close()
            log_wr(LF, "Successfully created the database'%s'." % dbFile)
        except Exception:
            log_fail(LF, "Err: failed to create the database '%s'." % dbFile)
        
    def write_to_DB(self, dbFile, tabName):
        """
        Write a table into the sqlite3 database on disk. 
        """

        conn = sqlite3.connect(dbFile)
        cursor = conn.cursor()
        
        # Delete existing entries in the table
        sql = "DELETE FROM %s" % (tabName)
        cursor.execute(sql)

        # Insert into the database
        sql = 'INSERT OR REPLACE INTO %s(%s) ' % \
              (tabName, ', '.join(self.tbl[tabName].dtype.names))
        sql += 'VALUES(%s) ' % (', '.join(['?']*len(self.tbl[tabName][0])))
        cursor.executemany(sql, irecarray_to_py(self.tbl[tabName]))
        
        # Clean up
        conn.commit()
        cursor.close()
        conn.close()
                    
