#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     util_DB.py                                                        #
#                                                                             #
# PURPOSE:  Functions to interface with a sqlite3 database.                   #
#                                                                             #
#                                                                             #
# MODIFIED: 28-Nov-2014 by C. Purcell                                         #
#                                                                             #
# CONTENTS:                                                                   #
#                                                                             #
#  register_sqlite3_numpy_dtypes ... setup sqlite3 to read numpy arrays       #
#  create_DB                ... create an empty DB using a dict of CREATE SQL #
#  insert_arr_db            ... insert recarray entries into the database     #
#  update_arr_db            ... update DB entries using a recarray            #
#  select_into_arr          ... run a SQL query and return a numpy recarray   #
#  schema_to_tabledef       ... parse the SQL table definitions               #
#                                                                             #
#=============================================================================#

import os
import sys
import shutil
import re
import numpy as np
import sqlite3

from util_rec import *


#-----------------------------------------------------------------------------#
def register_sqlite3_numpy_dtypes():
    """
    Map numpy data-types to the limited sqlite data-types. This must be called
    before using the sqlite3 database or INSERT statements will fail.
    """
    
    for t in (np.int8, np.int16, np.int32, np.int64,
          np.uint8, np.uint16, np.uint32, np.uint64):
        sqlite3.register_adapter(t, long)
    for t in (np.float16, np.float32, np.float64,
          np.float128, np.double):
        sqlite3.register_adapter(t, float)

    
#-----------------------------------------------------------------------------#
def create_db(dbFile, createSQLdict, LF=None):
    """
    Create an empty SQLite3 database on disk.
    """
    
    if os.path.exists(dbFile):
        os.remove(dbFile)
    try:
        conn = sqlite3.connect(dbFile)
        cursor = conn.cursor()
        for sql in createSQLdict.values():
            cursor.execute(sql)    
        conn.commit()
        cursor.close()
        conn.close()
    except Exception:
        return False
    
    return True

#-----------------------------------------------------------------------------#
def insert_arr_db(cursor, recArr, tabName, fieldNameLst=None,
                     insertSQL="INSERT OR REPLACE"):
    """
    Insert a numpy recarray into a database via a cursor object. It is assumed
    that the fields in the recarray and database have the same names and
    compatable datatypes. If the recarray contains fields NOT in the database
    the user must provide a list of field names to be to be inserted. These
    must be a subset of the field names in the table.
    """
    
    # Default to all fields
    if not fieldNameLst:
        fieldNameLst = recArr.dtype.names    
    
    sql = '%s INTO %s (%s) ' % (insertSQL, tabName, ', '.join(fieldNameLst))
    sql += 'VALUES(%s) ' % (', '.join(['?']*len(fieldNameLst)))
    cursor.executemany(sql, fields_view(recArr, fieldNameLst))


#-----------------------------------------------------------------------------#
def update_arr_db(cursor, recArr, tabName, keyName, fieldNameLst=None):
    """
    Do an UPDATE on existing rows
    """
    
    # Default to all fields
    if not fieldNameLst:
        fieldNameLst = list(recArr.dtype.names)

    # key must exist in list of field names
    if not keyName in fieldNameLst:
        print "ERR: Key '%s' not in column list" % keyName
        return

    # Remove the key from the list and format the SQL
    fieldNameLst.remove(keyName)        
    sql = 'UPDATE %s SET ' % tabName
    sql += '=?, '.join(fieldNameLst) + '=? '
    sql += 'WHERE %s = ?' % keyName

    # Attach the key to the end of the field list and fetch a view
    # Use fancy indexing to force the key to the last column
    fieldNameLst.append(keyName)
    a = fields_view(recArr, fieldNameLst)[fieldNameLst]
    cursor.executemany(sql, a)
    
    
#-----------------------------------------------------------------------------#
def select_into_arr(cursor, sql, args=[]):
    """
    Run a SQL query and return a numpy recordarray.
    """

    if args == []:
        cursor.execute(sql)
    else:
        cursor.execute(sql, tuple(args))
    try:
        rows = cursor.fetchall()
        if len(rows) == 0:
            rows = np.array([], dtype='i4')
        else:
            columnNameLst = zip(*cursor.description)[0]
            rows = np.rec.fromrecords(rows, names=columnNameLst)
        return rows
    except Exception:
        print "WARNING: failed to convert SQL result to a recarray!"
        return None


#-----------------------------------------------------------------------------#
def schema_to_tabledef(schemaFile, addColDict={}):
    """
    Parse the SQL file containing the CREATE TABLE definitions and convert to
    a format that python can understand. This is then used to initialise the
    numpy record arrays which store the catalogue in memory. Assumes the
    create statement has the form:
    
    'CREATE TABLE myTableName (entry1 double [args+], entry2 int(2), ... );'

    The 'CREATE TABLE' statement must have the same case. Currently recognises
    only 'float', 'double', 'int(n)' and 'varchar(n)' types in the SQL.
    """

    # Return these
    tableNameLst = []
    tableDefDict = {}
    tableSQLdict = {}
    
    # Compile a few useful regular expressions
    spaces = re.compile('\s+')
    commaAndSpaces = re.compile(',\s+')
    comment = re.compile('#.*')
    quotes = re.compile('\'[^\']*\'')
    brackets = re.compile('[\[|\]\(|\)|\{|\}]')
    createRe =  re.compile('^(CREATE TABLE|create table) (\w+)\s*\((.+)\)\s*$')
    
    # Translation function for data types: SQL to recarray
    # RecArray dtypes:
    #   bytes                = b<n>, e.g. b1
    #   ints                 = i<n>, e.g. i2, i4, i8,
    #   unsigned ints        = u<n>, e.g. u1, u2, u4, u8
    #   floats               = f<n>, e.g. f2, f4, f8
    #   complex              = c<n>, e.g. c8, c16
    #   fixed length strings = a<n>, e.g. a10, a100
    # where <n> is the number of bytes / chars, so float32=f4, float64=f8
    def trtype(dtype):
        floatRe = re.compile('^float')
        doubleRe = re.compile('^double')
        intRe = re.compile('^int')
        charRe = re.compile('^varchar\((\d+)\)')
        if floatRe.match(dtype):
            return 'f4'
        if doubleRe.match(dtype):
            return 'f8'
        if intRe.match(dtype):
            return 'i8'
        mch = charRe.match(dtype)
        if mch:
             return 'a' + mch.group(1)
        return 'f8'    # default to float64
        
    # Loop through the SQL statements
    fileStr=''
    FH = open(schemaFile)
    for line in FH:
        line = line.replace('\r', ' ')        # kill carriage-return
        line = line.replace('\n', ' ')        # kill newlines
        if not comment.match(line):
            fileStr += line
    sqlLst = fileStr.split(';')
    FH.close()
    for sql in sqlLst:

        # Simplify the SQL statement
        sql = sql.replace('\r', ' ')        # kill carriage-return
        sql = sql.replace('\n', ' ')        # kill newlines
        sql = sql.strip()                   # kill external whitespace
        sql = spaces.sub(' ', sql)          # shrink internal whitespaces
        sql = commaAndSpaces.sub(',', sql)  # kill ambiguous spaces

        mch = createRe.match(sql)
        if mch:            
            tableName = mch.group(2).strip()
            colDefStr = mch.group(3).strip()
            
            # Add in columns if required
            if tableName in addColDict:
                colDefStr += ",%s" % addColDict[tableName]
            tableSQLdict[tableName] = "CREATE TABLE %s (%s)" % \
                                      (tableName, colDefStr)
            tableNameLst.append(tableName)
            colDefLst = colDefStr.strip().split(',')
            colDefLst = [x.split(' ')[:2] for x in colDefLst]

            # Translate the data types into a python recarray dtype list
            for i in range(len(colDefLst)):
                colDefLst[i][1] = trtype(colDefLst[i][1])
                colDefLst[i] = tuple(colDefLst[i])

            # Add to the table definition dictionary
            tableDefDict[tableName] = colDefLst

    return tableDefDict, tableSQLdict
