#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     util_rec.py                                                       #
#                                                                             #
# PURPOSE:  Functions for operating on python record arrays.                  #
#                                                                             #
# MODIFIED: 24-September-2015 by C. Purcell                                   #
#                                                                             #
# CONTENTS:                                                                   #
#                                                                             #
#  pyify                    ... return type converters given type strings     #
#  irecarray_to_py          ... convert a recarray into a list                #
#  fields-view              ... return a view of chosen fields in a recarray  #
#                                                                             #
#=============================================================================#

import numpy as np


#-----------------------------------------------------------------------------#
def pyify(typestr):
    """
    Return a Python :class:'type' that most closely represents the
    type encoded by *typestr*
    """
    if typestr[1] in 'iu':
        return int
    elif typestr[1] == 'f':
        return float
    elif typestr[1] == 'S':
        return str
    return lambda x: x


#-----------------------------------------------------------------------------#
def irecarray_to_py(a):
    """
    Slow conversion of a recarray into a list of records with python types.
    Get the field names from :attr:'a.dtype.names'.
    :Returns: iterator so that one can handle big input arrays
    """
    pytypes = [pyify(typestr) for name,typestr in a.dtype.descr]
    def convert_record(r):
        return tuple([converter(value) for converter,
                      value in zip(pytypes,r)])
    return (convert_record(r) for r in a)


#-----------------------------------------------------------------------------#
def fields_view(arr, fieldNameLst=None):
    """
    Return a view of a numpy record array containing only the fields names in
    the fields argument. 'fields' should be a list of column names.
    """
    
    # Default to all fields
    if not fieldNameLst:
        fieldNameLst = arr.dtype.names
    dtype2 = np.dtype({name:arr.dtype.fields[name] for name in fieldNameLst})
    
    return np.ndarray(arr.shape, dtype2, arr, 0, arr.strides)
