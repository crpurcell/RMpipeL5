#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     verify_ascii_data.py                                              #
#                                                                             #
# USAGE:    ./verify_ascii_data.py PATH/TO/DATA-DIRECTORY [-m PATTERN]        #
#                                       [-f FILENAME] [-c COLUMN] [-h --help] #
#                                                                             #
# PURPOSE:  Read and verify Stokes I, Q & U spectra in ASCII '.dat' files.    #
#                                                                             #
# MODIFIED: 20-November-2015 by C. Purcell                                    #
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

# Default wildcard patterns to match the ASCII data files.
patDefault = '*.dat'

# END USER EDITS -------------------------------------------------------------#

import os
import sys
import re
import glob
import argparse
import numpy as np

from Imports.util_PPC import sort_nicely


#-----------------------------------------------------------------------------#
def main():
    """
    Start the verify_data function if called from the command line.
    """
    
    # Help string to be shown using the -h option
    descStr = """
    Scan a directory for ASCII files containing Stokes I, Q and U spectra.
    Each file should contain Stokes spectra in four columns:
    
       [frequency_Hz, StokesI_Jy, StokesQ_Jy, StokesU_Jy]
       
    The script reads in the spectra, sorts the channels by assending frequency
    and re-writes the vectors to the ASCII files. The script can use wildcard
    pattern matching to filter the files in the directory, or read a list of
    filenames from a text file. 

    Examples:
    
    ./1_verify_ascii_data.py -f *I.fits -Q *Q.fits -U *U.fits testData/
    """

    # Parse the command line options
    parser = argparse.ArgumentParser(description=descStr,
                                 formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('dataPath', metavar='PATH/TO/DATA', default='.',
                        nargs='?', help='Path to data directory [.]')
    parser.add_argument('-m', dest='pattern', default=patDefault,
                      help='Pattern to match the data files [%s]' % patDefault)
    parser.add_argument('-f', dest='listFile',
                        help='File containing list of data files')
    parser.add_argument('-c', dest='col', type=int, nargs='?', default=4,
                    help='Column in listFile containing the list of files [4]')
    args = parser.parse_args()
    dataPath = args.dataPath
    pattern = args.pattern
    listFile = args.listFile
    col = args.col
    
    # Call the verify function
    verify_image_data(dataPath, pattern, listFile, col)


#-----------------------------------------------------------------------------#
def verify_image_data(dataPath, pat, listFile=None, col=0):
    """
    Scan a data directory for ASCII files containing Stokes I, Q and U spectra.
    Loop through the files, read in the frequnecy vector and spectra. Sort the
    spectral channels by assending frequency and write back to the ASCII files.
    The script can use wildcard pattern matching to filter the files in the
    data directory, or read a list of filenames from a text file.
    """
    
    # Sanity checks
    if not os.path.exists(dataPath):
        print "\nErr: The directory '%s' does not exist." % dataPath
        sys.exit()
    dataPath = dataPath.rstrip('/')

    # Find all of the data files in the directory
    print "\nPattern = '%s'" % pat
    print "Scanning the directory '%s' ..." % dataPath
    dataLst = glob.glob(dataPath + '/' + pat)
    sort_nicely(dataLst)
    nFiles = len(dataLst)
    print "Found %d data files\n" % nFiles
    if nFiles==0:
        print "Err: No matching files found."
        sys.exit()
    
    # Open a catalogue file
    catFile = dataPath + "/catalogue.txt"
    catFH = open(catFile, "w")
    catFH.write("#Name fileName x_deg  y_deg\n")

    # Loop through the input data files
    print "Scanning through the files ..."
    for i in range(len(dataLst)):

        # Some feedback
        print os.path.split(dataLst[i])[-1]

        # Read in the spectra
        multiArr = np.loadtxt(dataLst[i], unpack=True)
        nCols = multiArr.shape[0]
        if nCols!=4:
            print "Err: Expecting 4 columns in ASCII file, found %d." % nCols
            continue
        
        # Sort the spectra into ascending frequency order
        idx = np.argsort(multiArr[0])
        if not np.all(idx==np.indices(multiArr[0].shape)[0]):
            multiArr = multiArr[:, idx]
            print "Warn: Overwriting '%s' with sorted version." % dataLst[i]
            np.savetxt(dataLst[i], multiArr.transpose())
        freqArr_Hz, IArr_Jy, QArr_Jy, UArr_Jy = multiArr

        # Add to the catalogue file
        catFH.write("Source%d  %s %f %f\n" % ((i+1), dataLst[i], 0.0, 0.0))
        
    catFH.close()

    # Save an ordered frequency list
    freqFile = dataPath + "/freqs_Hz.txt"
    print "\nSaving ascending frequency vector to '%s'." % freqFile
    if os.path.exists(freqFile):
        os.remove(freqFile)
    np.savetxt(freqFile, freqArr_Hz)

    # Note the type of data in a file
    typeFile = dataPath + "/dataType.txt"
    print "Noting dataType='ASCII_spectra' in file '%s'." % typeFile
    FH = open(typeFile, "w")
    FH.write("ASCII_spectra\n")
    FH.close()
    

#-----------------------------------------------------------------------------#
if __name__=="__main__":
    main()
