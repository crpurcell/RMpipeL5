#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     verify_image_data.py                                              #
#                                                                             #
# USAGE:    ./verify_image_data.py PATH/TO/DATA-DIRECTORY [-I patI] [-Q patQ] #
#                                                       [-U patU] [-h --help] #
#                                                                             #
# PURPOSE:  Read and verify the Stokes I, Q & U image-FITS data, perform      #
#           simple sanity checks and write vectors of frequency, channel and  #
#           file to disk.                                                     #
#                                                                             #
# MODIFIED: 21-Apr-2015 by C. Purcell                                         #
#                                                                             #
#=============================================================================#

# Default patterns to match the Stokes I Q and U files.
patIdefault = '*.iconvol.fits'
patQdefault = '*.qconvol.fits'
patUdefault = '*.uconvol.fits'

#-----------------------------------------------------------------------------#

import os
import sys
import re
import glob
import argparse
import numpy as np
import astropy.io.fits as pf

from Imports.util_PPC import sort_nicely


#-----------------------------------------------------------------------------#
def main():
    """
    Start the verify_data function if called from the command line.
    """
    
    # Help string to be shown using the -h option
    descStr = """
    Scan a directory for FITS format files containing Stokes I, Q and U data.
    Each FITS file should contain a single image (frequency plane) from a cube
    of data.  The script sorts the list of files by frequency (read from the
    CRVAL3 header key) and writes vectors of frequency, channel and filename to
    ASCII text files in the data directory. Files for each Stokes parameter
    should have a unique naming format, matched by the wildcard patterns.
    Default patterns are set at the top of the script and may be overridden
    using command line arguments.
    """

    # Parse the command line options
    parser = argparse.ArgumentParser(description=descStr)
    parser.add_argument('dataPath', metavar='PATH/TO/DATA', default='.',
                        nargs='?', help='Path to data directory [.]')
    parser.add_argument('-I', dest='patI', default=patIdefault,
                     help='Pattern to match Stokes I files [%s]' % patIdefault)
    parser.add_argument('-Q', dest='patQ', default=patQdefault,
                     help='Pattern to match Stokes Q files [%s]' % patQdefault)
    parser.add_argument('-U', dest='patU', default=patUdefault,
                     help='Pattern to match Stokes U files [%s]' % patUdefault)
    args = parser.parse_args()
    dataPath = args.dataPath
    patI = args.patI
    patQ = args.patQ
    patU = args.patU

    # Call the verify function
    verify_image_data(dataPath, patI, patQ, patU)


#-----------------------------------------------------------------------------#
def verify_image_data(dataPath, patI, patQ, patU):
    """
    Scan for Stokes I, Q & U files matching the wildcard patterns in the data
    directory. Loop through the files, read their frequencies from the CDELT3
    header card and write an asscending vector of frequency to an ASCII text
    file. Also write lists of channel numbers and separate lists of I, Q & U
    files in frequency order. The script assumes that each file is a single
    plane drawn from a larger cube.
    """
    
    # Sanity checks
    if not os.path.exists(dataPath):
        print "\nErr: The directory '%s' does not exist." % dataPath
        sys.exit()
    dataPath = dataPath.rstrip('/')

    # Find all of the Stokes I,Q,U files in the directory
    print "\nPatterns = '%s', '%s', '%s'" % (patI, patQ, patU)
    print "Scanning the directory '%s' ..." % dataPath
    dataILst = glob.glob(dataPath + '/' + patI)
    sort_nicely(dataILst)
    nIfiles = len(dataILst)
    dataQLst = glob.glob(dataPath + '/' + patQ)
    sort_nicely(dataQLst)
    nQfiles = len(dataQLst)
    dataULst = glob.glob(dataPath + '/' + patU)
    sort_nicely(dataULst)
    nUfiles = len(dataULst)

    # Check that the same number of files was found for each pattern
    print "Found %d I files, %d Q files and %d U files.\n" % (nIfiles,
                                                              nQfiles,
                                                              nUfiles)
    if nIfiles!=nQfiles or nIfiles!=nUfiles:
        print "Err: script cannot process unequal numbers of files."
        print "If necessary use NaN-filled files to replace missing data."
        sys.exit()
    if nIfiles==0:
        print "Err: No matching files found."
        sys.exit()
        
    # Lists to store chan, freq
    freqLst = []

    # Loop through the files in lock-step
    print "Scanning through the files ..."
    for i in range(len(dataILst)):

        # Some feedback
        print os.path.split(dataILst[i])[-1],
        print os.path.split(dataQLst[i])[-1],
        print os.path.split(dataULst[i])[-1]
        
        # Read the headers
        headI = pf.getheader(dataILst[i])
        headQ = pf.getheader(dataQLst[i])
        headU = pf.getheader(dataULst[i])

        # Check the number of axes are correct and >= 3
        naxisI = headI['NAXIS']
        naxisQ = headI['NAXIS']
        naxisU = headI['NAXIS']
        if naxisI!=naxisQ or naxisI!=naxisU:
            print "Err: The number of dimensions in each file do not match."
            print "[%s, %s, %s]" % (naxisI, naxisQ, naxisU)
            sys.exit()
        if naxisI<3:
            print "Err: Less than 3 data axes found [NAXIS = %d]." % naxisI
            sys.exit()

        # Check the shape of the images are the same
        shapeI = (headI['NAXIS2'], headI['NAXIS2'])
        shapeQ = (headQ['NAXIS2'], headQ['NAXIS2'])
        shapeU = (headU['NAXIS2'], headU['NAXIS2'])
        if shapeI!=shapeQ or shapeI!=shapeU:
            print "Err: The shape of the three images arrays do not match."
            print "[%s, %s, %s]" % (shapeI, shapeQ, shapeU)
            sys.exit()
        
        # Check the frequencies are the same (assume CRVAL3=freq)
        freqI = headI['CRVAL3']
        freqQ = headQ['CRVAL3']
        freqU = headU['CRVAL3']
        if freqI!=freqQ or freqI!=freqU:
            print "The frequencies of the three Stokes files do not match."
            print "[%s, %s, %s]" % (freqI, freqQ, freqU)
            sys.exit()

        # Record the parameters
        freqLst.append(freqI)

    # Sort the filenames, channels into frequency order
    multiLst = zip(freqLst,
                   [os.path.split(x)[-1] for x in dataILst],
                   [os.path.split(x)[-1] for x in dataQLst],
                   [os.path.split(x)[-1] for x in dataULst])
    multiLst.sort()
    freqLstS, dataILstS, dataQLstS, dataULstS = zip(*multiLst)

    # Save an ordered frequency list
    freqFile = dataPath + '/freqs_Hz.txt'
    print "\nSaving ascending frequency vector to '%s'." % freqFile
    if os.path.exists(freqFile):
        os.remove(freqFile)
    np.savetxt(freqFile, freqLstS)

    # Save the list of I,Q & U files
    dataFile = dataPath + '/fileLstI.txt'
    print "Saving list of Stokes I files to '%s'." % dataFile
    if os.path.exists(dataFile):
        os.remove(dataFile)
    np.savetxt(dataFile, dataILstS, fmt='%s')
    dataFile = dataPath + '/fileLstQ.txt'
    print "Saving list of Stokes Q files to '%s'." % dataFile
    if os.path.exists(dataFile):
        os.remove(dataFile)
    np.savetxt(dataFile, dataQLstS, fmt='%s')
    dataFile = dataPath + '/fileLstU.txt'
    print "Saving list of Stokes U files to '%s'." % dataFile
    if os.path.exists(dataFile):
        os.remove(dataFile)
    np.savetxt(dataFile, dataULstS, fmt='%s')

    # Note the type of data in a file
    typeFile = dataPath + '/dataType.txt'
    FH = open(typeFile, 'w')
    FH.write("FITS_planes\n")
    FH.close()

#-----------------------------------------------------------------------------#
if __name__=="__main__":
    main()
