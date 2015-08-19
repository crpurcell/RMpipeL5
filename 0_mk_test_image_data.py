#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     mk_test_image_data.py                                             #
#                                                                             #
# USAGE:    ./mk_test_data.py                                                 #
#                                                                             #
# PURPOSE:  Create a small FITS dataset for the purposes of testing the       #
#           RM-pipeline. The script outputs a set of IQU image planes         #
#           containing unresolved sources. Edit the values at the top of the  #
#           script and run.                                                   #
#                                                                             #
# MODIFIED: 19-August-2015 by C. Purcell                                      #
#                                                                             #
#=============================================================================#

# Session path/name
sessionPath = "testSession/"

# Data cube parameters
startFreq_Hz =  1.0e9
endFreq_Hz = 3.0e9
nChans = 30
rmsNoise_mJy = 0.1
beamMinFWHM_deg = 1.5/3600.0
beamMajFWHM_deg = 2.5/3600.0
beamPA_deg = 20.0
pixScale_deg = 0.3/3600.0
xCent_deg = 90.0
yCent_deg = 0.0
nPixRA = 120
nPixDec = 100

# Coordinate system ["EQU" or "GAL"]
coordSys = "GAL"

# Properties of injected Faraday thin point sources
#-------------------------------------------------------------#
# [0]        [1] [2]      [3]       [4]       [5]     [6]
# fluxI_mJy, SI, fracPol, psi0_deg, RM_radm2, x_deg, y_deg
#-------------------------------------------------------------#
srcInLst = [ [5.0, -0.7, 0.2, 30.0, 19.0, 89.9981, +0.0001],
             [10.0, -0.1, 0.6, 80.0, 10.0, 90.0010, -0.0020],
             [15.0, -0.5, 0.7, 10.0, -50.0, 90.0016, +0.0023],
             [9.0,  0.0, 0.1, 0.0, 60.0, 89.9979, +0.0025],
             [7.0, -0.2, 0.5, 45.0, -32.0, 89.9970, -0.0027],
             [2.3, +0.5, 0.7, 120.0, -90.0, 90.0026, +0.000],
             [0.3, +0.0, 0.1, 120.0, -90.0, 90.0, +0.000]]
#-------------------------------------------------------------#
freq0_Hz = startFreq_Hz  # Frequency at which the flux is specified

# END USER EDITS -------------------------------------------------------------#

import os
import sys
import argparse
import shutil
import math as m
import numpy as np
import astropy.io.fits as pf
import astropy.wcs.wcs as pw

from Imports.util_PPC import twodgaussian
from Imports.util_PPC import create_IQU_spectra_RMthin
from Imports.util_FITS import strip_fits_dims
from Imports.util_FITS import create_simple_fits_hdu

C = 2.99792458e8

#-----------------------------------------------------------------------------#
def main():
    """
    Start the create_IQU_fits_data function if called from the command line.
    """
    
    # Help string to be shown using the -h option
    descStr = """
    Create a new dataset directory and populate it with FITS files containing
    Stokes I, Q and U images. Each image corresponds to a single frequency
    channel in a data-cube. There are equal numbers of files for each of the
    3 Stokes parameters whose names are formatted as e.g.,'CH23_StokesQ.fits'.
    The data is populated with a number of unresolved, Faraday thin sources
    defined in a table at the the top of this script:

      [ [fluxI_mJy, SI, fracPol, psi0_deg, RM_radm2, x_deg, y_deg], ... ]
    
      fluxI_mJy ... flux of source in first channel (mJy)
      SI        ... frequency spectral index
      fracPol   ... fractional polarisation (constant across frequency range)
      psi0_deg  ... intrinsic polarisation angle )
      RM_radm2  ... rotation measure (rad/m^2)
      x_deg     ... X-position coordinate (deg)
      y_deg     ... Y-position coordinate (deg)

    The beam may be elliptical and is set to a constant angular size as a
    function of frequency (i.e., assumes all data has been convolved to the
    same resolution). Please edit the variables at the top of the script to
    change the properties of the output data.

    In addition to the FITS files, the script outputs a simple ASCII catalogue
    and a SQL description of that catalogue. The catalogue file is used to
    drive the pipeline and the SQL descripton file tells the pipeline the
    format of the catalogue. This allows the user to define custom columns in
    the input catalogue, which are then incorporated into the results database.

    Example:

    ./0_mk_test_image_data.py testData/
    """

    # Parse the command line options
    parser = argparse.ArgumentParser(description=descStr,
                                 formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("dataPath", metavar="PATH/TO/DATA",
                        default="testData/", nargs="?",
                        help="Path to new data directory [testData/]")
    args = parser.parse_args()
    dataPath = args.dataPath

    # Call the function to create FITS data
    create_IQU_fits_data(dataPath, startFreq_Hz, endFreq_Hz, nChans,
                         rmsNoise_mJy, beamMinFWHM_deg, beamMajFWHM_deg,
                         beamPA_deg, pixScale_deg, xCent_deg, yCent_deg,
                         nPixRA, nPixDec, coordSys)

    # Print summary to user
    catFile = dataPath.rstrip("/") + "/testCat.dat"
    sqlFile = dataPath.rstrip("/") + "/testCatDesc.sql"
    print
    print "-" * 80
    print ">>> How to run the RM pipeline:"
    print "-" * 80
    print "A test dataset had been created in the directory '%s/':" \
          % dataPath.rstrip("/")
    print "> %d FITS images for each Stokes parameter IQU (one per channel)" \
          % nChans
    print "> A simple catalogue in the file '%s' " % catFile
    print "> A SQL catalogue description in the file '%s'" % sqlFile
    print
    print "To run the RM-pipeline execute the following commands in order:"
    print
    print "./1_verify_image_data.py %s/" % dataPath.rstrip("/")
    print "./2_create_image_session.py %s/ %s/ %s %s" % \
          (dataPath.rstrip("/"), sessionPath.rstrip("/"), catFile, sqlFile)
    print "# Edit the file '%s/inputs.config' (optional)" \
          % sessionPath.rstrip("/")
    print "./3_extract_spectra.py %s/"  % sessionPath.rstrip("/")
    print "./4_do_RM-synthesis.py %s/" % sessionPath.rstrip("/")
    print "./5_do_RM-clean.py %s/" % sessionPath.rstrip("/")
    print
    print "NOTE: information and help on each script can be viewed by ",
    print "executing each\ncommand followed by a '-h' flag, e.g.: \n"
    print "./0_mk_test_image_data.py -h"
    print 


#-----------------------------------------------------------------------------#
def create_IQU_fits_data(dataPath, startFreq_Hz, endFreq_Hz, nChans,
                         rmsNoise_mJy, beamMinFWHM_deg, beamMajFWHM_deg,
                         beamPA_deg, pixScale_deg, xCent_deg, yCent_deg,
                         nPixRA, nPixDec, coordSys="EQU"):

    """
    Create a set of FITS images corresponding to a Stokes I Q & U data-cube.
    """

    # Create the output directory path
    print "Creating test dataset in '%s'" % dataPath
    dirs = dataPath.rstrip("/").split("/")
    for i in range(1, len(dirs)):
        dirStr = "/".join(dirs[:i])
        if not os.path.exists(dirStr):
            os.mkdir(dirStr)
    if os.path.exists(dataPath):
        shutil.rmtree(dataPath, True)
    os.mkdir(dataPath)

    # Create a catalogue file
    catFile = dataPath + "/testCat.dat"
    FH = open(catFile, "w")
    FH.write("#Name  x_deg  y_deg\n")
    for i in range(len(srcInLst)):
        FH.write("Source%d %f %f\n" % ((i+1), srcInLst[i][5], srcInLst[i][6]))
    FH.close()

    # Create a catalogue description file
    sqlFile = dataPath + "/testCatDesc.sql"
    FH = open(sqlFile, "w")
    descStr = """
CREATE TABLE sourceCat (
uniqueName varchar(20),
x_deg double,
y_deg double);
    """
    FH.write("%s\n" % descStr)
    FH.close()

    # Create a simple HDU template
    freqArr_Hz = np.linspace(startFreq_Hz, endFreq_Hz, nChans)
    dFreqArr_Hz = np.diff(freqArr_Hz)
    hduTmp = create_simple_fits_hdu(shape=(1, 1, nPixDec, nPixRA),
                                    freq_Hz=freqArr_Hz[0],
                                    dFreq_Hz=dFreqArr_Hz[0],
                                    xCent_deg=xCent_deg,
                                    yCent_deg=yCent_deg,
                                    beamMinFWHM_deg=beamMinFWHM_deg,
                                    beamMajFWHM_deg=beamMajFWHM_deg,
                                    beamPA_deg=beamPA_deg,
                                    pixScale_deg=pixScale_deg,
                                    stokes="I",
                                    system=coordSys)
    head2D = strip_fits_dims(header=hduTmp.header, minDim=2, forceCheckDims=4)
    wcs2D = pw.WCS(head2D)
    shape2D = (nPixDec, nPixRA)

    # Calculate some beam parameters
    gfactor = 2.0*m.sqrt(2.0*m.log(2.0))
    beamMinSigma_deg = beamMinFWHM_deg/gfactor
    beamMajSigma_deg = beamMajFWHM_deg/gfactor
    beamMinSigma_pix = beamMinSigma_deg/pixScale_deg
    beamMajSigma_pix = beamMajSigma_deg/pixScale_deg
    beamPA_rad = m.radians(beamPA_deg)
    
    # Loop through the sources, calculate the spectra and pixel position
    spectraILst = []
    spectraQLst = []
    spectraULst = []
    coordLst_pix = []
    for row in srcInLst:
        IArr_mJy, QArr_mJy, UArr_mJy =\
                  create_IQU_spectra_RMthin(freqArr_Hz,
                                            row[0]/1e3,  # fluxI_mJy -> Jy
                                            row[1],      # SI
                                            row[2],      # fracPol
                                            row[3],      # psi0_deg
                                            row[4],      # RM_radm2
                                            freq0_Hz)
        spectraILst.append(IArr_mJy)
        spectraQLst.append(QArr_mJy)
        spectraULst.append(UArr_mJy)
        [ (x_pix, y_pix) ] = wcs2D.wcs_world2pix([ (row[5], row[6]) ], 0)
        coordLst_pix.append([x_pix, y_pix])
        
    # Loop through the frequency channels & create IQU files for each 
    for iChan in range(len(freqArr_Hz)):
        for iSrc in range(len(srcInLst)):
            params = [spectraILst[iSrc][iChan],  # amplitude
                      coordLst_pix[iSrc][0],     # X centre (pix)
                      coordLst_pix[iSrc][1],     # Y centre
                      beamMinSigma_pix,          # width (sigma)
                      beamMajSigma_pix,          # height (sigma)
                      beamPA_rad]                # PA (rad) W of N (clockwise)
            planeI = twodgaussian(params, shape2D).reshape(hduTmp.data.shape)
            params[0] = spectraQLst[iSrc][iChan]
            planeQ = twodgaussian(params, shape2D).reshape(hduTmp.data.shape)
            params[0] = spectraULst[iSrc][iChan]
            planeU = twodgaussian(params, shape2D).reshape(hduTmp.data.shape)
            if iSrc==0:
                dataIArr = planeI
                dataQArr = planeQ
                dataUArr = planeU
            else:
                dataIArr += planeI
                dataQArr += planeQ
                dataUArr += planeU

        # Add the noise
        dataIArr += np.random.normal(scale=rmsNoise_mJy/1e3,  # mJy -> Jy
                                     size=hduTmp.data.shape)
        dataQArr += np.random.normal(scale=rmsNoise_mJy/1e3,  # mJy -> Jy
                                     size=hduTmp.data.shape)
        dataUArr += np.random.normal(scale=rmsNoise_mJy/1e3,  # mJy -> Jy
                                     size=hduTmp.data.shape)

        # Write to the FITS files
        print "Writing FITS files for channel %s ..." % (iChan+1),
        sys.stdout.flush()
        hduTmp.header["CRVAL3"] = freqArr_Hz[iChan]
        hduTmp.header["CRVAL4"] = 1.0
        fitsFileOut =  dataPath + "/CH%d_StokesI.fits" % (iChan+1)
        if os.path.exists(fitsFileOut):
            os.remove(fitsFileOut)
        pf.writeto(fitsFileOut, dataIArr, hduTmp.header)
        hduTmp.header["CRVAL4"] = 2.0
        fitsFileOut = dataPath + "/CH%d_StokesQ.fits" % (iChan+1)
        if os.path.exists(fitsFileOut):
            os.remove(fitsFileOut)
        pf.writeto(fitsFileOut, dataQArr, hduTmp.header)
        hduTmp.header["CRVAL4"] = 3.0
        fitsFileOut =  dataPath + "/CH%d_StokesU.fits" % (iChan+1)
        if os.path.exists(fitsFileOut):
            os.remove(fitsFileOut)
        pf.writeto(fitsFileOut, dataUArr, hduTmp.header)
        print "done."
        sys.stdout.flush()


#-----------------------------------------------------------------------------#
if __name__=="__main__":
    main()
