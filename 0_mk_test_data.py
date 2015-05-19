#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     mk_test_data.py                                                   #
#                                                                             #
# USAGE:    ./mk_test_data.py                                                 #
#                                                                             #
# PURPOSE:  Create a small FITS dataset for the purposes of testing the       #
#           RM-pipeline. Edit the values at the top of the script and run.    #
#                                                                             #
# MODIFIED: 19-May-2015 by C. Purcell                                         #
#                                                                             #
#=============================================================================#

# Dataset path/name
dataPath = "testData"

# Proposed session name
sessionPath = "testSession"

# Data cube parameters
startFreq_Hz =  1.0e9
endFreq_Hz = 3.0e9
nChans = 300
rmsNoise_mJy = 0.01
beamFWHM_deg = 1.5/3600.0
pixScale_deg = 0.3/3600.0
xCent_deg = 90.0
yCent_deg = 0.0
nPixRA = 100
nPixDec = 100

# Coordinate system ['EQU' or 'GAL']
coordSys = "EQU"

# Properties of injected Faraday thin point sources
#-------------------------------------------------------------#
# [0]        [1] [2]      [3]       [4]       [5]     [6]
# fluxI_mJy, SI, fracPol, psi0_deg, RM_radm2, x_deg, y_deg
#-------------------------------------------------------------#
#srcInLst = [ [5.0, -0.0, 0.2, 30.0, 19.0, 89.9981, +0.0001],
#             [10.0, -0.1, 0.6, 80.0, 10.0, 90.0010, -0.0020],
#             [15.0, -0.5, 0.7, 10.0, -50.0, 90.0016, +0.0023] ]
srcInLst = [ [5.0, -0.7, 0.2, 30.0, 19.0, 89.9981, +0.0001] ]

#-----------------------------------------------------------------------------#

import os
import sys
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
srcName varchar(20),
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
                                    beamFWHM_deg=beamFWHM_deg,
                                    pixScale_deg=pixScale_deg,
                                    stokes='I',
                                    system=coordSys)
    head2D = strip_fits_dims(header=hduTmp.header, minDim=2, forceCheckDims=4)
    wcs2D = pw.WCS(head2D)
    shape2D = (nPixDec, nPixRA)

    # Calculate some beam parameters
    gfactor = 2.0*m.sqrt(2.0*m.log(2.0))
    beamSigma_deg = beamFWHM_deg/gfactor
    beamSigma_pix = beamSigma_deg/pixScale_deg
    
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
                                            row[4])      # RM_radm2
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
                      beamSigma_pix,             # width (sigma)
                      beamSigma_pix,             # height (sigma)
                      0.0]                       # PA (rad) W of N (clockwise)
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
        hduTmp.header['CRVAL3'] = freqArr_Hz[iChan]
        hduTmp.header['CRVAL4'] = 1.0
        fitsFileOut =  dataPath + "/CH%d_StokesI.fits" % (iChan+1)
        if os.path.exists(fitsFileOut):
            os.remove(fitsFileOut)
        pf.writeto(fitsFileOut, dataIArr, hduTmp.header)
        hduTmp.header['CRVAL4'] = 2.0
        fitsFileOut = dataPath + "/CH%d_StokesQ.fits" % (iChan+1)
        if os.path.exists(fitsFileOut):
            os.remove(fitsFileOut)
        pf.writeto(fitsFileOut, dataQArr, hduTmp.header)
        hduTmp.header['CRVAL4'] = 3.0
        fitsFileOut =  dataPath + "/CH%d_StokesU.fits" % (iChan+1)
        if os.path.exists(fitsFileOut):
            os.remove(fitsFileOut)
        pf.writeto(fitsFileOut, dataUArr, hduTmp.header)
        print "done."
        sys.stdout.flush()

    # Print summary to user
    print
    print "-" * 80
    print ">>> Running the RM pipeline"
    print "-" * 80
    print "A test dataset had been created in the directory '%s/':" % dataPath
    print "> %d FITS images for each Stokes parameter IQU (one per channel)" \
          % len(freqArr_Hz)
    print "> A simple catalogue in the file '%s' " % catFile
    print "> A SQL catalogue description in the file '%s'" % sqlFile
    print
    print "To run the RM-pipeline execute the following commands in order:"
    print
    print "./1_verify_image_data.py %s/" % dataPath
    print "./2_create_image_session.py %s/ %s/ %s %s" % \
          (dataPath, sessionPath, catFile, sqlFile)
    print "# Edit the file '%s/inputs.config' (optional)" \
          % sessionPath
    print "./3_extract_spectra.py %s/"  % sessionPath
    print "./4_do_RM-synthesis.py %s/" % sessionPath
    print "./5_do_RM-clean.py %s/" % sessionPath


#-----------------------------------------------------------------------------#
if __name__=="__main__":
    main()
