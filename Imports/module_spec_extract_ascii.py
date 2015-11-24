#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     module_spec_extract_area.py                                       #
#                                                                             #
# PURPOSE:  Extract spectra from text files and copy into standard POSSUM     #
#           FITS files to be compatable with the FITS-porcessing steps.       #
#                                                                             #
# REQUIRED: Requires the numpy, scipy and astropy modules.                    #
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

import os
import sys
import copy
import math as m
import numpy as np
import traceback
import astropy.io.fits as pf
import astropy.wcs.wcs as pw

from mpfit import mpfit

from util_PPC import log_wr
from util_PPC import log_fail
from util_PPC import calc_sumbox_norm
from util_PPC import nanmedian
from util_PPC import MAD
from util_PPC import poly5

from util_FITS import mkWCSDict
from util_FITS import get_beam_from_header
from util_FITS import get_beam_area
from util_FITS import strip_fits_dims


#-----------------------------------------------------------------------------#
def mod_spec_extract_ascii(catRec, dataPath, polyOrd=2, outDataDir='./OUT',
                     doOverwrite=False, LF=None):
    """
    Run the procedure to extract data from ASCII files supplied by the user
    and copy to FITS files in the current session directory. For each catalogue
    entry and Stokes component a FITS file is created containing 1x1 pixel
    'cube', a dummy aperture mask and a spectra table with columns 'freq',
    'rms', 'src' and 'model'. A polynomial of order <=5 is fit to the Stokes I
    spectrum and stored in the 'model' column. Measurements made on each
    spectrum are written to a record array and returned.
    """

    # Default logging to STDOUT
    if LF is None:
        LF = sys.stdout

    # Create the output data directory
    if not os.path.exists(outDataDir):
        os.mkdir(outDataDir)
        
    # Check the catRec record array has the required fields
    reqCols = ['uniqueName','fileName']
    if not min([x in catRec.dtype.names for x in reqCols]):
        errStr = "Catalogue passed to the extraction module is " + \
                 "missing a required field: %s" % reqCols
        log_fail(LF, errStr)

    # Perform the extraction on the Stokes I data
    cat = (catRec['uniqueName'], catRec['fileName'])
    dataIFileLst, dataQFileLst, dataUFileLst  = \
                  extract_data_from_ascii(cat,
                                          dataPath,
                                          outDataDir,
                                          doOverwrite=doOverwrite,
                                          LF=LF)
    
    # Create a recArray to store the measurements
    dType = [('uniqueName', 'a20'),
             ('boxScale_pix', 'i8'),
             ('fluxMedI_Jybm', 'f8'),
             ('fluxMedQ_Jybm', 'f8'),
             ('fluxMedU_Jybm', 'f8'),
             ('rmsMedI_Jybm', 'f8'),
             ('rmsMedQ_Jybm', 'f8'),
             ('rmsMedU_Jybm', 'f8'),
             ('rmsMedQUAvg_Jybm', 'f8'),
             ('nFreqChan', 'i8'),
             ('deltaFreqChan_Hz', 'f8'),
             ('coeffPolyIspec', 'a500'),
             ('bmaj_deg', 'f8'),
             ('bmin_deg', 'f8'),
             ('bpa_deg', 'f8'),
             ('pixscale_deg', 'f8'),
             ('fNormSumbox', 'f8'),
             ('fitIredChiSq', 'f8'),
             ('fitIstatus', 'i8'),
             ('extractStatus', 'i8')]
    specRec = np.zeros(len(catRec), dtype=dType)

    # Loop through the sources and measure the spectra
    msgStr = "\n" + "-" * 80
    msgStr += "\nALL SPECTRA EXTRACTED, PROCEEDING TO MEASURE PROPERTIES"
    msgStr += "\n" + "-" * 80
    log_wr(LF, msgStr)
    for i in range(len(dataIFileLst)):
        log_wr(LF, " -> Measuring '%s'..." % catRec[i]['uniqueName'])
        try:
            HDUILst = pf.open(dataIFileLst[i], "update", memmap=True) 
            freqArr_Hz = HDUILst[2].data["freq"]
            specI  =  HDUILst[2].data["src"]
            rmsSpecI =   HDUILst[2].data["rms"]
            HDUQLst = pf.open(dataQFileLst[i], "update", memmap=True) 
            specQ  =  HDUQLst[2].data["src"]
            rmsSpecQ =   HDUQLst[2].data["rms"]
            HDUULst = pf.open(dataUFileLst[i], "update", memmap=True) 
            specU  =  HDUULst[2].data["src"]
            rmsSpecU =   HDUULst[2].data["rms"]
        except Exception:
            log_wr(LF, ">>> Warn: failed to load extracted spectra.")
            specRec[i]['extractStatus'] = 0
            continue

        # Calculate the average RMS and PI spectrum
        rmsSpecQUAvg = (rmsSpecQ + rmsSpecU) / 2.0
        
        # Fit a polynomial model to the Stokes I spectrum (order<=5)
        # Frequency axis must be in GHz to avoid overflow errors
        # Arrays must be float64
        try:
            p, fitIstatus, redChiSq = fit_spec_poly5(freqArr_Hz/1e9, specI,
                                                     rmsSpecI, polyOrd)
            specImodel = poly5(p)(freqArr_Hz/1e9)
        except Exception:
            fitIstatus= 0
            redChiSq=-1
            
        if fitIstatus<1:
            log_wr(LF, "> Warn: I fit failed on '%s'." % 
                   catRec[i]['uniqueName'] )
            specImodel = np.ones_like(specI)
            p = [0, 0, 0, 0, 0, 1]
            
        # Save the model I spectrum & coefficients
        HDUILst[2].data["model"] = specImodel
        
        # Add spectrum measurements 
        specRec[i]['uniqueName'] = catRec[i]['uniqueName']
        specRec[i]['fluxMedI_Jybm'] = nanmedian(specI)
        specRec[i]['fluxMedQ_Jybm'] = nanmedian(specQ)
        specRec[i]['fluxMedU_Jybm'] = nanmedian(specU)
        specRec[i]['rmsMedI_Jybm'] = nanmedian(rmsSpecI)
        specRec[i]['rmsMedQ_Jybm'] = nanmedian(rmsSpecQ)
        specRec[i]['rmsMedU_Jybm'] = nanmedian(rmsSpecU)
        specRec[i]['rmsMedQUAvg_Jybm'] = nanmedian(rmsSpecQUAvg)
        
        # Add other metadata to the catalogue
        specRec[i]['boxScale_pix'] = 1
        specRec[i]['nFreqChan'] = len(freqArr_Hz)
        specRec[i]['deltaFreqChan_Hz'] = np.nanmin(np.diff(freqArr_Hz))
        specRec[i]['bmaj_deg'] = 1.0
        specRec[i]['bmin_deg'] = 1.0
        specRec[i]['bpa_deg'] = 0.0
        specRec[i]['pixscale_deg'] = 1.0
        specRec[i]['fNormSumbox'] = 1.0
        specRec[i]['fitIredChiSq'] = redChiSq
        specRec[i]['fitIstatus'] = fitIstatus
        specRec[i]['coeffPolyIspec'] = ','.join([str(x) for x in p])
        specRec[i]['extractStatus'] = 1

        # Clean up
        HDUILst.close()
        HDUQLst.close()
        HDUULst.close()
        
    return specRec



#-----------------------------------------------------------------------------#
def extract_data_from_ascii(cat, dataPath, outDataDir='./', doOverwrite=False,
                            LF=None):
    """
    Extract data from a list of ASCII files and save to FITS files on disk.

    All data is saved as a single FITS file in the following extensions:
    EXT 1: Dummy cube covering with 1 pixel to mimic the POSSUM FITS format.
    EXT 2: Dummy extraction mask image.
    EXT 3: A FITS table containing the columns [frequency_Hz, src, rms].
    
    """

    # Unpack the catalogue
    nameArr, fileArr = cat

    # Create a list of root path names for output files
    saveRootLst = map(''.join, zip([outDataDir + '/']*len(nameArr), nameArr))

    # Loop through the list of ASCII files
    outFileLOL = [[],[],[]]
    for i in range(len(nameArr)):
        log_wr(LF, "> Reading ASCII file '%s' ..." % nameArr[i])
        multiArr = np.loadtxt(dataPath + "/" + fileArr[i], unpack=True)
        nCols = multiArr.shape[0]
        if nCols==4:
            freqArr_Hz, IArr_Jy, QArr_Jy, UArr_Jy = multiArr
            dIArr_Jy = np.zeros_like(IArr_Jy)
            dQArr_Jy = np.zeros_like(QArr_Jy)
            dUArr_Jy = np.zeros_like(UArr_Jy)
        elif nCols==7:
            [freqArr_Hz, IArr_Jy, QArr_Jy, UArr_Jy,
             dIArr_Jy, dQArr_Jy, dUArr_Jy] = multiArr
        else:
            log_write(LF, "Err: expecting 4 or 7 columns in ASCII file.")
            continue

        # Loop through the Stokes parameters
        specLst = [IArr_Jy, QArr_Jy, UArr_Jy]
        dSpecLst = [dIArr_Jy, dQArr_Jy, dUArr_Jy]
        suffixLst = ["I", "Q", "U"]
        for j in range(3):
            specArr = specLst[j]
            dSpecArr = dSpecLst[j]
            suffix = suffixLst[j]

            # Create a 1-pixel cube to store the spectrum
            cubeSpec = specArr.reshape((specArr.shape[0],1,1))
            hdu0 = pf.PrimaryHDU(cubeSpec)

            # Create a mask image to show the source pixels extracted
            maskData = np.ones_like(cubeSpec)
            hdu1 = pf.ImageHDU(maskData)
        
            # Create a record array table to store the spectrum
            col1 = pf.Column(name="freq", format="f4",
                             array=freqArr_Hz)
            col2 = pf.Column(name="rms", format="f4", 
                             array=dSpecArr)
            col3 = pf.Column(name="src", format="f4", 
                             array=specArr)
            col4 = pf.Column(name="model", format="f4", 
                             array=np.zeros_like(freqArr_Hz))
            hdu2 = pf.new_table([col1, col2, col3, col4])
                
            # Write the file to disk and clean up
            hduLst = pf.HDUList([hdu0, hdu1, hdu2])
            outFits = saveRootLst[i] + "_spec" + suffix  + ".fits"
            hduLst.writeto(outFits, output_verify="fix", clobber=True)
            hduLst.close()
            outFileLOL[j].append(outFits)

    return outFileLOL


#-----------------------------------------------------------------------------#
def get_ell_indices(nAxisX, nAxisY, xArr_pix, yArr_pix, minAx_pix, majAx_pix,
                    pa_rad):
    """
    Get the indices of pixels within an ellipse drawn on an image of shape
    (nAxisX, nAxisY), e.g., the footprint of a Gaussian fit. Works in pixel
    coordinates assuming no significant field curvature over the source area.

    TODO: Operate on a small square box instead of the whole image each time.
    """

    mskIndxLst = []

    # Create a pixel index array 
    indArr = np.indices((nAxisY, nAxisX), dtype="f4")

    # Loop through the input positions creating a mask
    for i in range(len(xArr_pix)):

        # Shift the zero position to the source coordinate
        dX = indArr[-1] -xArr_pix[i]
        dY = indArr[-2] -yArr_pix[i]
        paMath_rad = m.pi/2.0 - pa_rad[i]      # Astro to math PA convention

        # Construct an ellipse mask
        ellArr = ( (dX*np.cos(paMath_rad) + dY*np.sin(paMath_rad))**2.0
                   / majAx_pix[i]**2.0 + 
                   (dX*np.sin(paMath_rad) - dY*np.cos(paMath_rad))**2.0
                   / minAx_pix[i]**2.0 )

        # Save the indices of pixels within the ellipse
        mskIndxLst.append(np.where(ellArr<=1.0))
        
    return mskIndxLst


#-----------------------------------------------------------------------------#
def calc_cutout_bounds(nAxisX, nAxisY, xArr_pix, yArr_pix, cutoutRadius_pix):
    """
    Given an image of shape (nAxisX,xaxisY), determine the absolute pixel
    bounds of a square cutout box for each (xArr_pix[i], yArr_pixY[i]) entry
    in a catalogue. If the cutout box overlaps an edge it is shifted so
    that it butts against the edge instead. This function assumes that the
    cutout box is smaller than the image array.
    """
    
    xRndArr_pix = np.round(xArr_pix).astype('i4')
    yRndArr_pix = np.round(yArr_pix).astype('i4')

    lowLim = np.zeros_like(xRndArr_pix)
    highLim = np.zeros_like(xRndArr_pix) + nAxisX -1
    xMinArr_pix = xRndArr_pix - cutoutRadius_pix
    xMinArr_pix = np.max(zip(xMinArr_pix, lowLim), 1)
    xMaxArr_pix = xMinArr_pix + cutoutRadius_pix*2
    xMaxArr_pix = np.where(xMaxArr_pix>highLim,
                           highLim,
                           xMaxArr_pix)
    xMinArr_pix = np.where(xMaxArr_pix>=highLim,
                           xMaxArr_pix - cutoutRadius_pix*2,
                           xMinArr_pix)
    edgeFlag = np.zeros_like(xRndArr_pix)    
    edgeFlag = np.where(xMinArr_pix<=lowLim, 1, edgeFlag)
    edgeFlag= np.where(xMaxArr_pix>=highLim, 1, edgeFlag)

    highLim = np.zeros_like(yRndArr_pix) + nAxisY -1
    yMinArr_pix = yRndArr_pix - cutoutRadius_pix
    yMinArr_pix = np.max(zip(yMinArr_pix, lowLim), 1)
    yMaxArr_pix = yMinArr_pix + cutoutRadius_pix*2
    yMaxArr_pix = np.where(yMaxArr_pix>highLim,
                             highLim,
                             yMaxArr_pix)
    yMinArr_pix = np.where(yMaxArr_pix>=highLim,
                             yMaxArr_pix - cutoutRadius_pix*2,
                             yMinArr_pix)
    edgeFlag = np.where(yMinArr_pix<=lowLim, 1, edgeFlag)
    edgeFlag= np.where(yMaxArr_pix>=highLim, 1, edgeFlag)
    
    return [xMinArr_pix, xMaxArr_pix, yMinArr_pix, yMaxArr_pix, edgeFlag]


#-----------------------------------------------------------------------------#
def fit_spec_poly5(xData, yData, dyData = None, order=5):
    """
    Fit a <=5th order polynomial to a spectrum. To avoid overflow errors the
    X-axis data should not be large numbers (e.g.: x10^9 Hz; use GHz instead).
    """
    
    xData = np.array(xData, dtype='f8')
    yData = np.array(yData, dtype='f8')
    if dyData is None:
        dyData = np.ones_like(yData)
    else:
        dyData = np.array(dyData, dtype='f8')
    if np.all(dyData==0.0):
        dyData = np.ones_like(yData)
        
    # Limit the order between 1 and 5
    if order<1:
        order = 1
    if order>5:
        order = 5
        
    # Estimate initial coefficients
    C1 = nanmedian(np.diff(yData)) / nanmedian(np.diff(xData))
    ind = int(np.median(np.where(~np.isnan(yData))))
    C0 = yData[ind] - (C1 * xData[ind])
    C5 = 0.0
    C4 = 0.0
    C3 = 0.0
    C2 = 0.0
    inParms=[ {'value': C5, 'parname': 'C5'},
              {'value': C4, 'parname': 'C4'},
              {'value': C3, 'parname': 'C3'},
              {'value': C2, 'parname': 'C2'},
              {'value': C1, 'parname': 'C1'},
              {'value': C0, 'parname': 'C0'} ]

    # Set the polynomial order
    for i in range(len(inParms)):
        if len(inParms)-i-1>order:
            inParms[i]['fixed'] = True
        else:
            inParms[i]['fixed'] = False
            
    # Function to evaluate the difference between the model and data.
    # This is minimised in the least-squared sense by the fitter
    def errFn(p, fjac=None):
        status = 0
        return status, (poly5(p)(xData) - yData) / dyData

    # Use mpfit to perform the fitting
    mp = mpfit(errFn, parinfo=inParms, quiet=True)
    return mp.params, mp.status, mp.fnorm/mp.dof
