#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     module_spec_extract_area.py                                       #
#                                                                             #
# PURPOSE:  Extract emission regions from FITS files for POSSUM pipeline.     #
#                                                                             #
# REQUIRED: Requires the numpy, scipy and astropy modules.                    #
#                                                                             #
# MODIFIED: 24-November-2015 by C. Purcell                                    #
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
def mod_spec_extract_area(catRec, fitsLstI, fitsLstQ, fitsLstU, fitsLstV=None,
                          freqArr_Hz=None, extractMode="box", sumBox_pix=3,
                          gaussAperture_sigma=2.3548, polyOrd=2,
                          outDataDir='./OUT', nBeams=50.0, saveCubeMode=1,
                          doOverwrite=False, LF=None):
    """
    Run the procedure to extract data from each of the Stokes I, Q and U image
    planes, stored in FITS files. For each catalogue entry and Stokes component
    a FITS file is created containing a cutout cube (area=nBeams), an extraction
    aperture mask and a spectra table with columns 'freq', 'rms', 'src' and
    'model'. A polynomial of order <=5 is fit to the Stokes I spectrum and
    stored in the 'model' column. Measurements made on each source spectrum
    are written to a record array and returned.
    """

    # Default logging to STDOUT
    if LF is None:
        LF = sys.stdout

    # Read sample header from the first FITS file (WCS and beam parameters)
    # Determine the size of the postage stamp cutout box (N beam-areas)
    try:
        head = pf.getheader(fitsLstQ[0])
        wcsDict = mkWCSDict(head)
        bmaj_deg, bmin_deg, bpa_deg =  get_beam_from_header(head)
        log_wr(LF, 'Successfully read a sample FITS header (%s).' %
               fitsLstQ[0])        
        beamArea_deg = get_beam_area(bmaj_deg, bmin_deg, 1)
        cutoutSide_deg = m.sqrt(nBeams * beamArea_deg)
    except Exception:
        log_wr(LF, 'Err: Failed to read beam from sample header (%s).' %
               fitsLstQ[0])        
        log_fail(LF, traceback.format_exc())

    # Create the output data directory
    if not os.path.exists(outDataDir):
        os.mkdir(outDataDir)
        
    # Check the catRec record array has the required fields
    reqCols = ['uniqueName','x_deg', 'y_deg']
    if not min([x in catRec.dtype.names for x in reqCols]):
        errStr = "Catalogue passed to the extraction module is " + \
                 "missing a required field: %s" % reqCols
        log_fail(LF, errStr)

    # Perform the extraction on the Stokes I data
    cat = (catRec['uniqueName'], catRec['x_deg'], catRec['y_deg'])
    dataIFileLst = extract_data_from_planes(cat,
                                            fitsLstI,
                                            cutoutSide_deg,
                                            freqArr_Hz,
                                            extractMode,
                                            gaussAperture_sigma,
                                            sumBox_pix,
                                            outDataDir,
                                            'I',
                                            saveCubeMode=saveCubeMode,
                                            save1stPlane=True,
                                            doOverwrite=doOverwrite,
                                            LF=LF)
    dataQFileLst = extract_data_from_planes(cat,
                                            fitsLstQ,
                                            cutoutSide_deg,
                                            freqArr_Hz,
                                            extractMode,
                                            gaussAperture_sigma,
                                            sumBox_pix,
                                            outDataDir,
                                            'Q',
                                            saveCubeMode=saveCubeMode,
                                            doOverwrite=doOverwrite,
                                            LF=LF)
    dataUFileLst = extract_data_from_planes(cat,
                                            fitsLstU,
                                            cutoutSide_deg,
                                            freqArr_Hz,
                                            extractMode,
                                            gaussAperture_sigma,
                                            sumBox_pix,
                                            outDataDir,
                                            'U',
                                            saveCubeMode=saveCubeMode,
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
        specRec[i]['boxScale_pix'] = sumBox_pix
        specRec[i]['nFreqChan'] = len(freqArr_Hz)
        specRec[i]['deltaFreqChan_Hz'] = np.nanmin(np.diff(freqArr_Hz))
        specRec[i]['bmaj_deg'] = bmaj_deg
        specRec[i]['bmin_deg'] = bmin_deg
        specRec[i]['bpa_deg'] = bpa_deg
        specRec[i]['pixscale_deg'] = wcsDict['pixscale']
        specRec[i]['fNormSumbox'] = HDUILst[0].header["fNorm"]
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
def extract_data_from_planes(cat, fitsLst, cutoutSide_deg, freqArr_Hz=None, 
                             extractMode="ellipse", gaussAperture_sigma=2.3548,
                             sumBox_pix=3, outDataDir='./', suffix='I',
                             saveCubeMode=2, save1stPlane=False,
                             doOverwrite=False, LF=None):
    """
    Extract data from a list of FITS planes of identical size and save them to
    disk. The data extracted are:
    * A one dimentional spectrum from the source:
     - if 'sumBox_pix': summed and normalised assuming an unresolved source.
     - if 'gaussAperture_sigma': summed and divided by the beam area.
    * A spectrum of the RMS noise in each plane.
    Optional:
    * A cube centred on each source, or offset if abutting an edge.
    * A single-plane mask image showing which pixels the source spectra have
      been extracted from.
    Two different source extraction aperture may be used:
      - a 'sumBox_pix' square grid of pixels (typically 3x3)
      - an elliptical beam aperture of size 'gaussAperture_sigma'
        -- i.e., the FWHM beam is gaussAperture_sigma=2.3548

    All data is saved as a single FITS file in the following extensions:
    EXT 1: The cube covering the whole cutout (saveCubeMode=1), or just the
           source aperture (saveCubeMode=2), or save only the 1st plane
           (saveCubeMode=3) or blank (saveCubeMode=0).
    EXT 2: A mask image of the source extraction aperture.
    EXT 3: A FITS table containing the columns [frequency_Hz, src, rms].
    
    """

    # Default to channel numbers if no frequency array has been provided
    nChan = len(fitsLst)
    if freqArr_Hz is None:
        freqArr_Hz = np.arange(nChan, dtype="f4") + 1
    if len(freqArr_Hz)!=nChan:
        msg = ">>> Err: Number of planes and frequency channels do not match!"
        log_fail(LF, msg)
                
    # Read sample FITS header and parse WCS information
    head = pf.getheader(fitsLst[0])
    head2D = strip_fits_dims(header=head, minDim=2, forceCheckDims=5)
    wcsDict = mkWCSDict(head)
    wcs2D = pw.WCS(head2D)

    # Parse beam information and calculate area in pixels
    bmaj_deg, bmin_deg, bpa_deg = get_beam_from_header(head)
    beamArea_pix = get_beam_area(bmaj_deg, bmin_deg, wcsDict['pixscale'])
    log_wr(LF, '> Beam major=%.2f", Beam minor=%.2f", Beam area=%.2f px' %
           (bmaj_deg*3600.0, bmin_deg*3600.0, beamArea_pix))
    
    # Unpack the catalogue, convert to pixel coords and calculate
    # positions rounded to the nearest pixel centre
    nameArr, xArr_deg, yArr_deg = cat
    xArr_pix, yArr_pix = zip(*wcs2D.wcs_world2pix(zip(xArr_deg, yArr_deg), 0))
    xRndArr_pix = np.round(xArr_pix).astype('i4')
    yRndArr_pix = np.round(yArr_pix).astype('i4')
    
    # Determine the absolute pixel bounds of the cutout box for each entry
    log_wr(LF, "> Determining cutout bounds ...")
    cutoutRadius_pix =  m.floor(cutoutSide_deg / (2.0*wcsDict['pixscale']))
    [xMinCutArr_pix, xMaxCutArr_pix, yMinCutArr_pix, yMaxCutArr_pix,
     cutoutEdgeFlag] = calc_cutout_bounds(head2D['NAXIS1'],
                                          head2D['NAXIS2'],
                                          xRndArr_pix,
                                          yRndArr_pix,
                                          cutoutRadius_pix)
    
    # If extracting spectra under the elliptical beam footprint
    if extractMode=="ellipse":
        log_wr(LF, "> Extraction mode = 'ellipse'")
        
        # Use the beam footprint as a source extraction aperture (fill the
        # catalogue with beam ellipses). Default makes ellipse cover the FWHM
        gfac = 2.0 * m.sqrt(2.0 * m.log(2.0))  # sigma = FWHM/gfactor
        minAxArr_deg = (np.ones_like(xArr_deg) * bmin_deg
                        / gfac * gaussAperture_sigma)
        majAxArr_deg = (np.ones_like(xArr_deg) * bmaj_deg
                        / gfac * gaussAperture_sigma)
        paArr_deg = np.ones_like(xArr_deg) * bpa_deg
        minAxArr_pix = minAxArr_deg / wcsDict['pixscale']
        majAxArr_pix = majAxArr_deg / wcsDict['pixscale']
        paArr_rad = np.radians(paArr_deg)

        # Calculate indices of the pixels under the Gaussian fit footprints
        # These are later used to mask the small cutout of the source pixels
        # saved to a FITS file on disk
        log_wr(LF, "> Determining extraction mask indices ...")
        mskIndxLst = get_ell_indices(head2D["NAXIS1"],
                                     head2D["NAXIS2"],
                                     xArr_pix,
                                     yArr_pix,
                                     minAxArr_pix/2.0,
                                     majAxArr_pix/2.0,
                                     paArr_rad)
        
        # Determine the absolute pixel bounds of the source extraction aperture
        xMinSrcArr_pix = np.zeros_like(xRndArr_pix)
        xMaxSrcArr_pix = np.zeros_like(xRndArr_pix)
        yMinSrcArr_pix = np.zeros_like(xRndArr_pix)
        yMaxSrcArr_pix = np.zeros_like(xRndArr_pix)
        for i in range(len(mskIndxLst)):
            xiArr = mskIndxLst[i][1]
            yiArr = mskIndxLst[i][0]
            xMinSrcArr_pix[i] = np.min(xiArr)
            xMaxSrcArr_pix[i] = np.max(xiArr)
            yMinSrcArr_pix[i] = np.min(yiArr)
            yMaxSrcArr_pix[i] = np.max(yiArr)
        
        # The normalisation is the beam area in pixels
        fNorm =1.0/beamArea_pix

    # If extracting spectra from a square box [default]
    else:
        log_wr(LF, "> Extraction mode = 'box'")
        
        # Determine the absolute pixel bounds of the source extraction box
        halfSideSrc_pix = int(m.floor(sumBox_pix/2.0))
        xMinSrcArr_pix = xRndArr_pix - halfSideSrc_pix
        xMaxSrcArr_pix = xRndArr_pix + halfSideSrc_pix
        yMinSrcArr_pix = yRndArr_pix - halfSideSrc_pix
        yMaxSrcArr_pix = yRndArr_pix + halfSideSrc_pix
        
        # Calculate the indices of the pixels under the extraction box
        log_wr(LF, "> Determining extraction mask indixes ...")
        mskIndxLst = []
        indArr = np.where(np.ones((sumBox_pix, sumBox_pix)))
        for i  in range(len(xRndArr_pix)):
            xiArr = indArr[-1] + xMinSrcArr_pix[i]
            yiArr = indArr[-2] + yMinSrcArr_pix[i]
            mskIndxLst.append((yiArr, xiArr))
        
        # Calculate the normalisation to be applied to box-extracted spectra
        beamFWHM_deg = m.sqrt(bmaj_deg * bmin_deg)          # Geometric mean
        beamFWHM_pix = beamFWHM_deg / wcsDict['pixscale']
        fNorm = calc_sumbox_norm(beamFWHM_pix, sumBox_pix)
        log_wr(LF, '> Beam FWHM=%.2f pixels, sumBox=%d, fNorm=%.5f' %
               (beamFWHM_pix, sumBox_pix, fNorm))

    # Calculate source mask indices with respect to the cutout
    mskIndxWrtCutLst = []
    for j in range(len(mskIndxLst)):
        xIndxArr = mskIndxLst[j][-1] - xMinCutArr_pix[j]
        xIndxArr = xIndxArr.astype(int)
        yIndxArr = mskIndxLst[j][-2] - yMinCutArr_pix[j]
        yIndxArr = yIndxArr.astype(int)
        mskIndxWrtCutLst.append((yIndxArr, xIndxArr))

    # ... and with respect to the source box
    mskIndxWrtSrcLst = []
    for j in range(len(mskIndxLst)):
        xIndxArr = mskIndxLst[j][-1] - xMinSrcArr_pix[j]
        xIndxArr = xIndxArr.astype(int)
        yIndxArr = mskIndxLst[j][-2] - yMinSrcArr_pix[j]
        yIndxArr = yIndxArr.astype(int)
        mskIndxWrtSrcLst.append((yIndxArr, xIndxArr))

    # Set the pixel limits depending on the saveCubeMode
    if saveCubeMode==1:        # Save whole cutout
        mskIndxSaveLst = mskIndxWrtCutLst
        xMinSaveArr_pix = xMinCutArr_pix
        xMaxSaveArr_pix = xMaxCutArr_pix
        yMinSaveArr_pix = yMinCutArr_pix
        yMaxSaveArr_pix = yMaxCutArr_pix
        
    else:                      # Save only the source cutout
        mskIndxSaveLst = mskIndxWrtSrcLst
        xMinSaveArr_pix = xMinSrcArr_pix
        xMaxSaveArr_pix = xMaxSrcArr_pix
        yMinSaveArr_pix = yMinSrcArr_pix
        yMaxSaveArr_pix = yMaxSrcArr_pix

    # Create a list of root path names for output files
    saveRootLst = map(''.join, zip([outDataDir + '/']*len(xArr_pix), nameArr))
    skipLst = [False] * len(xArr_pix)

    # Loop through the list of FITS file (planes = frequency channels)
    for i in range(nChan):
        log_wr(LF, "> Reading FITS file '%s' ..." % fitsLst[i])
        planeHDULst = pf.open(fitsLst[i], 'readonly', memmap=True)
        
        # Loop through the catalogue positions
        for j in range(len(xArr_pix)):
            log_wr(LF, "> Extracting '%s' ..." % nameArr[j])
            outFits = saveRootLst[j] + "_spec" + suffix + ".fits"
            if i==0 and not doOverwrite and os.path.exists(outFits):
                skipLst[j] = True
            if skipLst[j]:
                print ">>> Older version exists: skipping ..."
                continue

            # Extract the cutout data from the sub-cube
            try:
                data = planeHDULst[0].data
                nDim = len(data.shape)
                if nDim == 2:
                    dataCut = data[yMinCutArr_pix[j]:yMaxCutArr_pix[j]+1,
                                   xMinCutArr_pix[j]:xMaxCutArr_pix[j]+1]
                elif nDim == 3:
                    dataCut = data[0,
                                   yMinCutArr_pix[j]:yMaxCutArr_pix[j]+1,
                                   xMinCutArr_pix[j]:xMaxCutArr_pix[j]+1]
                elif nDim == 4:
                    dataCut = data[0, 0,
                                   yMinCutArr_pix[j]:yMaxCutArr_pix[j]+1,
                                   xMinCutArr_pix[j]:xMaxCutArr_pix[j]+1]
                else:
                    errStr = "Err: FITS file '%s' has %d dimensions." % \
                             (fitsLst[i], nDim) + \
                             "Only files with 2 - 4 dimensions supported."
                    log_fail(LF, errStr)

            except Exception:
                log_wr(LF, ">>> Warn: failed to read FITS array.")
                log_wr(LF, traceback.format_exc())
                dataCut = None

            # If this is the first channel create the empty FITS file
            if i==0:

                # Create a zero-filled cube to store the full cutout or
                # the smaller cube containing the source
                shape = (nChan,
                         yMaxSaveArr_pix[j] - yMinSaveArr_pix[j] +1,
                         xMaxSaveArr_pix[j] - xMinSaveArr_pix[j] +1)
                cubeSpec = np.zeros(shape, dtype="f4")
                headCut = head2D.copy()
                headCut["DATAMIN"] = 0.0
                headCut["DATAMAX"] = 0.01
                headCut["fNorm"] = fNorm
                headCut['CRPIX1'] = head2D['CRPIX1'] - xMinSaveArr_pix[j]
                headCut['CRPIX2'] = head2D['CRPIX2'] - yMinSaveArr_pix[j]
                hdu0 = pf.PrimaryHDU(cubeSpec, headCut)

                # Create a mask image to show the source pixels extracted
                shape = (yMaxSaveArr_pix[j] - yMinSaveArr_pix[j] +1,
                         xMaxSaveArr_pix[j] - xMinSaveArr_pix[j] +1)
                maskData = np.zeros(shape, dtype="f4")
                maskData[mskIndxSaveLst[j]] = 1.0
                hdu1 = pf.ImageHDU(maskData, headCut)

                # Create a record array table to store the summary spectrum
                col1 = pf.Column(name="freq", format="f4", 
                                 array=freqArr_Hz)
                col2 = pf.Column(name="rms", format="f4", 
                                 array=np.zeros_like(freqArr_Hz))
                col3 = pf.Column(name="src", format="f4", 
                                 array=np.zeros_like(freqArr_Hz))
                col4 = pf.Column(name="model", format="f4", 
                                 array=np.zeros_like(freqArr_Hz))
                hdu2 = pf.new_table([col1, col2, col3, col4])
                
                # Write the file to disk and clean up
                hduLst = pf.HDUList([hdu0, hdu1, hdu2])
                hduLst.writeto(outFits, output_verify="fix", clobber=True)
                hduLst.close()
            
            # Measure the noise in the current plane, excluding the source
            dataCutMasked = dataCut.copy()
            dataCutMasked[mskIndxWrtCutLst[j]] = np.nan
            dataRMS = MAD(dataCutMasked.flatten())
            
            # Sum the source spectra in X and Y directions
            dataSum = np.nansum(dataCut[mskIndxSaveLst[j]]) * fNorm
            
            # Open the FITS file and write the extracted spectrum and rms
            outHDULst = pf.open(outFits, "update", memmap=True) 
            outHDULst[2].data["src"][i] = dataSum
            outHDULst[2].data["rms"][i] = dataRMS               
            if saveCubeMode>0:
                xMin_pix = xMinSaveArr_pix[j] - xMinCutArr_pix[j]
                xMax_pix = xMaxSaveArr_pix[j] - xMinCutArr_pix[j]
                yMin_pix = yMinSaveArr_pix[j] - yMinCutArr_pix[j]
                yMax_pix = yMaxSaveArr_pix[j] - yMinCutArr_pix[j]
                dataSave = dataCut[yMin_pix:yMax_pix+1, xMin_pix:xMax_pix+1]
                outHDULst[0].data[i,:,:] = dataSave            
            outHDULst.close()
            del outHDULst

    # List of files containing the spectra
    outFileLst = map(''.join,zip(saveRootLst,
                                ['_spec' + suffix + '.fits']*len(xRndArr_pix)))

    return outFileLst


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
