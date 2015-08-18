#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     module_spec_extract.py                                            #
#                                                                             #
# PURPOSE:  Function to extract the spectra from FITS files for the PPC.      #
#                                                                             #
# REQUIRED: Requires the numpy, scipy and astropy modules.                    #
#                                                                             #
# MODIFIED: 22-July-2015 by C. Purcell                                        #
#                                                                             #
# TODO:                                                                       #
#           * Catch errors with extraction, fitting and propagate into        #
#             returned catalogue in a sensible way (flags?).                  #
#           * Return chi-sq from the Stokes I model fit                       #
#                                                                             #
#=============================================================================#

import os
import sys
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
def mod_spec_extract(catRec, fitsLstI, fitsLstQ, fitsLstU, freqArr_Hz,
                     sumBox_pix=3, polyOrd=2, outDataDir='./OUT',
                     nBeams=50.0, LF=None):
    """
    Extract on-source spectra and noise-spectra from Stokes I, Q and U FITS
    files for each entry in a catalogue. Files must contain single planes for
    one frequency (channel) and must be listed in frequency order.
    
    The spectra are summed over a square extraction box and normalised so that
    the spectrum extracted from an unresolved source would have the same
    amplitude as the peak spectrum.
    
    Each channel in the noise spectrum is measured using the MADFM statistic
    from a square aperture of 50 beam areas. Determing the aperture size
    requires the header keywords BMAJ and BMIN in the FITS files.
    
    This routine also returns a model Stokes I spectrum created by fitting a
    polynomal (<=5th order) to the Stokes I data. Requires the corresponding
    frequency vector to the Stokes I amplitudes.
    
    ARGUMENTS:
               catRec               ... record array catalogue table
                                        required: [uniqueName, x_deg, y_deg]
               fitsLstI             ... ordered list of Stokes I FITS files
               fitsLstQ             ... ordered list of Stokes Q FITS files
               fitsLstU             ... ordered list of Stokes U FITS files
               freqArr_Hz           ... ordered array of frequencies in Hz
               sumBox_pix           ... side of source box in pixels [3]
               polyOrd=2            ... order of polynomial fit [2]
               outDataDir           ... output directory for spectra [./OUT]
               nBeams               ... area (in beams) of noise aperture
               LF                   ... file-handle for log file [None]
    
    """

    # Default logging to STDOUT
    if LF is None:
        LF = sys.stdout

    # Read sample header from the first FITS file (WCS and beam parameters)
    try:
        head = pf.getheader(fitsLstQ[0])
        wcsDict = mkWCSDict(head)
        bmaj, bmin, bpa =  get_beam_from_header(head)
        log_wr(LF, 'Successfully read a sample FITS header (%s).' %
               fitsLstQ[0])
    except Exception:
        log_wr(LF, 'Err: Failed to read a sample FITS file header (%s).' %
                 fitsLstQ[0])        
        log_fail(LF, traceback.format_exc())
        
    # Calculate the normalisation to be applied to extracted spectra
    # Currently use geometric mean FWHM = sqrt(bmaj*bmin)
    try:
        beamFWHM_deg = m.sqrt(bmaj * bmin)
        beamFWHM_pix = beamFWHM_deg / wcsDict['pixscale']
        fNormSumbox = calc_sumbox_norm(beamFWHM_pix, sumBox_pix)
        log_wr(LF, '> Beam FWHM=%.2f pixels, sumBox=%d, fNorm=%.5f' % 
               (beamFWHM_pix, sumBox_pix, fNormSumbox))
    except Exception:
        log_wr(LF, 'Err: Failed to calculate beam/f_norm from header.')
        log_fail(LF, traceback.format_exc())

    # Determine the size of the cutout box to use
    beamArea_pix = get_beam_area(bmaj, bmin, wcsDict['pixscale'])
    halfSideNoise_pix = m.floor(m.sqrt(nBeams * beamArea_pix) / 2.0)
    halfSideSrc_pix = m.floor(sumBox_pix/2.0)

    # Create the output data directory
    if not os.path.exists(outDataDir):
        os.mkdir(outDataDir)

    # Check the catRec record array has the required fields
    reqCols = ['uniqueName','x_deg', 'y_deg']
    if not min([x in catRec.dtype.names for x in reqCols]):
        log_fail(LF, "Catalogue is missing a required field: %s" % reqCols)

    # Extract the I spectra from the list of FITS files
    cat = (catRec['uniqueName'], catRec['x_deg'], catRec['y_deg'])
    specIFileLst, rmsIFileLst = extract_spec_planes(cat,
                                                    fitsLstI,
                                                    freqArr_Hz,
                                                    halfSideSrc_pix,
                                                    halfSideNoise_pix,
                                                    fNormSumbox,
                                                    outDataDir,
                                                    'I',
                                                    False,
                                                    LF)
    
    # Extract the Q spectra from the list of FITS files
    specQFileLst, rmsQFileLst = extract_spec_planes(cat,
                                                    fitsLstQ,
                                                    freqArr_Hz,
                                                    halfSideSrc_pix,
                                                    halfSideNoise_pix,
                                                    fNormSumbox,
                                                    outDataDir,
                                                    'Q',
                                                    False,
                                                    LF)
    
    # Extract the U spectra from the list of FITS files
    specUFileLst, rmsUFileLst = extract_spec_planes(cat,
                                                    fitsLstU,
                                                    freqArr_Hz,
                                                    halfSideSrc_pix,
                                                    halfSideNoise_pix,
                                                    fNormSumbox,
                                                    outDataDir,
                                                    'U',
                                                    False,
                                                    LF)

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
             ('fitIstatus', 'i8'),
             ('extractStatus', 'i8')]
    specRec = np.zeros(len(catRec), dtype=dType)

    # Loop through the sources and measure the spectra
    log_wr(LF, "> Measuring the spectra ...")
    for i in range(len(specIFileLst)):
        log_wr(LF, " -> Processing '%s'..." % catRec[i]['uniqueName'])

        try:
            freqArr_Hz, specI = np.loadtxt(specIFileLst[i], unpack=True,
                                           dtype='f8')
            freqArr_Hz, rmsSpecI = np.loadtxt(rmsIFileLst[i], unpack=True,
                                              dtype='f8')
            freqArr_Hz, specQ = np.loadtxt(specQFileLst[i], unpack=True,
                                           dtype='f8')
            freqArr_Hz, rmsSpecQ = np.loadtxt(rmsQFileLst[i], unpack=True,
                                              dtype='f8')
            freqArr_Hz, specU = np.loadtxt(specUFileLst[i], unpack=True,
                                           dtype='f8')
            freqArr_Hz, rmsSpecU = np.loadtxt(rmsQFileLst[i], unpack=True,
                                              dtype='f8')
        except Exception:
            log_wr(LF, ">>> Warn: failed to load extracted spectrum.")
            specRec[i]['extractStatus'] = 0
            continue

        # Calculate the average RMS and PI spectrum
        rmsSpecQUAvg = (rmsSpecQ + rmsSpecU) / 2.0
        np.savetxt(outDataDir + '/' + catRec[i]['uniqueName'] +
                   '_rmsSpecQUavg.dat', zip(freqArr_Hz, rmsSpecQUAvg))
        
        # Fit a polynomial model to the Stokes I spectrum (order<=5)
        # Frequency axis must be in GHz to avoid overflow errors
        # Arrays must be float64
        try:
            p, fitIstatus = fit_spec_poly5(freqArr_Hz/1e9, specI, rmsSpecI,
                                           polyOrd)
            specImodel = poly5(p)(freqArr_Hz/1e9)
        except Exception:
            fitIstatus= 0
            
        if fitIstatus<1:
            log_wr(LF, "> Warn: I fit failed on '%s'." %
                   catRec[i]['uniqueName'] )
            specImodel = np.ones_like(specI)
            p = [0, 0, 0, 0, 0, 1]

        # Save the model I spectrum and coefficient to simple text files
        np.savetxt(outDataDir + '/' + catRec[i]['uniqueName'] +
                   '_specImodel.dat', zip(freqArr_Hz, specImodel))
        np.savetxt(outDataDir + '/' + catRec[i]['uniqueName'] +
                   '_polyCoeff.dat', p)
    
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
        specRec[i]['bmaj_deg'] = bmaj
        specRec[i]['bmin_deg'] = bmin
        specRec[i]['bpa_deg'] = bpa
        specRec[i]['pixscale_deg'] = wcsDict['pixscale']
        specRec[i]['fNormSumbox'] = fNormSumbox
        specRec[i]['fitIstatus'] = fitIstatus
        specRec[i]['coeffPolyIspec'] = ','.join([str(x) for x in p])
        specRec[i]['extractStatus'] = 1

    return specRec


#-----------------------------------------------------------------------------#
def extract_spec_planes(cat, fitsLst, freqArr_Hz, halfSideSrc_pix=0,
                        halfSideNoise_pix=10, fNormSumbox=1.0,
                        outDataDir='./', suffix='I',
                        saveAvgFits=False, LF=None):
    """
    This function extracts a spectrum from a list of 2D FITS images.
    Spectra are saved to ASCII files in the specified output directory and
    a list of files returned.
    """
    
    # dtype conversions
    halfSideSrc_pix = int(halfSideSrc_pix)
    halfSideNoise_pix = int(halfSideNoise_pix)
    
    # Read sample FITS header & and parse WCS information
    head0 = pf.getheader(fitsLst[0])
    head2D = strip_fits_dims(header=head0, minDim=2, forceCheckDims=5)
    wcs2D = pw.WCS(head2D)
    
    # Unpack the catalogue & convert to nearest pixel coords
    names, x_deg, y_deg = cat
    x_pix, y_pix = zip(*wcs2D.wcs_world2pix(zip(x_deg, y_deg), 0))
    x_pix = np.round(x_pix).astype('i4')
    y_pix = np.round(y_pix).astype('i4')

    # Determine the absolute pixel bounds of the noise box for each entry
    # Noise-box is offset if it overlaps an edge. Assumes box will
    # always be smaller than the image array
    lowLim = np.zeros_like(x_pix)
    highLim = np.zeros_like(x_pix) + head2D['NAXIS1'] -1
    xMinNoise_pix = x_pix - halfSideNoise_pix
    xMinNoise_pix = np.max(zip(xMinNoise_pix, lowLim), 1)
    xMaxNoise_pix = xMinNoise_pix + halfSideNoise_pix*2
    xMaxNoise_pix = np.where(xMaxNoise_pix>highLim,
                             highLim,
                             xMaxNoise_pix)
    xMinNoise_pix = np.where(xMaxNoise_pix>=highLim,
                             xMaxNoise_pix - halfSideNoise_pix*2,
                             xMinNoise_pix)
    noiseEdgeFlag = np.zeros_like(x_pix)
    noiseEdgeFlag = np.where(xMinNoise_pix<=lowLim, 1, noiseEdgeFlag)
    noiseEdgeFlag= np.where(xMaxNoise_pix>=highLim, 1, noiseEdgeFlag)

    highLim = np.zeros_like(y_pix) + head2D['NAXIS2'] -1
    yMinNoise_pix = y_pix - halfSideNoise_pix
    yMinNoise_pix = np.max(zip(yMinNoise_pix, lowLim), 1)
    yMaxNoise_pix = yMinNoise_pix + halfSideNoise_pix*2
    yMaxNoise_pix = np.where(yMaxNoise_pix>highLim,
                             highLim,
                             yMaxNoise_pix)
    yMinNoise_pix = np.where(yMaxNoise_pix>=highLim,
                             yMaxNoise_pix - halfSideNoise_pix*2,
                             yMinNoise_pix)
    noiseEdgeFlag = np.where(yMinNoise_pix<=lowLim, 1, noiseEdgeFlag)
    noiseEdgeFlag= np.where(yMaxNoise_pix>=highLim, 1, noiseEdgeFlag)
    
    # Determine the source-box bounds     
    xMinSrc_pix = x_pix - xMinNoise_pix - halfSideSrc_pix
    xMaxSrc_pix = x_pix - xMinNoise_pix + halfSideSrc_pix
    yMinSrc_pix = y_pix - yMinNoise_pix - halfSideSrc_pix
    yMaxSrc_pix = y_pix - yMinNoise_pix + halfSideSrc_pix
    
    # Set flags if the source box overlaps the edge
    srcEdgeFlag = np.zeros_like(x_pix)    
    lowLim = np.zeros_like(x_pix)
    highLim = np.zeros_like(x_pix) + head2D['NAXIS1'] -1
    srcEdgeFlag = np.where(xMinSrc_pix<=lowLim, 1, srcEdgeFlag)
    srcEdgeFlag= np.where(xMaxSrc_pix>=highLim, 1, srcEdgeFlag) 
    highLim = np.zeros_like(y_pix) + head2D['NAXIS2'] -1
    srcEdgeFlag = np.where(yMinSrc_pix<=lowLim, 1, srcEdgeFlag)
    srcEdgeFlag= np.where(yMaxSrc_pix>=highLim, 1, srcEdgeFlag) 

    # TODO:
    # If srcEdgeFlag skip source extraction
    
    # Root names of save files
    saveRootLst = map(''.join,zip([outDataDir + '/']*len(x_pix), names))
    
    # Loop through the FITS planes and extract each postage stamp
    for i in range(len(fitsLst)):

        log_wr(LF, "> Processing fits file '%s' ..." % fitsLst[i])
        
        # Open the FITS file in memory-mapping mode
        planeHDULst = pf.open(fitsLst[i], 'readonly', memmap=True)
        
        # Loop through the catalogue positions
        for j in range(len(x_pix)):
            log_wr(LF, "> Extracting '%s' ..." % names[j])

            # Open save files on first call
            if i==0:
                wMode = 'w'
                saveRootLst[j] = outDataDir + '/' + names[j]
            else:
                wMode = 'a'
                
            # Extract the noise data from the sub-cube
            try:
                data = planeHDULst[0].data
                nDim = len(data.shape)
                if nDim == 2:
                    dataSub = data[yMinNoise_pix[j]:yMaxNoise_pix[j]+1,
                                   xMinNoise_pix[j]:xMaxNoise_pix[j]+1]
                elif nDim == 3:
                    dataSub = data[0,
                                   yMinNoise_pix[j]:yMaxNoise_pix[j]+1,
                                   xMinNoise_pix[j]:xMaxNoise_pix[j]+1]
                elif nDim == 4:
                    dataSub = data[0, 0,
                                   yMinNoise_pix[j]:yMaxNoise_pix[j]+1,
                                   xMinNoise_pix[j]:xMaxNoise_pix[j]+1]
                else:
                    errStr = "Err: FITS file '%s' has %d dimensions." % \
                             (fitsLst[i], nDim) + \
                             'Only files with 2 - 4 dimensions supported.'
                    log_fail(LF, errStr)
                dataSpec = dataSub[yMinSrc_pix[j]:yMaxSrc_pix[j]+1,
                                   xMinSrc_pix[j]:xMaxSrc_pix[j]+1]
                rms = MAD(dataSub.flatten())                
                
                # Average the spectra in X and Y directions
                dataSpec = dataSpec.sum(-1).sum(-1)
                dataSpec *= fNormSumbox
                
            except Exception:
                log_wr(LF, ">>> Warn: failed to read FITS array.")
                log_wr(LF, traceback.format_exc())
                rms = None
                dataSpec = None
                dataSub = None

            # Save a sub-plane FITS file of the first channel
            if i==0 and dataSub is not None:
                headSub = head2D.copy()
                headSub['CRPIX1'] = head2D['CRPIX1'] - xMinNoise_pix[j]
                headSub['CRPIX2'] = head2D['CRPIX2'] - yMinNoise_pix[j]
                outFits = saveRootLst[j] + "_stamp" + suffix + ".fits"
                pf.writeto(outFits, dataSub, headSub, clobber=True)

            # Write the value for the current plane to disk
            # If extraction failed then write a NaN
            outSpecDat = saveRootLst[j] + '_spec' + suffix + '.dat'
            if not rms is None:
                open(outSpecDat, wMode).write("%e %e\n" % (freqArr_Hz[i],
                                                           dataSpec))
            else:
                open(outSpecDat, wMode).write("%e nan\n" % freqArr_Hz[i])
            outRMSDat = saveRootLst[j] + '_rmsSpec' + suffix + '.dat'
            if not dataSpec is None:
                open(outRMSDat, wMode).write("%e %e\n" % (freqArr_Hz[i], rms))
            else:
                open(outRMSDat, wMode).write("%e nan\n" % freqArr_Hz[i])

    # List of files containing the spectra
    specFileLst = map(''.join,zip(saveRootLst,
                                  ['_spec' + suffix + '.dat']*len(x_pix)))
    rmsFileLst = map(''.join,zip(saveRootLst,
                                 ['_rmsSpec' + suffix + '.dat']*len(x_pix)))

    return specFileLst, rmsFileLst
    

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
    coeffs = mp.params

    return mp.params, mp.status
