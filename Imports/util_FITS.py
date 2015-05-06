#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     util_FITS.py                                                      #
#                                                                             #
# PURPOSE:  Utility functions to operate on FITS data.                        #
#                                                                             #
# MODIFIED: 12-Sep-2014 by C. Purcell                                         #
#                                                                             #
# CONTENTS:                                                                   #
#                                                                             #
#  mkWCSDict         ... parse fits header using wcslib                       #
#  world2pix         ... simple linear conversion world->pix                  #
#  pix2world         ... simple linear conversion pix->world                  #
#  chan2world        ... simple linear conversion chan->world                 #
#  world2chan        ... simple linear conversion world->chan                 #
#  strip_fits_dims   ... strip header and / or data dimensions                #
#  get_beam_from_header ... fetch the beam parameters from a FITS header      #
#  get_beam_area     ... calculate the effective beam area in px              #
#  get_subfits       ... cut a sub portion from a FITS cube                   #
#                                                                             #
#=============================================================================#

import os
import sys
import shutil
import re
import math as m
import numpy as np
import astropy.io.fits as pf
import astropy.wcs.wcs as pw


#-----------------------------------------------------------------------------#
def mkWCSDict(header, wrapX180=False, forceCheckDims=5):
    """
    Parse and calculate key WCS parameters using python WCSLIB.
    """
    w = {}
    
    header2D = strip_fits_dims(header=header, minDim=2,
                               forceCheckDims=forceCheckDims)
    wcs2D = pw.WCS(header2D)
    w['header2D'] = header2D
    if header['NAXIS']>=3:
        header3D = strip_fits_dims(header=header, minDim=3,
                                   forceCheckDims=forceCheckDims)
        wcs3D = pw.WCS(header3D,)
        w['header3D'] = header3D

    # Header values
    w['xrpix']  = header['CRPIX1']
    w['xrval']  = header['CRVAL1']
    w['xnaxis'] = header['NAXIS1']
    w['yrpix']  = header['CRPIX2']
    w['yrval']  = header['CRVAL2']
    w['ynaxis'] = header['NAXIS2']
    if header['NAXIS']>=3:
        w['zrpix']  = header['CRPIX3']
        w['zrval']  = header['CRVAL3']
        w['znaxis'] = header['NAXIS3']

    # Calculate the centre coordinates
    xCent_pix = float(header['NAXIS1']) / 2.0 + 0.5
    yCent_pix = float(header['NAXIS2']) / 2.0 + 0.5
    if header['NAXIS']>=3:
        zCent_pix = float(header['NAXIS3']) /2.0 + 0.5
        [[w['xcent'], w['ycent'], w['zcent']]] = \
               wcs3D.wcs_pix2world([(xCent_pix, yCent_pix, zCent_pix)], 1)
        
    else:
        [[w['xcent'], w['ycent']]] = \
                       wcs2D.wcs_pix2world([(xCent_pix,yCent_pix)],1)
    if wrapX180:
        if w['xcent']>180.0:
            w['xcent'] = w['xcent'] - 360.0
    
    # Calculate the image bounds in world coords
    try:
        a = wcs2D.calc_footprint()
    except Exception:
        # Deprecated version of the function
        a = wcs2D.calcFootprint()
    w['xmin'] = np.hsplit(a, 2)[0].min()
    w['xmax'] = np.hsplit(a, 2)[0].max()
    w['ymin'] = np.hsplit(a, 2)[1].min()
    w['ymax'] = np.hsplit(a, 2)[1].max()
    if header['NAXIS']>=3:
        [[dummy1, dummy2, z1]] = \
             wcs3D.wcs_pix2world([(xCent_pix, yCent_pix, 1.0)], 1)
        [[dummy1, dummy2, z2]] = \
             wcs3D.wcs_pix2world([(xCent_pix, yCent_pix, header['NAXIS3'])], 1)
        w['zmin'] = min([z1, z2])
        w['zmax'] = max([z1, z2])
    if wrapX180:
        if w['xmin']>180.0:
            w['xmin'] = w['xmin'] - 360.0
        if w['xmax']>180.0:
            w['xmax'] = w['xmax'] - 360.0
    
    # Set the type of position coordinate axes
    if wcs2D.wcs.lngtyp == 'RA':
        w['coord_type'] = 'EQU'
    elif wcs2D.wcs.lngtyp == 'GLON':
        w['coord_type'] = 'GAL'
    else:
        w['coord_type'] = ''
        
    # Determine the pixel scale at the refpix
    if header['NAXIS']>=3:
        crpix = np.array(wcs3D.wcs.crpix)
        [[x1, y1, z1]] = wcs3D.wcs_pix2world([crpix], 1)
        crpix += 1.0
        [[x2, y2, z2]] = wcs3D.wcs_pix2world([crpix], 1)
        w['zdelt'] = z2 - z1
    else:
        crpix = np.array(wcs2D.wcs.crpix)
        [[x1, y1]] = wcs2D.wcs_pix2world([crpix], 1)
        crpix -= 1.0
        [[x2, y2]] = wcs2D.wcs_pix2world([crpix], 1)
    cosy = m.cos( m.radians((y1 + y2) / 2.0) )
    w['xdelt'] = (x2 - x1) * cosy
    w['ydelt'] = y2 - y1
    w['pixscale'] = abs( (w['xdelt']+w['xdelt']) / 2.0 )

    return w
    

#-----------------------------------------------------------------------------#
def world2pix(w, x_deg, y_deg):
    """
    Simple linear conversion from deg to pixel (1-based pix numbering).
    """
    
    x_deg = float(x_deg)
    y_deg = float(y_deg)    
    cosF = m.cos( m.radians(y_deg) )
    x_pix = w['xrpix'] + (x_deg - w['xrval']) / (w['xdelt'] / cosF )
    y_pix = w['yrpix'] + (y_deg - w['yrval']) / w['ydelt'] 
    
    return x_pix, y_pix


#-----------------------------------------------------------------------------#
def pix2world(w, x_pix, y_pix):
    """
    Simple linear conversion from pixel to deg (1-based pix numbering).
    """
    
    x_pix = float(x_pix)
    y_pix = float(y_pix)   
    y_deg = w['yrval'] + ( (y_pix - w['yrpix']) * w['ydelt'] )
    cosF  = m.cos( m.radians(y_deg) )
    x_deg = w['xrval'] + ( (x_pix - w['xrpix']) * w['xdelt'] / cosF)

    return x_deg, y_deg


#-----------------------------------------------------------------------------#
def chan2world(w, z_pix):
    """
    Simple linear conversion from channel to world coordinates (1-based).
    """
    
    z_deg = w['zrval'] + ( (z_pix - w['zrpix']) * w['zdelt'] )

    return z_deg


#-----------------------------------------------------------------------------#
def world2chan(w, z_deg):
    """
    Simple linear conversion from channel to world coordinates (1-based).
    """
    
    z_pix = w['zrpix'] + ( (z_deg - w['zrval']) / w['zdelt'] )
    
    return z_pix


#-----------------------------------------------------------------------------#
def strip_fits_dims(data=None, header=None, minDim=2, forceCheckDims=0):
    """
    Strip array and / or header dimensions from a FITS data-array or header.
    """
    
    xydata = None

    # Strip unused dimensions from the header
    if not data is None:
        
        naxis = len(data.shape)
        extraDims = naxis - minDim
        if extraDims < 0:
            print "Too few dimensions in data. "
            sys.exit(1)

        # Slice the data to strip the extra dims
        if extraDims == 0:
            xydata = data.copy()
        elif extraDims == 1:
            xydata = data[0].copy()
        elif extraDims == 2:
            xydata = data[0][0].copy()
        elif extraDims == 3:
            xydata = data[0][0][0].copy()
        else:
            print "Data array contains %s axes" % naxis 
            print "This script supports up to 5 axes only."
            sys.exit(1)
        del data
        
    # Strip unused dimensions from the header
    if not header is None:
        
        header = header.copy()
        naxis = header['NAXIS']
        
        stripHeadKeys = ['NAXIS','CRVAL','CRPIX','CDELT','CTYPE','CROTA',
                         'CD1_','CD2_', 'CUNIT']

        # Force a check on all relevant keywords
        if forceCheckDims>0:
            for key in stripHeadKeys:
                for i in range(forceCheckDims+1):
                    if key+str(i) in header:
                        if i > naxis:
                            naxis = i
        
        extraDims = naxis - minDim
        if extraDims < 0:
            print "Too few dimensions in data. "
            sys.exit(1)

        # Delete the entries
        for i in range(minDim+1,naxis+1):
            for key in stripHeadKeys:
                if key+str(i) in header:
                    try:
                        del header[key+str(i)]
                    except Exception:
                        pass
        header['NAXIS'] = minDim
        header['WCSAXES'] = minDim

    # Return the relevant object(s)
    if not xydata is None and not header is None:
        return [xydata,header]
    elif not xydata is None and header is None:
        return xydata
    elif xydata is None and not header is None:
        return header
    else:
        print "Both header and data are 'Nonetype'."
        sys.exit(1)


#-----------------------------------------------------------------------------#
def get_beam_from_header(header):
    """
    Get the beam FWHM and PA from the FITS header. Will return 'None' for
    missing cards.
    """
    
    # First try standard header cards
    try:
        bmaj = header['BMAJ']
    except Exception:
        bmaj = None
    try:
        bmin = header['BMIN']
    except Exception:
        bmin = None
    try:
        bpa = header['BPA']
    except Exception:
        bpa = None

    # Or try history for AIPS CLEAN Beam
    beamHistStr = 'AIPS\s+CLEAN\sBMAJ=\s+(\S+)\s+BMIN=\s+(\S+)\s+BPA=\s+(\S+)'
    bmHistPat = re.compile(beamHistStr)
    if bmaj is None or bmin is None:
        try:
            history = header.get_history()
            history.reverse()
            for i in range(len(history)):
                mch = bmHistPat.match(history[i])
                if mch:
                    bmaj = float(mch.group(1))
                    bmin = float(mch.group(2))
                    bpa  = float(mch.group(3))
                    break
        except Exception:
            pass

    return bmaj, bmin, bpa


#-----------------------------------------------------------------------------#
def get_beam_area(bmaj, bmin, pixscale):
    """
    Calculate the effective beam area in pixels given the major and minor
    Gaussian beam FWHM and the pixel-scale.
    """
    
    # Convert beam FWHMs to Gaussian sigma
    gfactor = 2.0 * m.sqrt( 2.0 * m.log(2.0) )  # FWHM = gfactor*sigma
    sigmaMaj_deg = bmaj / gfactor
    sigmaMin_deg = bmin / gfactor
    
    # Calculate the beam area in the input units and in pixels
    beamArea_degsq = 2.0 * m.pi * sigmaMaj_deg * sigmaMin_deg
    pixArea_degsq = abs(pixscale**2.0)
    beamArea_pix = beamArea_degsq / pixArea_degsq
    
    return beamArea_pix


#-----------------------------------------------------------------------------#
def get_subfits(inFileName, x_deg, y_deg, radius_deg, zMin_w=None, zMax_w=None,
                do3D=False, wrapX180=False):
    """
    Cut out a sub-section of a FITS file (cube or map).
    """

    # Enforce the sign convention whereby the coordinates range from
    # 180 to -180 degrees, rather than 0 to 360.This is necessary for files
    # which allow negative Galactic coordinates.
    if wrapX180 and  x_deg>180.0:
        x_deg -= 360.0

    # Some sanity checks
    if not os.path.exists(inFileName): 
        return None, None

    # Read the FITS file
    hduLst = pf.open(inFileName, 'readonly', memmap=True)       
    w = mkWCSDict(hduLst[0].header, wrapX180, 5)
    wcs2D = pw.WCS(w['header2D'])
    
    # Find the pixel at the center of the new sub-image
    [ [xCent_pix, yCent_pix] ] = wcs2D.wcs_world2pix([(x_deg, y_deg)], 0)
    xCentInt_pix = round(xCent_pix)
    yCentInt_pix = round(yCent_pix)

    # Determine the X and Y bounds of the subarray within the main array
    radX_pix = m.ceil( radius_deg / abs(w['xdelt']) )
    radY_pix = m.ceil( radius_deg / abs(w['ydelt']) )
    xMin_pix = max(xCentInt_pix - radX_pix, 0)
    xMax_pix = min(xCentInt_pix + radX_pix, w['header2D']['NAXIS1']-1 )
    yMin_pix = max(yCentInt_pix - radY_pix, 0)
    yMax_pix = min(yCentInt_pix + radY_pix, w['header2D']['NAXIS2']-1 )
    zMin_chan = None
    zMax_chan = None
    if zMin_w and zMax_w and w['header3D']['NAXIS']>2:
        zMin_chan = max(m.ceil( world2chan(w, zMin_w) ), 0)
        zMax_chan = min(m.floor( world2chan(w, zMax_w) ) ,
                                w['header3D']['NAXIS3']-1)

    # Shift the refpix to reflect the new offset
    headSub = hduLst[0].header.copy()
    headSub['CRPIX1'] = hduLst[0].header['CRPIX1'] - xMin_pix
    headSub['CRPIX2'] = hduLst[0].header['CRPIX2'] - yMin_pix
    if zMin_chan and zMax_chan:
        headSub['CRPIX3'] = hduLst[0].header['CRPIX3'] - zMin_chan
      
    # Slice the data
    nAxis = len(hduLst[0].data.shape)
    if zMin_chan and zMax_chan:
        if nAxis == 2:
            dataSub = hduLst[0].data[yMin_pix:yMax_pix+1,
                                     xMin_pix:xMax_pix+1]
        elif nAxis == 3:
            dataSub = hduLst[0].data[zMin_chan:zMax_chan+1,
                                     yMin_pix:yMax_pix+1,
                                     xMin_pix:xMax_pix+1]
        elif nAxis == 4: 
            dataSub = hduLst[0].data[:,
                                     zMin_chan:zMax_chan+1,
                                     yMin_pix:yMax_pix+1,
                                     xMin_pix:xMax_pix+1]
        else:
            return None, None
    else:
        if nAxis == 2:
            dataSub = hduLst[0].data[yMin_pix:yMax_pix+1,
                                     xMin_pix:xMax_pix+1]
        elif nAxis == 3:
            dataSub = hduLst[0].data[:,
                                     yMin_pix:yMax_pix+1,
                                     xMin_pix:xMax_pix+1]
        elif nAxis == 4: 
            dataSub = hduLst[0].data[:,
                                     :,
                                     yMin_pix:yMax_pix+1,
                                     xMin_pix:xMax_pix+1]
        else:
            return None, None

    # Free the memory by closing the file
    hduLst.close()

    # Update the sub header
    try:
        headSub["DATAMIN"] = np.nanmin(dataSub.flatten())
        headSub["DATAMAX"] = np.nanmax(dataSub.flatten())
        headSub["NAXIS1"] = dataSub.shape[-1]
        headSub["NAXIS2"] = dataSub.shape[-2]
        if nAxis > 2:
            headSub["NAXIS3"] = dataSub.shape[-3]
    except Exception:
        pass
    
    return dataSub, headSub
