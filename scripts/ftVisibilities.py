#!/usr/bin/env python
"""
Perform a Fourier Transform (standard or fast) on LOFAR ACC/XST data or widefield MS data (e.g. PAPER) to form a complex or Stokes dirty image dirty image, single files only
"""

import numpy as np
from matplotlib import pyplot as plt
import datetime
import ephem
import sys,os
import SWHT
try:
    import casacore.tables as tbls
except ImportError:
    print 'Warning: could not import casacore.tables, will not be able to read measurement sets'

#import scipy.constants
#cc = scipy.constants.c
cc = 299792458.0 #speed of light, m/s

if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] ACC/XST/MS FILE')
    o.set_description(__doc__)
    o.add_option('--station', dest='station', default=None,
        help = 'LOFAR ONLY: station name, e.g. SE607, if this is used then the ant_field and ant_array options are not required, default: None')
    o.add_option('-F', '--ant_field', dest='ant_field', default=None,
        help = 'LOFAR ONLY: AntennaField.conf file for the LOFAR station of the ACC files, default: None')
    o.add_option('-A', '--ant_array', dest='ant_array', default=None,
        help = 'LOFAR ONLY: AntennaArray.conf file for the LOFAR station geographical coordinates, default: None')
    o.add_option('-D', '--deltas', dest='deltas', default=None,
        help = 'LOFAR ONLY: iHBADeltas.conf file, only required for HBA imaging, default: None')
    o.add_option('-r', '--rcumode', dest='rcumode', default=3, type='int',
        help = 'LOFAR ONLY: Station RCU Mode, usually 3,5,6,7, for XST it will override filename metadata default: 3(LBA High)')
    o.add_option('-s', '--subband', dest='subband', default='0',
        help = 'Select which subband(s) to image, for ACC and MS it will select, for multiple subbands use X,Y,Z and for range use X_Y notation, for XST it will override filename metadata, default:0')
    o.add_option('-p', '--pixels', dest='pixels', default=64, type='int',
        help = 'Width of image in pixels, default: 64')
    o.add_option('-C', '--cal', dest='calfile', default=None,
        help = 'LOFAR ONLY: Apply a calibration soultion file to the data.')
    o.add_option('-S', '--save', dest='savefig', default=None,
        help = 'Save the figure using this name, type is determined by extension')
    o.add_option('--conv', dest='conv', default='fast',
        help = 'If using FFT, choose a convolution function: fast(nearest neighbor), rectangle, gaussian, prolate spheroid. default:fast')
    o.add_option('--dft', dest='dft', action='store_true',
        help = 'Form image with a direct FT instead of an FFT')
    o.add_option('--nodisplay', dest='nodisplay', action='store_true',
        help = 'Do not display the generated image')
    o.add_option('--pkl', dest='pkl', default=None,
        help = 'Save complex images in a numpy array in a pickle file using this name (include .pkl extention), default: tempImage.pkl')
    o.add_option('-i', '--int', dest='int_time', default=1., type='float',
        help = 'LOFAR ONLY: Integration time, used for accurate zenith pointing, for XST it will override filename metadata, default: 1 second')
    o.add_option('-c', '--column', dest='column', default='CORRECTED_DATA', type='str',
        help = 'MS ONLY: select which data column to image, default: CORRECTED_DATA')
    o.add_option('--override', dest='override', action='store_true',
        help = 'LOFAR XST ONLY: override filename metadata for RCU, integration length, and subband')
    o.add_option('--autos', dest='autos', action='store_true',
        help = 'Include the auto-correlation in the image, by default they are blanked')
    o.add_option('--weight', dest='weighting', default='natural',
        help = 'Weighting mode, natural (default), uniform')
    o.add_option('--fov', dest='fov', default=180., type='float',
        help = 'Field of View in degrees, default: 180 (all-sky)')
    o.add_option('--psf', dest='psf', action='store_true',
        help='Plot the PSF instead of the image')
    o.add_option('-t', '--times', dest='times', default='0',
        help = 'KAIRAXST ONLY: Select which integration(s) to image, can use a[seconds] to average, d[step size] to decimate, of a specific range of integrations similar to the subband selection option, default:0 (select the first integration of the file)')
    o.add_option('--uvplot', dest='uvplot', action='store_true',
        help='Display a 2D UV coverage/sampling plot of the projected baselines')
    opts, args = o.parse_args(sys.argv[1:])

    # parse subbands
    sbs = np.array(SWHT.util.convert_arg_range(opts.subband))

    # setup variables for combined visibilities and uvw samples
    visComb = np.array([]).reshape(4, 0, len(sbs))
    uvwComb = np.array([]).reshape(0, 3, len(sbs))

    dataFmt = None
    if (not (opts.station is None)) or (not (opts.ant_field is None)): # If using LOFAR data, get station information
        lofarStation = SWHT.lofarConfig.getLofarStation(name=opts.station, affn=opts.ant_field, aafn=opts.ant_array, deltas=opts.deltas)
        antGains = None # Setup variable so that the gain table isn't re-read for every file if used
        if lofarStation.name=='KAIRAXST': dataFmt='KAIRA'

    ####################
    ## Read Visibilities
    ####################
    visFiles = args # filenames to image
    for vid,visFn in enumerate(visFiles):
        print 'Using %s (%i/%i)'%(visFn, vid+1, len(visFiles))
        fDict = SWHT.fileio.parse(visFn, fmt=dataFmt)

        # Pull out the visibility data in a (u,v,w) format
        if fDict['fmt']=='acc': # LOFAR station all subbands ACC file visibilities
            decomp = True
            fDict['rcu'] = opts.rcumode # add the RCU mode to the meta data of an ACC file
            fDict['sb'] = sbs # select subbands to use
            fDict['int'] = opts.int_time # set integration length (usually 1 second)

            vis, uvw, freqs, obsInfo = SWHT.fileio.readACC(visFn, fDict, lofarStation, sbs, calTable=opts.calfile)
            [obsLat, obsLong, LSTangle] = obsInfo

            # add visibilities to previously processed files
            visComb = np.concatenate((visComb, vis), axis=1)
            uvwComb = np.concatenate((uvwComb, uvw), axis=0)

        elif fDict['fmt']=='xst': # LOFAR XST format visibilities
            decomp = True
            if opts.override: # Override XST filename metadata
                fDict['rcu'] = opts.rcumode
                fDict['sb'] = sbs
                fDict['int'] = opts.int_time
            else:
                sbs = fDict['sb']

            vis, uvw, freqs, obsInfo = SWHT.fileio.readXST(visFn, fDict, lofarStation, sbs, calTable=opts.calfile)
            [obsLat, obsLong, LSTangle] = obsInfo

            # add visibilities to previously processed files
            visComb = np.concatenate((visComb, vis), axis=1)
            uvwComb = np.concatenate((uvwComb, uvw), axis=0)

        elif fDict['fmt']=='KAIRA': # KAIRA LOFAR XST format visibilities
            decomp = True
            if opts.override: # Override XST filename metadata
                fDict['rcu'] = opts.rcumode
                fDict['sb'] = sbs
                fDict['int'] = opts.int_time
            else:
                sbs = fDict['sb']

            vis, uvw, freqs, obsInfo = SWHT.fileio.readKAIRAXST(visFn, fDict, lofarStation, sbs, times=opts.times)
            [obsLat, obsLong, LSTangle] = obsInfo

            # add visibilities to previously processed files
            visComb = np.concatenate((visComb, vis), axis=1)
            uvwComb = np.concatenate((uvwComb, uvw), axis=0)

        elif fDict['fmt']=='ms': # MS-based visibilities
            decomp = True

            fDict['sb'] = sbs

            vis, uvw, freqs, obsInfo = SWHT.fileio.readMS(visFn, sbs, column=opts.column)
            [obsLat, obsLong, LSTangle] = obsInfo

            # add visibilities to previously processed files
            visComb = np.concatenate((visComb, vis), axis=1)
            uvwComb = np.concatenate((uvwComb, uvw), axis=0)

        else:
            print 'ERROR: unknown data format, exiting'
            exit()

    print 'AUTO-CORRELATIONS:', opts.autos
    if not opts.autos: # remove auto-correlations
        autoIdx = np.argwhere(uvwComb[:,0]**2. + uvwComb[:,1]**2. + uvwComb[:,2]**2. == 0.)
        visComb[:,autoIdx] = 0.

    # prepare for Fourier transform
    print 'Performing Fourier Transform'
    pixels = opts.pixels
    px = [pixels,pixels]
    fov = opts.fov * (np.pi/180.) #Field of View in radians
    res = fov / px[0] #pixel resolution
    print 'Resolution(deg):', res*180./np.pi

    # convert uvw to wavelength, reduce dimensions of uvw and vis arrays for Fourier Transform
    visComb = np.reshape(visComb, (4, visComb.shape[1]*visComb.shape[2])) # 4 pols x (number of samples * number of subbands)
    uvwComb[:,0] *= freqs/cc
    uvwComb[:,1] *= freqs/cc
    uvwComb[:,2] *= freqs/cc
    uvwComb = np.swapaxes(uvwComb, 1, 2)
    uvwComb = np.reshape(uvwComb, ( uvwComb.shape[0]*uvwComb.shape[1], 3))

    # uvw coordinates returned by the SWHT functions are in reference to a (HA=LST, dec=+90) phase centre, need to rotate uvw coordinates to (HA=0, dec=obsLat)
    ha = float(LSTangle) # Hour Angle is set to LST, we need to reverse the rotation
    haRotMat = np.array([   [    np.sin(ha), np.cos(ha), 0.],
                            [-1.*np.cos(ha), np.sin(ha), 0.],
                            [0.,             0.,         1.]]) #rotate about z-axis
    dec = obsLat # observatory latitude
    decRotMat = np.array([  [1.,              0.,          0.],
                            [0.,     np.sin(dec), np.cos(dec)],
                            [0., -1.*np.cos(dec), np.sin(dec)]]) #rotate about x-axis
    ha0 = 0.
    haRotMat0 = np.array([  [    np.sin(ha0), np.cos(ha0), 0.],
                            [-1.*np.cos(ha0), np.sin(ha0), 0.],
                            [0.,              0.,         1.]]) #rotate about z-axis
    rotMatrix = np.dot( decRotMat, np.dot( haRotMat0, np.dot(haRotMat.T, np.identity(3).T))) # reverse the hour angle rotation in the fileio.read function, the dec rotation is the identity, so that does not need to be reversed
    uvwComb = np.dot(rotMatrix, uvwComb.T).T

    if opts.uvplot: # display the projected UV coverage
        fig, ax = SWHT.display.dispVis2D(uvwComb)
        plt.show()

    # Plot the PSF
    if opts.psf:
        visComb = np.ones_like(visComb)

    # perform DFT or FFT
    if opts.dft:
        print 'DFT'
        xxIm = SWHT.ft.dftImage(visComb[0], uvwComb, px, res, mask=False, rescale=False, stokes=False)
        xyIm = SWHT.ft.dftImage(visComb[1], uvwComb, px, res, mask=False, rescale=False, stokes=False)
        yxIm = SWHT.ft.dftImage(visComb[2], uvwComb, px, res, mask=False, rescale=False, stokes=False)
        yyIm = SWHT.ft.dftImage(visComb[3], uvwComb, px, res, mask=False, rescale=False, stokes=False)
    else:
        print 'FFT'
        conv = opts.conv #rotate about z-axisv
        xxIm = SWHT.ft.fftImage(visComb[0], uvwComb, px, res, mask=False, wgt=opts.weighting, conv=conv)
        xyIm = SWHT.ft.fftImage(visComb[1], uvwComb, px, res, mask=False, wgt=opts.weighting, conv=conv)
        yxIm = SWHT.ft.fftImage(visComb[2], uvwComb, px, res, mask=False, wgt=opts.weighting, conv=conv)
        yyIm = SWHT.ft.fftImage(visComb[3], uvwComb, px, res, mask=False, wgt=opts.weighting, conv=conv)

    #save complex image to pickle file
    if opts.pkl is None: outPklFn = 'tempImage.pkl'
    else: outPklFn = opts.pkl
    if opts.dft: fttype = 'dft'
    else: fttype = opts.conv
    print 'Writing image to file %s ...'%outPklFn,
    SWHT.fileio.writeImgPkl(outPklFn, np.array([xxIm,xyIm,yxIm,yyIm]), fDict, res=res, fttype=fttype, imtype='complex')
    print 'done'
    
    #display stokes plots
    if not opts.nodisplay or not (opts.savefig is None):
        fig, ax = SWHT.display.disp2DStokes(xxIm, xyIm, yxIm, yyIm)
    if not (opts.savefig is None): plt.savefig(opts.savefig)
    if not opts.nodisplay: plt.show()

