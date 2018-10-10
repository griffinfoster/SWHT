#!/usr/bin/env python
"""
Perform a Spherical Wave Harmonic Transform on LOFAR ACC/XST data or widefield MS data (e.g. PAPER) to form a complex or Stokes dirty image dirty image
"""

import numpy as np
import matplotlib.pyplot as plt
import datetime
import sys,os
import SWHT

try:
    import healpy as hp
    healpyEnabled = True
except ImportError:
    healpyEnabled = False

#try:
#    import casacore.tables as tbls
#except ImportError:
#    print 'Warning: could not import casacore.tables, will not be able to read measurement sets'

#import scipy.constants
#cc = scipy.constants.c
cc = 299792458.0 #speed of light, m/s

if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] ACC/XST/MS/PKL FILE')
    o.set_description(__doc__)
    o.add_option('--station', dest='station', default=None,
        help = 'LOFAR ONLY: station name, e.g. SE607, if this is used then the ant_field and ant_array options are not required, default: None')
    o.add_option('-F', '--ant_field', dest='ant_field', default=None,
        help = 'LOFAR ONLY: AntennaField.conf file for the LOFAR station of the ACC files, default: None')
    o.add_option('-A', '--ant_array', dest='ant_array', default=None,
        help = 'LOFAR ONLY(NOT REQUIRED): AntennaArray.conf file for the LOFAR station geographical coordinates, default: None')
    o.add_option('-D', '--deltas', dest='deltas', default=None,
        help = 'LOFAR ONLY: iHBADeltas.conf file, only required for HBA imaging, default: None')
    o.add_option('-r', '--rcumode', dest='rcumode', default=3, type='int',
        help = 'LOFAR ONLY: Station RCU Mode, usually 3,5,6,7, for XST it will override filename metadata default: 3(LBA High)')
    o.add_option('-s', '--subband', dest='subband', default='0',
        help = 'Select which subband(s) to image, for ACC and MS it will select, for multiple subbands use X,Y,Z and for range use X_Y notation, for XST it will override filename metadata, default:0')
    o.add_option('-p', '--pixels', dest='pixels', default=64, type='int',
        help = 'Width of 2D image in pixels, or number of steps in 3D image, or NSIDE in HEALPix image, default: 64')
    o.add_option('-C', '--cal', dest='calfile', default=None,
        help = 'LOFAR ONLY: Apply a calibration soultion file to the data.')
    o.add_option('-S', '--save', dest='savefig', default=None,
        help = 'Save the figure using this name, type is determined by extension')
    o.add_option('--nodisplay', dest='nodisplay', action='store_true',
        help = 'Do not display the generated image')
    o.add_option('--of', dest='of', default=None,
        help = 'Save complex images in a numpy array in a pickle file or HEALPix map using this name (include .pkl or .hpx extention), default: tempImage.pkl')
    o.add_option('-i', '--int', dest='int_time', default=1., type='float',
        help = 'LOFAR ONLY: Integration time, used for accurate zenith pointing, for XST it will override filename metadata, default: 1 second')
    o.add_option('-c', '--column', dest='column', default='CORRECTED_DATA', type='str',
        help = 'MS ONLY: select which data column to image, default: CORRECTED_DATA')
    o.add_option('--override', dest='override', action='store_true',
        help = 'LOFAR XST ONLY: override filename metadata for RCU, integration length, and subband')
    o.add_option('--autos', dest='autos', action='store_true',
        help = 'Include the auto-correlations in the image')
    o.add_option('--fov', dest='fov', default=180., type='float',
        help = '2D IMAGING MODE ONLY: Field of View in degrees, default: 180 (all-sky)')
    o.add_option('-l', '--lmax', dest='lmax', default=32, type='int',
        help = 'Maximum l spherical harmonic quantal number, rule-of-thumb: lmax ~ (pi/longest baseline resolution), default: 32')
    o.add_option('--lmin', dest='lmin', default=0, type='int',
        help = 'Minimum l spherical harmonic quantal number, usually left as 0, default: 0')
    o.add_option('--ocoeffs', dest='ocoeffs', default=None,
        help = 'Save output image coefficients to a pickle file using this name (include .pkl extention), default: tempCoeffs.pkl')
    o.add_option('-I', '--image', dest='imageMode', default='2D',
        help='Imaging mode: 2D (hemisphere flattened), 3D, healpix, coeff (coefficients) default: 2D')
    o.add_option('--uvwplot', dest='uvwplot', action='store_true',
        help='Display a 3D UVW coverage/sampling plot')
    o.add_option('--psf', dest='psf', action='store_true',
        help='Plot the PSF instead of the image')
    o.add_option('-t', '--times', dest='times', default='0',
        help = 'KAIRA ONLY: Select which integration(s) to image, can use a[seconds] to average, d[step size] to decimate, of a specific range of integrations similar to the subband selection option, default:0 (select the first integration of the file)')
    o.add_option('--pol', dest='polMode', default='I',
        help='Polarization selection: I, Q, U, V, XX, YY, XY, YX, default: I')
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
        if lofarStation.name=='KAIRA': dataFmt='KAIRA'

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

        elif fDict['fmt']=='pkl':
            print 'Loading Image Coefficients file:', visFn
            coeffDict = SWHT.fileio.readCoeffPkl(visFn)
            imgCoeffs = coeffDict['coeffs']
            LSTangle = coeffDict['lst']
            obsLong = coeffDict['phs'][0]
            obsLat = coeffDict['phs'][1]
            decomp = False

        else:
            print 'ERROR: unknown data format, exiting'
            exit()

    if opts.uvwplot: # display the total UVW coverage
        fig, ax = SWHT.display.dispVis3D(uvwComb)
        plt.show()

    ####################
    ## Decompose the input visibilities into spherical harmonics visibility coefficeints
    ####################
    if decomp:
        # compute the ideal l_max given the average solid angle angular resolution of an l-mode is Omega ~ 4pi / 2l steradian, and if the PSF is circular theta ~ pi / l radians
        blLen = np.sqrt(uvwComb[:,0,:]**2. + uvwComb[:,1,:]**2. + uvwComb[:,2,:]**2.) # compute the baseline lengths (in meters)
        maxBl = np.max(blLen) # maximum baseline length (in meters)
        meanWl = cc / np.mean(freqs) # mean observing wavelength
        maxRes = 1.22 * meanWl / maxBl
        print 'MAXIMUM RES: %f (radians) %f (deg)'%(maxRes, maxRes * (180. / np.pi))
        idealLmax = int(np.pi / maxRes)
        print 'SUGGESTED L_MAX: %i, %i (oversample 3), %i (oversample 5)'%(idealLmax, idealLmax*3, idealLmax*5)

        if opts.psf:
           visComb = np.ones_like(visComb) 

        print 'AUTO-CORRELATIONS:', opts.autos
        if not opts.autos: # remove auto-correlations
            autoIdx = np.argwhere(uvwComb[:,0]**2. + uvwComb[:,1]**2. + uvwComb[:,2]**2. == 0.)
            visComb[:,autoIdx] = 0.

        # prepare for SWHT
        print 'Performing Spherical Wave Harmonic Transform'
        print 'LMAX:', opts.lmax

        polMode = opts.polMode.upper()
        print 'Polarization Mode:', polMode
        if polMode=='I': polVisComb = visComb[0] + visComb[3]
        elif polMode=='Q': polVisComb = visComb[0] - visComb[3]
        elif polMode=='U': polVisComb = visComb[1] + visComb[2]
        elif polMode=='V': polVisComb = 1j * np.conj(visComb[1] - visComb[2]) #flip imaginary and real
        elif polMode=='XX': polVisComb = visComb[0]
        elif polMode=='XY': polVisComb = visComb[1]
        elif polMode=='YX': polVisComb = visComb[2]
        elif polMode=='YY': polVisComb = visComb[3]

        imgCoeffs = SWHT.swht.swhtImageCoeffs(polVisComb, uvwComb, freqs, lmax=opts.lmax, lmin=opts.lmin) # perform SWHT

        # save image coefficients to file
        if opts.ocoeffs is None: outCoeffPklFn = 'tempCoeffs.pkl'
        else: outCoeffPklFn = opts.pkl
        SWHT.fileio.writeCoeffPkl(outCoeffPklFn, imgCoeffs, [float(obsLong), float(obsLat)], float(LSTangle))

    ####################
    ## Imaging
    ####################
    if opts.of is None:
        if opts.imageMode.startswith('heal'): 
            outFn = 'tempImage.hpx'
            if os.path.exists(outFn):
                os.remove(outFn)
        else: outFn = 'tempImage.pkl'
    else: outFn = opts.of

    #TODO: not doing the correct projection
    if opts.imageMode.startswith('2'): # Make a 2D hemispheric image
        fov = opts.fov * (np.pi/180.) # Field of View in radians
        px = [opts.pixels, opts.pixels]
        res = fov/px[0] # pixel resolution
        print 'Generating 2D Hemisphere Image of size (%i, %i)'%(px[0], px[1])
        print 'Resolution(deg):', res*180./np.pi
        #img = SWHT.swht.make2Dimage(imgCoeffs, res, px, phs=[0., 0.]) #TODO: 0 because the positions have already been rotated to the zenith RA of the first snapshot, if multiple snaphsots this needs to be reconsidered
        img = SWHT.swht.make2Dimage(imgCoeffs, res, px, phs=[float(LSTangle), obsLat]) #TODO: 0 because the positions have already been rotated to the zenith RA of the first snapshot, if multiple snaphsots this needs to be reconsidered
        #img = SWHT.swht.make2Dimage(imgCoeffs, res, px, phs=[0., float(obsLat)]) #TODO: 0 because the positions have already been rotated to the zenith RA of the first snapshot, if multiple snaphsots this needs to be reconsidered
        fig, ax = SWHT.display.disp2D(img, dmode='abs', cmap='jet')

        # save complex image to pickle file
        print 'Writing image to file %s ...'%outFn,
        SWHT.fileio.writeSWHTImgPkl(outFn, img, fDict, mode='2D')
        print 'done'

    elif opts.imageMode.startswith('3'): # Make a 3D equal stepped image
        print 'Generating 3D Image with %i steps in theta and %i steps in phi'%(opts.pixels, opts.pixels)
        img, phi, theta = SWHT.swht.make3Dimage(imgCoeffs, dim=[opts.pixels, opts.pixels])
        fig, ax = SWHT.display.disp3D(img, phi, theta, dmode='abs', cmap='jet')

        # save complex image to pickle file
        print 'Writing image to file %s ...'%outFn,
        SWHT.fileio.writeSWHTImgPkl(outFn, [img, phi, theta], fDict, mode='3D')
        print 'done'

    elif opts.imageMode.startswith('heal'): # plot healpix and save healpix file using the opts.pkl name
        print 'Generating HEALPix Image with %i NSIDE'%(opts.pixels)
        # use the healpy.alm2map function as it is much faster, there is a ~1% difference between the 2 functions, this is probably due to the inner workings of healpy
        #m = SWHT.swht.makeHEALPix(imgCoeffs, nside=opts.pixels) # TODO: a rotation issue
        m = hp.alm2map(SWHT.util.array2almVec(imgCoeffs), opts.pixels) # TODO: a rotation issue

        # save complex image to HEALPix file
        print 'Writing image to file %s ...'%outFn,
        hp.write_map(outFn, m.real, coord='C') # only writing the real component, this should be fine, maybe missing some details, but you know, the sky should be real.
        print 'done'
    
    elif opts.imageMode.startswith('coeff'): # plot the complex coefficients
        fig, ax = SWHT.display.dispCoeffs(imgCoeffs, zeroDC=True, vis=False)

    if not (opts.savefig is None): plt.savefig(opts.savefig)
    if not opts.nodisplay:
        if opts.imageMode.startswith('heal'): hp.mollview(m.real, coord='CG')
        plt.show()

