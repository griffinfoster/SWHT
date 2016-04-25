#!/usr/bin/env python
"""
Produce an image from a set of SWHT brightness coefficients
"""

import numpy as np
from matplotlib import pyplot as plt
import healpy as hp
import sys,os
import SWHT

#import scipy.constants
#cc = scipy.constants.c
cc = 299792458.0 #speed of light, m/s

if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] PKL FILE')
    o.set_description(__doc__)
    o.add_option('-p', '--pixels', dest='pixels', default=64, type='int',
        help = 'Width of 2D image in pixels, or number of steps in 3D image, or NSIDE in HEALPix image, default: 64')
    o.add_option('-S', '--save', dest='savefig', default=None,
        help = 'Save the figure using this name, type is determined by extension')
    o.add_option('--nodisplay', dest='nodisplay', action='store_true',
        help = 'Do not display the generated image')
    o.add_option('--of', dest='of', default=None,
        help = 'Save complex images in a numpy array in a pickle file or HEALPix map using this name (include .pkl or .hpx extention), default: tempImage.pkl')
    o.add_option('--fov', dest='fov', default=180., type='float',
        help = '2D IMAGING MODE ONLY: Field of View in degrees, default: 180 (all-sky)')
    o.add_option('-I', '--image', dest='imageMode', default='2D',
        help='Imaging mode: 2D (hemisphere flattened), 3D, healpix, coeff (coefficients) default: 2D')
    o.add_option('--vis', dest='viscoeffs', action='store_true',
        help='If plotting coefficients, convert them to visibility coefficients')
    opts, args = o.parse_args(sys.argv[1:])

    #get filenames to image
    #TODO: can add together image coefficients, just need to consider the meta data setup
    coeffFiles = args
    for cid,coeffFn in enumerate(coeffFiles):
        print 'Using %s (%i/%i)'%(coeffFn, cid+1, len(coeffFiles))
        fDict = SWHT.fileio.parse(coeffFn)

        if fDict['fmt']=='pkl':
            print 'Loading Image Coefficients file:', coeffFn
            coeffDict = SWHT.fileio.readCoeffPkl(coeffFn)
            iImgCoeffs = coeffDict['coeffs']
            LSTangle = coeffDict['lst']
            obsLong = coeffDict['phs'][0]
            obsLat = coeffDict['phs'][1]
            decomp = False
        else:
            print 'ERROR: unknown data format, exiting'
            exit()

    #Imaging
    if opts.of is None:
        if opts.imageMode.startswith('heal'): outFn = 'tempImage.hpx'
        else: outFn = 'tempImage.pkl'
    else: outFn = opts.of

    if opts.imageMode.startswith('2'): #Make a 2D hemispheric image
        fov = opts.fov * (np.pi/180.) #Field of View in radians
        px = [opts.pixels, opts.pixels]
        res = fov/px[0] #pixel resolution
        print 'Generating 2D Hemisphere Image of size (%i, %i)'%(px[0], px[1])
        print 'Resolution(deg):', res*180./np.pi
        img = SWHT.swht.make2Dimage(iImgCoeffs, res, px, phs=[0., float(obsLat)]) #0 because the positions have already been rotated to the zenith RA of the first snapshot, if multiple snaphsots this needs to be reconsidered
        fig, ax = SWHT.display.disp2D(img, dmode='abs', cmap='jet')

        #save complex image to pickle file
        print 'Writing image to file %s ...'%outFn,
        SWHT.fileio.writeSWHTImgPkl(outFn, img, fDict, mode='2D')
        print 'done'

    elif opts.imageMode.startswith('3'): #Make a 3D equal stepped image
        print 'Generating 3D Image with %i steps in theta and %i steps in phi'%(opts.pixels, opts.pixels)
        img, phi, theta = SWHT.swht.make3Dimage(iImgCoeffs, dim=[opts.pixels, opts.pixels])
        fig, ax = SWHT.display.disp3D(img, phi, theta, dmode='abs', cmap='jet')

        #save complex image to pickle file
        print 'Writing image to file %s ...'%outFn,
        SWHT.fileio.writeSWHTImgPkl(outFn, [img, phi, theta], fDict, mode='3D')
        print 'done'

    elif opts.imageMode.startswith('heal'): #plot healpix and save healpix file using the opts.pkl name
        print 'Generating HEALPix Image with %i NSIDE'%(opts.pixels)
        #use the healpy.alm2map function as it is much faster, there is a ~1% difference between the 2 functions, this is probably due to the inner workings of healpy
        #m = SWHT.swht.makeHEALPix(iImgCoeffs, nside=opts.pixels)
        m = hp.alm2map(SWHT.util.array2almVec(iImgCoeffs), opts.pixels)

        #save complex image to HEALPix file
        print 'Writing image to file %s ...'%outFn,
        hp.write_map(outFn, np.abs(m), coord='C') #TODO: should this be abs or real?
        print 'done'
    
    elif opts.imageMode.startswith('coeff'): #plot the complex coefficients
        fig, ax = SWHT.display.dispCoeffs(iImgCoeffs, zeroDC=True, vis=opts.viscoeffs)

    if not (opts.savefig is None): plt.savefig(opts.savefig)
    if not opts.nodisplay:
        if opts.imageMode.startswith('heal'): hp.mollview(np.abs(m))
        plt.show()

