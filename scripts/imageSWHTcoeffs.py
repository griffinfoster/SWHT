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
        img = np.fliplr(img)
        plt.imshow(np.abs(img), interpolation='nearest')
        plt.colorbar()

        #save complex image to pickle file
        print 'Writing image to file %s ...'%outFn,
        SWHT.fileio.writeSWHTImgPkl(outFn, img, fDict, mode='2D')
        print 'done'

    elif opts.imageMode.startswith('3'): #Make a 3D equal stepped image
        print 'Generating 3D Image with %i steps in theta and %i steps in phi'%(opts.pixels, opts.pixels)
        img, phi, theta = SWHT.swht.make3Dimage(iImgCoeffs, dim=[opts.pixels, opts.pixels])
        img = np.abs(img)
        #X = np.cos(theta-(np.pi/2.)) * np.cos(phi)
        #Y = np.cos(theta-(np.pi/2.)) * np.sin(phi)
        #Z = np.sin(theta-(np.pi/2.))
        X, Y, Z = SWHT.util.sph2cart(theta, phi)

        #http://stackoverflow.com/questions/22175533/what-is-the-equivalent-of-matlabs-surfx-y-z-c-in-matplotlib
        from mpl_toolkits.mplot3d import Axes3D
        from matplotlib import cm
        from matplotlib.colors import Normalize
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        imin = img.min()
        imax = img.max()
        scalarMap = cm.ScalarMappable(norm=Normalize(vmin=imin, vmax=imax), cmap=cm.jet)
        #scalarMap = cm.ScalarMappable(norm=Normalize(vmin=imin, vmax=imax), cmap=cm.gist_earth_r)
        scalarMap.set_array(img)
        C = scalarMap.to_rgba(img)
        surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors=C, antialiased=True)
        fig.colorbar(scalarMap)

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
        hp.write_map(outFn, np.abs(m), coord='C')
        print 'done'
    
    elif opts.imageMode.startswith('coeff'): #plot the complex coefficients
        if opts.viscoeffs: #convert brightness coefficients to visibility coefficients
            iImgCoeffs = SWHT.swht.computeblm(iImgCoeffs, reverse=True)
        #iImgCoeffs[0,0] = 0 #zero out DC offset component

        plt.subplot(231)
        plt.title('Real Components')
        plt.imshow(iImgCoeffs.real, interpolation='nearest')
        plt.colorbar()

        plt.subplot(232)
        plt.title('Imaginary Components')
        plt.imshow(iImgCoeffs.real, interpolation='nearest')
        plt.imshow(iImgCoeffs.imag, interpolation='nearest')
        plt.colorbar()

        plt.subplot(234)
        plt.title('Amplitude (dB)')
        plt.imshow(iImgCoeffs.real, interpolation='nearest')
        plt.imshow(10.*np.log10(np.abs(iImgCoeffs)), interpolation='nearest')
        plt.colorbar()

        plt.subplot(235)
        plt.title('Phase')
        plt.imshow(iImgCoeffs.real, interpolation='nearest')
        plt.imshow(np.angle(iImgCoeffs), interpolation='nearest')
        plt.colorbar()

        plt.subplot(233)
        coeffsFlat = []
        mms = []
        lls = []
        for ll in np.arange(iImgCoeffs.shape[0]):
            for mm in np.arange(-1*ll, ll+1):
                mms.append(mm)
                lls.append(ll)
                coeffsFlat.append(iImgCoeffs[ll,ll+mm])
        coeffsFlat = np.array(coeffsFlat)
        plt.ylabel('Amplitude (dB)')
        plt.xlabel('l')
        plt.plot(lls, 10.*np.log10(np.abs(coeffsFlat)), '.')

        plt.subplot(236)
        plt.ylabel('Amplitude (dB)')
        plt.xlabel('m')
        plt.plot(mms, 10.*np.log10(np.abs(coeffsFlat)), '.')

    if not (opts.savefig is None): plt.savefig(opts.savefig)
    if not opts.nodisplay:
        if opts.imageMode.startswith('heal'): hp.mollview(np.abs(m))
        plt.show()
    
