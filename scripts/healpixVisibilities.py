#!/usr/bin/env python
"""
Convert a HEALPIX map into visibilities based on input Measurement Sets of LOFAR ACC, XST files
"""

import sys,os
import numpy as np
import healpy as hp
import SWHT

#TODO: LOFAR option: station id
#TODO: LOFAR: generate new XST based on input subband option
#TODO: MS option: column to write to
#TODO: include beam
#TODO: using dec_min/max creates sharp discontinuities which results in a poor decomposition
import matplotlib.pyplot as plt

if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] -i FITS_MAP MS/ACC/XST FILES')
    o.set_description(__doc__)
    o.add_option('-i', '--imap', dest='imap', default=None,
        help='REQUIRED: Input HEALPIX map or spherical harmonics coefficient file')
    o.add_option('-l', '--lmax', dest='lmax', default=32, type='int',
        help = 'Maximum l spherical harmonic quantal number, default: 32')
    o.add_option('--dec_min', dest='dec_min', default=None, type='float',
        help='Min declination to plot, in degrees, default: None')
    o.add_option('--dec_max', dest='dec_max', default=None, type='float',
        help='Max declination to plot, in degrees, default: None')
    opts, args = o.parse_args(sys.argv[1:])

    if opts.imap is None:
        print "ERROR: no input HEALPIX map set with -i/--imap option"
        exit()

    #Get HEALPIX map
    m = None
    w = None
    print 'Opening:', opts.imap
    hpMap = hp.read_map(opts.imap, field=None, h=True)
    if len(hpMap)==2: #no weight map
        m,hdr = hpMap
    elif len(hpMap)==3: #weight map
        m,w,hdr = hpMap

    if w is not None: m /= w #divide by the pixel weights
    
    print 'Map :: min=%f :: max=%f'%(np.nanmin(m), np.nanmax(m))

    #hp.mollview(m)
    #plt.show()

    #Convert HEALPIX map into Alm spherical harmonics coefficients
    alms = hp.sphtfunc.map2alm(m, lmax=opts.lmax, mmax=opts.lmax)
    blm = SWHT.util.almVec2array(alms, opts.lmax)
    #print alms.shape
    #print hp.Alm.getlm(opts.lmax)
    alms0 = SWHT.util.array2almVec(blm)
    print alms-alms0

    exit()

    #plt.imshow(10.*np.log10(np.abs(blm)), interpolation='nearest')
    #plt.imshow(blm.imag, interpolation='nearest')
    #plt.colorbar()
    #plt.show()
    m = SWHT.swht.makeHEALPix(blm, nside=32)
    m0 = hp.alm2map(alms, 32)
    #print np.max(np.abs(m-m0))
    #hp.mollview(np.abs(m-m0))
    hp.mollview(np.abs(m))
    plt.show()

    #m0 = hp.alm2map(alms, 64)
    #hp.mollview(m0)
    #plt.show()

    #Convert brightness coefficients into visibility coefficients
    vislm = SWHT.swht.computeblm(blm, reverse=True)
    #plt.imshow(10.*np.log10(np.abs(vislm)), interpolation='nearest')
    #plt.imshow(blm.imag, interpolation='nearest')
    #plt.colorbar()
    #plt.show()

    #Generate visibilities with SWHT.swht.computeVisSamples()
    #if LOFAR station: get station information
    #for each visibility file:
    #   if MS:
    #       open MS, get XYZ positions, convert to r,theta,phi
    #       compute visibilities
    #       save visibilities to column
    #   if LOFAR:
    #       get file metadata
    #       get XYZ positions, convert to r,theta,phi
    #       compute visibilities
    #       generate a new binary file and save visibilities

