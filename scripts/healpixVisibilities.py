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

#Function taken from healpy/sphtfunc.py
def getlm(lmax,i=None):
    """Get the l and m from index and lmax.
    
    Parameters
    ----------
    lmax : int
      The maximum l defining the alm layout
    i : int or None
      The index for which to compute the l and m.
      If None, the function return l and m for i=0..Alm.getsize(lmax)
    """
    if i is None:
        i=np.arange(Alm.getsize(lmax))
    m=(np.ceil(((2*lmax+1)-np.sqrt((2*lmax+1)**2-8*(i-lmax)))/2)).astype(int)
    l = i-m*(2*lmax+1-m)//2
    return (l,m)

#Function taken from healpy/sphtfunc.py
def getidx(lmax,l,m):
    """Returns index corresponding to (l,m) in an array describing alm up to lmax.
    
    Parameters
    ----------
    lmax : int
      The maximum l, defines the alm layout
    l : int
      The l for which to get the index
    m : int
      The m for which to get the index
    Returns
    -------
    idx : int
      The index corresponding to (l,m)
    """
    return m*(2*lmax+1-m)//2+l

def almVec2array(vec, lmax):
    """Convert the vector output of healpy.map2alm into a 2-D array of the same format as used in the SWHT
    healpy.map2alm returns coefficients for 0=<l<=lmax and 0<=m<=l
    vec: output of healpy.map2alm
    lmax: maximum l number used in healpy.map2alm"""
    lmaxp1 = lmax + 1 #account for the 0 mode
    coeffs = np.zeros((lmaxp1, 2*lmaxp1+1), dtype='complex')

    for l in np.arange(lmaxp1):
        for m in np.arange(l+1):
            coeffs[l,l+m] = vec[getidx(lmax,l,m)]
            coeffs[l,l-m] = ((-1.)**m) * np.conj(vec[getidx(lmax,l,m)]) #since the map is real, the a_l,-m = (-1)**m * a_l,m.conj
            #these lines make the array symetric
            #coeffs[l,l+m + (lmax-l+1)] = vec[getidx(lmax,l,m)]
            #coeffs[l,l-m + (lmax-l+1)] = ((-1.)**m) * np.conj(vec[getidx(lmax,l,m)]) #since the map is real, the a_l,-m = (-1)**m * a_l,m.conj

    return coeffs

if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] -i FITS_MAP MS/ACC/XST FILES')
    o.set_description(__doc__)
    o.add_option('-i', '--imap', dest='imap', default=None,
        help='Input HEALPIX map or spherical harmonics coefficient file')
    o.add_option('-l', '--lmax', dest='lmax', default=32, type='int',
        help = 'Maximum l spherical harmonic quantal number, default: 32')
    o.add_option('--dec_min', dest='dec_min', default=None, type='float',
        help='Min declination to plot, in degrees, default: None')
    o.add_option('--dec_max', dest='dec_max', default=None, type='float',
        help='Max declination to plot, in degrees, default: None')
    opts, args = o.parse_args(sys.argv[1:])

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

    #mask out declination regions
    nside = hp.pixelfunc.get_nside(m)
    if opts.dec_min is None: dec_min = 180.
    else: dec_min = 90. - opts.dec_min
    if opts.dec_max is None: dec_max = 0.
    else: dec_max = 90. - opts.dec_max
    theta_min = (dec_min / 180.) * np.pi
    theta_max = (dec_max / 180.) * np.pi
    ring_min = hp.pixelfunc.ang2pix(nside, theta_min, 0.)
    ring_max = hp.pixelfunc.ang2pix(nside, theta_max, 0.)
    m[ring_min:] = hp.pixelfunc.UNSEEN
    m[:ring_max] = hp.pixelfunc.UNSEEN

    #hp.mollview(m)
    #plt.show()

    #Convert HEALPIX map into Alm spherical harmonics coefficients
    alms = hp.sphtfunc.map2alm(m, lmax=opts.lmax, mmax=opts.lmax)
    blm = almVec2array(alms, opts.lmax)
    #plt.imshow(10.*np.log10(np.abs(blm)), interpolation='nearest')
    #plt.imshow(blm.imag, interpolation='nearest')
    #plt.colorbar()
    #plt.show()

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

