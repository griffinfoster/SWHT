#!/usr/bin/env python
"""
Simple script to convert the output of GSM (http://space.mit.edu/~angelica/gsm/index.html) into a HEALPIX FITS file
"""

#try to us pandas to read in the text file as it has a faster parser than numpy
useNumpy = False
try:
    import pandas as pd
except ImportError:
    useNumpy = True

import numpy as np
import healpy as hp
import sys,os
import time

if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] GSM_OUTPUT_FILE')
    o.set_description(__doc__)
    o.add_option('-C', '--celestial', dest='celestial', action='store_true',
        help = 'Output the HEALPIX map in Celestial coordinates instead of Galactic')
    opts, args = o.parse_args(sys.argv[1:])

    print 'Using pandas:', not useNumpy
    for fid,fn in enumerate(args):
        outfn = fn + '.hpx'
        print 'Converting: %s to %s (%i/%i)'%(fn, outfn, fid+1, len(args))

        start_time = time.time()
        if useNumpy:
            m = np.genfromtxt(fn) #slightly faster than: m = np.loadtxt(fn)
        else:
            m = pd.read_csv(fn, sep=' ', header=None, engine='c')[3].values
        
        nside = hp.npix2nside(m.shape[0])
        npix = m.shape[0]
        print '\t(%f s) NSIDE: %i PIXELS: %i'%(time.time() - start_time, nside, npix)

        if opts.celestial:
            print '\tConverting to Celestial Coordinates'
            #the output of gsm is in galactic coordinates, but it will be more useful in celestial coordinates for interferometry
            r = hp.Rotator(coord=['C','G'])  # Transforms celestial coordinates to galactic
            celeTheta, celePhi = hp.pix2ang(nside, np.arange(npix)) #get celestial coordinates from pixel numbers
            galTheta, galPhi = r(celeTheta, celePhi) #convert celestial coordinates to galactic
            galPix = hp.ang2pix(nside, galTheta, galPhi) #get pixel numbers from galactic coordinates
            #Now, this next bit is just a little bit dodgy as there is not a direct 1-to-1 conversion between pixels, some galactic pixels are repeated and some are dropped
            mCele = m[galPix]
            hp.write_map(outfn, mCele, coord='C')
        else:
            hp.write_map(outfn, m, coord='G')

    print 'Finished'

