#!/usr/bin/env python
"""
Plot a HEALPIX FITS map
"""

import sys,os
import numpy as np
#import matplotlib
#matplotlib.use('GTKAgg')
import matplotlib.pyplot as plt
import healpy as hp

if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] FITS_MAP')
    o.set_description(__doc__)
    o.add_option('-m', '--mode', dest='mode', default='I',
        help='Scale mode: H(histogram), I(linear), L(log), S(sqrt); default: I')
    o.add_option('-p', '--proj', dest='proj', default='CG',
        help='Projection of map, E(cliptic), G(alactic), C(equatorial), default: CG')
    o.add_option('-v', '--visu', dest='visu', default='moll',
        help='Visualization projection: Mollweide(moll), Cartesian(cart), Orthographic(ortho) default: moll')
    o.add_option('-t', '--title', dest='title', default='',
        help='Map title, default: Sky Map')
    o.add_option('-s', '--savefig', dest='savefig', default=None,
        help='Save map to file')
    o.add_option('--max', dest='max', default=None, type='float',
        help='Max flux value, default: None')
    o.add_option('--min', dest='min', default=None, type='float',
        help='Min flux value, default: None')
    o.add_option('--dec_min', dest='dec_min', default=None, type='float',
        help='Min declination to plot, in degrees, only works propertly if the input map is in Celestial coordinates, default: None')
    o.add_option('--dec_max', dest='dec_max', default=None, type='float',
        help='Max declination to plot, in degrees, only works propertly if the input map is in Celestial coordinates, default: None')
    o.add_option('-w', '--weight', dest='weight_map', action='store_true',
        help='Plot the sample weighting as a map instead of the data')
    o.add_option('--cmap', dest='cmap', default='jet',
        help='Set a matplotlib colormap(jet,hot,gist_heat,gist_earth,PRGn,RdYlBu,spectral), default:jet')
    o.add_option('--cbar', dest='cbar', default=True, action='store_false',
        help='Disable colorbar')
    o.add_option('--fill', dest='fill', default=False, action='store_true',
        help='Fill UNSEEN/NaN values with zeros')
    o.add_option('--rot', dest='rot', default=None, type='str',
        help='Rotate the map such that the center is at (lon,lat)')
    o.add_option('--size', dest='xsize', default=800, type='int',
        help='Figure size, pick around 5000 for a large figure, default: 800')
    opts, args = o.parse_args(sys.argv[1:])

    cbar = opts.cbar
    unit = 'Flux (Jy/Beam)'
    if opts.weight_map: unit = 'Samples'
    if opts.rot:
        lon,lat = opts.rot.split(',')
        rot = (float(lon),float(lat),0.)
    else: rot = None

    xsize = opts.xsize

    m = None
    w = None
    coordSys = 'C' #default to celestial coordinates
    for fid,fn in enumerate(args):
        print 'Opening: %s (%i/%i)'%(fn, fid+1, len(args))
        hpMap = hp.read_map(fn, field=None, h=True)
        hdr = hpMap[-1]
        if len(hpMap)==2: #no weight map
            if m is None: m = hpMap[0]
            else: m += hpMap[0]
        elif len(hpMap)==3: #weight map
            if m is None: m,w,hdr = hpMap
            else:
                m += hpMap[0]
                w += hpMap[1]

    #set map projection
    for item in hdr: #get coordinate system from the header
        if item[0].startswith('COORDSYS'): coordSys = item[1]
    if len(opts.proj) == 1:
        coord_proj = coordSys + opts.proj.upper()
    else: #override the projection
        coord_proj = opts.proj.upper()
    print 'Using coordinate projection: %s --> %s'%(coord_proj[0], coord_proj[1])

    if w is not None: m /= w #divide by the pixel weights
    
    print 'Map :: min=%f :: max=%f'%(np.nanmin(m), np.nanmax(m))

    ma_m = np.ma.masked_array(m, np.isnan(m))

    if opts.weight_map: m = w

    if opts.mode.lower()=='s': m = np.sqrt(m - np.min(np.nan_to_num(m)))

    #mask out declination regions, only works properly when the input map is in Celestial coordinates
    nside = hp.pixelfunc.get_nside(m)
    if opts.dec_min is None: dec_min = 180.
    else:
        dec_min = 90. - opts.dec_min
        theta_min = (dec_min / 180.) * np.pi
        ring_min = hp.pixelfunc.ang2pix(nside, theta_min, 0.)
        m[ring_min:] = hp.pixelfunc.UNSEEN
    if opts.dec_max is None: dec_max = 0.
    else:
        dec_max = 90. - opts.dec_max
        theta_max = (dec_max / 180.) * np.pi
        ring_max = hp.pixelfunc.ang2pix(nside, theta_max, 0.)
        m[:ring_max] = hp.pixelfunc.UNSEEN

    #replace nan values with zero
    if opts.fill: m = np.nan_to_num(m)

    fig = plt.figure(1, figsize=(12,8))

    if opts.mode.lower()=='h':
        norm = 'hist'
    elif opts.mode.lower()=='i':
        norm = None
    elif opts.mode.lower()=='l':
        norm = 'log'
    elif opts.mode.lower()=='s':
        norm = None

    if opts.visu.lower().startswith('moll'): #Mollweide projection
        hp.mollview(m, coord=coord_proj, norm=norm, max=opts.max, min=opts.min, cmap=plt.get_cmap(opts.cmap), cbar=cbar, fig=1, unit=unit, title=opts.title, xsize=xsize, rot=rot)
    elif opts.visu.lower().startswith('cart'): #Cartesian
        hp.cartview(m, coord=coord_proj, norm=norm, max=opts.max, min=opts.min, cmap=plt.get_cmap(opts.cmap), cbar=cbar, fig=1, unit=unit, title=opts.title, xsize=xsize, rot=rot)
    elif opts.visu.lower().startswith('ortho'): #Orthographic
        hp.orthview(m, coord=coord_proj, norm=norm, max=opts.max, min=opts.min, cmap=plt.get_cmap(opts.cmap), cbar=cbar, fig=1, unit=unit, title=opts.title, xsize=xsize, rot=rot)

    #set grid lines
    hp.graticule(dpar=5, dmer=5)

    if not(opts.savefig is None):
        plt.savefig(opts.savefig)
    else: plt.show()
    
