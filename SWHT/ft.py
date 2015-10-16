"""
Fourier Transform Functions

functions phsCenterSrc, eq2top_m, get_baseline, gen_uvw, xyz2uvw are taken from AIPY (https://github.com/AaronParsons/aipy), used to compute (U,V,W) from ITRF (X,Y,Z)
"""

import numpy as np
import ephem
import sys,os
import struct
import time

def phsCenterSrc(obs, t):
    """return an ephem FixedBody source based on the time offset from the obs"""
    src = ephem.FixedBody()
    t0 = obs.date
    obs.date = t
    src._ra = obs.sidereal_time()
    src._dec = obs.lat
    obs.date = t0
    return src

def eq2top_m(ha, dec):
    """Return the 3x3 matrix converting equatorial coordinates to topocentric
    at the given hour angle (ha) and declination (dec)."""
    sin_H, cos_H = np.sin(ha), np.cos(ha)
    sin_d, cos_d = np.sin(dec), np.cos(dec)
    zero = np.zeros_like(ha)
    map =  np.array([[    sin_H    ,       cos_H  ,       zero  ],
                        [ -sin_d*cos_H,   sin_d*sin_H,      cos_d  ],
                        [  cos_d*cos_H,  -cos_d*sin_H,      sin_d  ]])
    if len(map.shape) == 3: map = map.transpose([2, 0, 1])
    return map

def get_baseline(i, j, src, obs):
    """Return the baseline corresponding to i,j""" 
    bl = j - i
    try:
        if src.alt < 0:
            raise PointingError('Phase center below horizon')
        m=src.map
    except(AttributeError):
        ra,dec = src._ra,src._dec
        #opposite HA since the we want to move the source at zenith away to phase to the original zenith source
        m = eq2top_m(ra-obs.sidereal_time(), dec)
        #normal HA
        #m = eq2top_m(obs.sidereal_time() - ra, dec)
    return np.dot(m, bl).transpose()

def gen_uvw(i, j, src, obs, f):
    """Compute uvw coordinates of baseline relative to provided FixedBody"""
    x,y,z = get_baseline(i,j,src,obs)
    afreqs = np.reshape(f, (1,f.size))
    afreqs = afreqs/ephem.c #1/wavelength
    if len(x.shape) == 0: return np.array([x*afreqs, y*afreqs, z*afreqs]).T
    x.shape += (1,); y.shape += (1,); z.shape += (1,)
    return np.array([np.dot(x,afreqs), np.dot(y,afreqs), np.dot(z,afreqs)]).T

def xyz2uvw(xyz, src, obs, f):
    """Return an array of UVW values"""
    uvw = np.zeros((f.shape[0], xyz.shape[0], xyz.shape[0], 3))
    for i in range(xyz.shape[0]):
        for j in range(xyz.shape[0]):
            if i==j: continue
            uvw[:, i, j, :] = gen_uvw(xyz[i], xyz[j], src, obs, f)[:,0,:]
    return uvw

def dft2(d, l, m, u, v, psf=False):
    """compute the 2d DFT for position (m,l) based on (d,uvw)"""
    if psf: return np.sum(np.exp(2.*np.pi*1j*((u*l) + (v*m))))/u.size
    else: return np.sum(d*np.exp(2.*np.pi*1j*((u*l) + (v*m))))/u.size

def dftImage(d, uvw, px, res, mask=False, rescale=False, stokes=False):
    """return a DFT image
    d: complex visibilities [F, Q] F frequency subbands, Q samples
    uvw: visibility sampling in units of wavelengths [Q, 3]
    px: [int, int], number of pixels in image
    res: float, resolution of central pixel in radians
    rescale: account for missing np.sqrt(1-l^2-m^2) in flat-field approximation
    """
    if stokes: im = np.zeros((px[0], px[1], 4),dtype=complex)
    else: im = np.zeros((px[0], px[1]), dtype=complex)
    maskIm = np.zeros((px[0], px[1]), dtype=bool)

    mid_m = int(px[0]/2.) #middle pixel number in m direction
    mid_l = int(px[1]/2.) #middle pixel number in l direction

    u = np.array(uvw[:,0])
    v = np.array(uvw[:,1])
    w = np.array(uvw[:,2])

    #fov = [px[0]*res*(180./np.pi), px[1]*res*(180./np.pi)] #Field of View in degrees
    #set (l,m) range based on the number of pixels and resolution
    lrange = np.linspace(-1.*px[0]*res/2., px[0]*res/2., num=px[0], endpoint=True)/(np.pi/2.) #m range (-1,1)
    mrange = np.linspace(-1.*px[1]*res/2., px[1]*res/2., num=px[1], endpoint=True)/(np.pi/2.) #l range (-1,1)

    start_time = time.time()
    for mid,m in enumerate(mrange):
        for lid,l in enumerate(lrange):
            #rescale to account for missing np.sqrt(1-l^2-m^2) in flat-field approximation
            if rescale: scale = np.sqrt(1.-(l**2.)-(m**2.))
            else: scale = 1.
            
            if stokes:
                im[lid,mid,0] = dft2(d[0], l, m, u, v) * scale
                im[lid,mid,1] = dft2(d[1], l, m, u, v) * scale
                im[lid,mid,2] = dft2(d[2], l, m, u, v) * scale
                im[lid,mid,3] = dft2(d[3], l, m, u, v) * scale
            else: im[lid,mid] = dft2(d, m, l,  u, v) * scale
            if mask: #mask out region beyond field of view
                rad = (m**2 + l**2)**.5
                if rad > 1.: maskIm[lid,mid] = True
    print time.time() - start_time

    im = np.flipud(np.fliplr(im)) #make top-left corner (0,0) the south-east point
    maskIm = np.flipud(np.fliplr(maskIm))
    
    if mask: return im, maskIm
    else: return im

def fftImage(d, uvw, px, res, mask=False, conv='fast', wgt='natural'):
    """Grid visibilities and perform an FFT to return an image
    d: complex visibilities
    uvw: visibility sampling in units of wavelengths
    px: [int, int], number of pixels in image
    res: float, resolution of central pixel in radians
    """
    start_time = time.time()

    im = np.zeros((px[0], px[1]), dtype=complex)
    maskIm = np.zeros((px[0], px[1]), dtype=bool)

    mid_m = int(px[0]/2.) #middle pixel number in m direction
    mid_l = int(px[1]/2.) #middle pixel number in l direction

    u = np.array(uvw[:,0])
    v = np.array(uvw[:,1])
    w = np.array(uvw[:,2])

    gridVis = np.zeros((px[0], px[1]), dtype=complex) #array for the gridded visibilities
    gridWgt = np.ones((px[0], px[1]), dtype=float) #array for the grid weights

    #u,v grid spacing based on the number of pixels and resolution of the desired image
    deltau = (np.pi/2.) * 1./(px[0]*res)
    deltav = (np.pi/2.) * 1./(px[1]*res)

    if conv.startswith('fast'):
        for did,dd in enumerate(d):
            #simple, rectangular convolution function (nearest neighbor convolution)
            uu = int(u[did]/deltau)
            vv = int(v[did]/deltav)
            gridVis[(uu+(px[0]/2))%px[0], (vv+(px[1]/2))%px[1]] += dd
    else:
        gridUV = np.mgrid[-.5*px[0]*deltau:.5*px[0]*deltau:deltau, -.5*px[1]*deltav:.5*px[1]*deltav:deltav]

        #choose a convolution function to grid with
        if conv.startswith('rect'):
            convFunc = convRect(deltau, deltav)
            truncDist = deltau/2. #truncate box car to within a single pixel
        if conv.startswith('gauss'):
            convFunc = convGauss(deltau/2., deltav/2.) #half-power point at 1/2 a UV pixel distance
            truncDist = deltau*4. #truncate the convolution function to a 4 pixel radius
        if conv.startswith('prolate'):
            convFunc = convProlate(deltau, deltav)
            truncDist = deltau #only values of sqrt(u**2 + v**2) < deltau are valid

        #Grid visibilities
        for uid in range(px[0]):
            for vid in range(px[1]):
                ucentre,vcentre = gridUV[:,uid,vid]
                udiff = u-ucentre #distance between the ungridded u positon and the gridded u position
                vdiff = v-vcentre #distance between the ungridded v positon and the gridded v position
                idx = np.argwhere(np.sqrt((udiff**2.)+(vdiff**2.)) < truncDist) #convolution function should be truncated at a reasonable kernel size
                if idx.size > 0:
                    gridWgt[uid,vid] = np.sum(convFunc(udiff[idx],vdiff[idx]))
                    gridVis[uid,vid] = np.sum(convFunc(udiff[idx],vdiff[idx]) * d[idx])

    if wgt.startswith('uni'): gridVis /= gridWgt #uniform weighting, default is natural weighting
            
    gridVis = np.fft.fftshift(gridVis) #(0,0) position is in the middle, need to shift it to a corner
    im = np.fft.fftshift(np.fft.fft2(gridVis)) #shift (0,0) back to the middle
    im = np.rot90(np.fliplr(im)) #make top-left corner (0,0) the south-east point
    
    print time.time() - start_time
    if mask: return im, maskIm
    else: return im

def convGauss(ures, vres, alpha=1.):
    """Return a Gaussian convolution function
    ures,vres: distance from centre to half power point in uv distance
    alpha: scaling factor"""
    return lambda uu,vv: ((1./(alpha*np.sqrt(ures*vres*np.pi)))**2.)*np.exp(-1.*(uu/(alpha*ures))**2.)*np.exp(-1.*(vv/(alpha*vres))**2.)
    
def convRect(ures, vres):
    """Return a boxcar/rectangle convolution function"""
    return lambda uu,vv: np.ones_like(uu)

def convProlate(ures, vres, aa=1., cc=1.):
    """Return a prolate spheroid function which returns the function z(uu, vv) = sqrt( cc**2. * ( 1. - (((uu*ures)**2. + (vv*vres)**2.)/aa**2.))), c > a for a prolate function"""
    return lambda uu,vv: np.sqrt( cc**2. * (1. - (((uu/ures)**2. + (vv/vres)**2.)/aa**2.)))

if __name__ == '__main__':
    print 'Running test cases'

    #TODO: add tests

    print 'Made it through without any errors.'


