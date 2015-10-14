"""
Functions for producing SWHT-based dirty images
Based on T. Carozzi MATLAB code
"""

#TODO: speedup spharm: multi-core, write function instead of scipy.special, use Ylm
#TODO: tests

import numpy as np
import scipy.special
import math
import time
import sys

import Ylm
import util

cc = 299792458. #speed of light, m/s

def sphBj(l, r):
    """spherical Bessel function of first kind
    l: int, order
    r: positive float array, radius"""
    jl = np.sqrt(np.pi/(2.*r)) * scipy.special.jv( l+0.5, r)

    rorgInd = np.argwhere(r == 0.)
    if len(rorgInd) > 0:
        if l==0.: jl[rorgInd]=1.
        else: jl[rorgInd]=0.

    return jl

def spharm(l, m, theta, phi):
    """Spherical harmonic function with l,m polar and azimuthal quantal numbers
    theta: float array, polar/elevation angle, range (-pi/2, pi/2)
    phi: float array, azimith angle, range (0, 2pi)
    l: positive int
    m: int, -l <= m <= l
    this is just a wrapper scipy.special.sph_harm, we use the convention that phi is the azimuthal angle [0, 2pi], theta is the colatitude/elevation angle [0,pi] where 0 is the norht pole and pi is the south pole
    """
    #for some reason, if m<0 the call is much slower, but Y_l,-m is the complex conjugate of -1.*Y_l,m
    if m<0:
        return -1.*np.conjugate(scipy.special.sph_harm(n=l, m=np.abs(m), theta=phi, phi=theta)) #scipy uses non-standard notation
    else:
        return scipy.special.sph_harm(n=l, m=m, theta=phi, phi=theta) #scipy uses non-standard notation

def computeVislm(lmax, k, r, theta, phi, vis):
    """Compute the spherical wave harmonics visibility coefficients, Eq. 16 of Carozzi 2015
    lmax: positive int, maximum spherical harmonic l number
    k: [N, 1] float array, wave number, observing frequencies/c (1/meters)
    r, theta, phi: [Q, 1] float arrays of visibility positions transformed from (u,v,w) positions, r (meters)
    vis: [Q, N] complex array, observed visibilities

    returns: [lmax+1, 2*lmax+1, nfreq] array of coefficients, only partially filled, see for loops in this function
    """
    #vis *= 2. #Treat the conjugate baslines as doubling the nonconjugate visibilities

    kr = np.dot(r, k.T) #compute the radii in wavelengths for each visibility sample
    vislm = np.zeros((lmax+1, 2*lmax+1, vis.shape[1]), dtype='complex')

    print 'L:',
    for l in np.arange(lmax+1): #increase lmax by 1 to account for starting from 0
        print l,
        sys.stdout.flush()
        for m in np.arange(-1*l, l+1):
            #Compute visibility spherical harmonic coefficients according to SWHT, i.e. multiply visibility by spherical wave harmonics for each L&M and sum over all baselines.
            #Note that each non-zero baseline is effectively summed twice in the preceding formula. (In the MNRAS letter image the NZ baselines were only weighted once, i.e. their conjugate baselines were not summed.)
            spharmlm = np.repeat(np.conj(Ylm.Ylm(l,m,phi,theta)), vis.shape[1], axis=1) #spherical harmonics only needed to be computed once for all baselines, independent of observing frequency, TODO: spharm is by far the slowest call
            vislm[l, l+m] = (((2.*(k**2.))/np.pi) * np.sum(vis * sphBj(l, kr) * spharmlm, axis=0)[np.newaxis].T).flatten() #sum visibilites of same obs frequency
    print 'done'

    return vislm

def computeblm(vislm, reverse=False):
    """Compute the spherical wave harmonics brightness coefficients from the spherical wave harmonics visibility coefficients, Eq. 11 of Carozzi 2015
    vislm: complex array, spherical wave harmonics visibility coefficients computed from computeVislm()
    reverse: if true, convert vislm from input blm
    """
    blm = np.zeros_like(vislm)
    for l in np.arange(blm.shape[0]):
        for m in np.arange(-1*l, l+1):
            if reverse: blm[l, l+m] = vislm[l, l+m] * (4.*np.pi*((-1.*1j)**float(l)))
            else: blm[l, l+m] = vislm[l, l+m] / (4.*np.pi*((-1.*1j)**float(l)))
    return blm

def swhtImageCoeffs(vis, uvw, freqs, lmax):
    """Generate brightness coefficients based converting visibilities with the SWHT
    vis: complex array [Q, F], Q observed visibilities at F frequencies, can be size [Q] if only using 1 frequency
    uvw: float array [Q, 3], meters
    freqs: float array [F, 1] or [F], observing frequencies in Hz
    lmax: positive int, maximum spherical harmonic l number
    """
    start_time = time.time()

    if vis.ndim==1: vis = vis[np.newaxis].T

    if freqs.ndim==1: freqs = freqs[np.newaxis].T
    k = 2. * np.pi * freqs/cc #obs freq/c

    #convert u,v,w to r,phi,theta
    r, phi, theta = util.cart2sph(uvw[:,0], uvw[:,1], uvw[:,2])
    #make arrays 2D
    r = r[np.newaxis].T
    phi = np.pi - phi[np.newaxis].T #make range -pi to pi
    theta = theta[np.newaxis].T
    #r = np.sqrt(uvw[:,0]**2. + uvw[:,1]**2. + uvw[:,2]**2.)[np.newaxis].T
    #phi = np.arctan2(uvw[:,1], uvw[:,0])[np.newaxis].T
    #theta = (np.pi/2.) - np.arctan2(uvw[:,2], np.sqrt(uvw[:,0]**2. + uvw[:,1]**2.))[np.newaxis].T #make range -pi/2 to pi/2

    #compute the SWHT visibility coefficients
    vislm = computeVislm(lmax, k, r, theta, phi, vis)
    #compute the SWHT brightness coefficients
    blm = computeblm(vislm)

    print time.time() - start_time

    return blm

#TODO: make2Dimage: FoV
def make2Dimage(coeffs, dim=[64, 64]):
    """Make a flat image of a single hemisphere from SWHT image coefficients
    coeffs: SWHT brightness coefficients
    dim: [int, int], number of pixels, note these are equivalent to the l,m coordinates in FT imaging
    """
    start_time = time.time()

    #start from a regular Cartesian grid
    xx,yy = np.meshgrid(np.linspace(-1., 1., num=dim[0]), np.linspace(-1., 1., num=dim[1]))
    img = np.zeros(xx.shape, dtype='complex')
    
    #convert to polar positions
    r = np.sqrt(xx**2. + yy**2.)
    phi = np.arctan2(yy, xx)
    #zero out undefiined regions of the image where r>0
    #numpy is super cunty and makes it fucking hard to do simple things, so the next few lines are bullshit that should only take 2 lines
    idx = np.argwhere(r.flatten()>1)
    rflat = r.flatten()
    phiflat = phi.flatten()
    rflat[idx] = 0.
    phiflat[idx] = 0.
    r = np.reshape(rflat, r.shape)
    phi = np.reshape(phiflat, phi.shape)

    #convert to unit sphere coordinates
    thetap = np.pi/2 - np.arccos(r) #zenith is at pi in spherical coordinates
    phip = phi + np.pi #azimuth range [0, 2pi]

    lmax = coeffs.shape[0]
    print 'L:',
    for l in np.arange(lmax):
        print l,
        sys.stdout.flush()
        for m in np.arange(-1*l, l+1):
            img += coeffs[l, l+m] * Ylm.Ylm(l, m, phip, thetap) #TODO: a slow call
    print 'done'

    print time.time() - start_time

    return img

#TODO: make3Dimage: masking
def make3Dimage(coeffs, dim=[64, 64]):
    """Make a 3D sphere from SWHT image coefficients
    coeffs: SWHT brightness coefficients
    dim: [int, int] number of steps in theta and phi
    """
    #equal-spaced sample of theta and phi, not ideal as equal pixel (i.e. HEALPIX)
    [theta, phi] = np.meshgrid(np.linspace(0, np.pi, num=dim[0], endpoint=True), np.linspace(0, 2.*np.pi, num=dim[1], endpoint=True))
    img = np.zeros(theta.shape, dtype='complex')

    lmax = coeffs.shape[0]
    print 'L:',
    for l in np.arange(lmax):
        print l,
        sys.stdout.flush()
        for m in np.arange(-1*l, l+1):
            img += coeffs[l, l+m] * Ylm.Ylm(l, m, phi, theta) #TODO: a slow call
    print 'done'

    return img, phi, theta

if __name__ == "__main__":
    print 'Running test cases'

    import matplotlib.pyplot as plt

    #theta, phi = np.meshgrid(np.linspace(0, np.pi, 10), np.linspace(0,2.*np.pi, 10))
    #Y = spharm(l=1, m=-1, theta=theta, phi=phi)
    #print Y
    #exit()

    jl0 = sphBj(0, np.linspace(0, 50, num=256))
    jl1 = sphBj(1, np.linspace(0, 50, num=256))
    jl2 = sphBj(2, np.linspace(0, 50, num=256))
    #plt.plot(jl0)
    #plt.plot(jl1)
    #plt.plot(jl2)

    theta, phi = np.meshgrid(np.linspace(0, np.pi, 100), np.linspace(0,2.*np.pi, 100))
    #theta, phi = np.meshgrid(np.linspace(-1.*np.pi/2., np.pi/2., 100), np.linspace(0,2.*np.pi, 100))

    l = 1
    m = -1
    #Yp = scipy.special.sph_harm(n=l, m=m, theta=phi, phi=theta)
    Y = spharm(l=l, m=m, theta=theta, phi=phi)
    #Yp = scipy.special.sph_harm(n=l, m=m, theta=phi, phi=theta)
    #Yp = 0.5 * np.sqrt(3./np.pi) * np.cos(theta) #1,0
    #Yp = -0.5 * np.sqrt(3./(2.*np.pi)) * np.sin(theta) * np.exp(phi * 1j) #1,1
    Yp = 0.5 * np.sqrt(3./(2.*np.pi)) * np.sin(theta) * np.exp(phi * -1j) #1,-1
    print np.allclose(Y, Yp, atol=1e-08)

    #plt.subplot(231)
    #plt.imshow(Y.real, interpolation='nearest', extent=[0., np.pi, 0, 2.*np.pi])
    #plt.xlabel('theta')
    #plt.ylabel('phi')
    #plt.colorbar()

    #plt.subplot(234)
    #plt.imshow(Y.imag, interpolation='nearest', extent=[0., np.pi, 0, 2.*np.pi])
    #plt.xlabel('theta')
    #plt.ylabel('phi')
    #plt.colorbar()

    #plt.subplot(232)
    #plt.imshow(Yp.real, interpolation='nearest', extent=[0., np.pi, 0, 2.*np.pi])
    #plt.xlabel('theta')
    #plt.ylabel('phi')
    #plt.colorbar()

    #plt.subplot(235)
    #plt.imshow(Yp.imag, interpolation='nearest', extent=[0., np.pi, 0, 2.*np.pi])
    #plt.xlabel('theta')
    #plt.ylabel('phi')
    #plt.colorbar()

    #plt.subplot(233)
    #plt.imshow((Y-Yp).real, interpolation='nearest', extent=[0., np.pi, 0, 2.*np.pi])
    #plt.xlabel('theta')
    #plt.ylabel('phi')
    #plt.colorbar()

    #plt.subplot(236)
    #plt.imshow((Y-Yp).imag, interpolation='nearest', extent=[0., np.pi, 0, 2.*np.pi])
    #plt.xlabel('theta')
    #plt.ylabel('phi')
    #plt.colorbar()

    theta = np.random.rand(100) * np.pi
    theta = theta[np.newaxis].T
    phi = np.random.rand(100) * np.pi*2.
    phi = phi[np.newaxis].T
    r = np.random.rand(100) * 60. #meters
    r = r[np.newaxis].T
    k = np.array([100., 110., 120.])*1e6/299792458.0 #obs freq/c
    k = k[np.newaxis].T
    vis = np.random.rand(3*100)*2. -1. + 1j*(np.random.rand(3*100)*2.-1.)
    vis = np.reshape(vis, (100,3))
    lmax = 5
    vislm = computeVislm(lmax, k, r, theta, phi, vis)
    print vislm.shape

    blm = computeblm(vislm)
    print blm.shape

    freqs = np.array([100., 110., 120.])*1e6
    uvw = ((np.random.rand(300) * 60.) - 30.).reshape(100, 3)
    blm = swhtImageCoeffs(vis, uvw, freqs, lmax)
    print blm.shape

    #plt.show()

    print 'Made it through without any errors.'
