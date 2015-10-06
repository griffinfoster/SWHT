"""
Functions for producing SWHT-based dirty images
Based on T. Carozzi MATLAB code
"""

import numpy as np
import scipy.special
import math
import time

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
    vis *= 2. #Treat the conjugate baslines as doubling the nonconjugate visibilities

    kr = np.dot(r, k.T) #compute the radii in wavelengths for each visibility sample
    vislm = np.zeros((lmax+1, 2*lmax+1, vis.shape[1]), dtype='complex')
    for l in np.arange(lmax+1): #increase lmax by 1 to account for starting from 0
        print l
        for m in np.arange(-1*l, l+1):
            #Compute visibility spherical harmonic coefficients according to SWHT, i.e. multiply visibility by spherical wave harmonics for each L&M and sum over all baselines.
            #Note that each non-zero baseline is effectively summed twice in the preceding formula. (In the MNRAS letter image the NZ baselines were only weighted once, i.e. their conjugate baselines were not summed.)
            #print l,m
            spharmlm = np.repeat(np.conj(spharm(l, m, theta, phi)), vis.shape[1], axis=1) #spherical harmonics only needed to be computed once for all baselines, independent of observing frequency, TODO: spharm is by far the slowest call
            vislm[l, l+m] = (((2.*(k**2.))/np.pi) * np.sum(vis * sphBj(l, kr) * spharmlm, axis=0)[np.newaxis].T).flatten() #sum visibilites of same obs frequency
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
    k = freqs/cc #obs freq/c

    #convert u,v,w to r,phi,theta
    r = np.sqrt(uvw[:,0]**2. + uvw[:,1]**2. + uvw[:,2]**2.)[np.newaxis].T
    phi = np.arctan2(uvw[:,1], uvw[:,0])[np.newaxis].T
    theta = (np.pi/2.) - np.arctan2(uvw[:,2], np.sqrt(uvw[:,0]**2. + uvw[:,1]**2.))[np.newaxis].T

    #compute the SWHT visibility coefficients
    vislm = computeVislm(lmax, k, r, theta, phi, vis)
    #compute the SWHT brightness coefficients
    blm = computeblm(vislm)

    print time.time() - start_time

    return blm

if __name__ == "__main__":
    print 'Running test cases'

    import matplotlib.pyplot as plt

    jl0 = sphBj(0, np.linspace(0, 50, num=256))
    jl1 = sphBj(1, np.linspace(0, 50, num=256))
    jl2 = sphBj(2, np.linspace(0, 50, num=256))
    #plt.plot(jl0)
    #plt.plot(jl1)
    #plt.plot(jl2)

    theta, phi = np.meshgrid(np.linspace(0, np.pi, 100), np.linspace(0,2.*np.pi, 100))
    #theta, phi = np.meshgrid(np.linspace(-1.*np.pi/2., np.pi/2., 100), np.linspace(0,2.*np.pi, 100))

    l = 2
    m = -1
    #Yp = scipy.special.sph_harm(n=l, m=m, theta=phi, phi=theta)
    Y = spharm(l=l, m=m, theta=theta, phi=phi)
    #Yp = scipy.special.sph_harm(n=l, m=m, theta=phi, phi=theta)
    #Yp = 0.5 * np.sqrt(3./np.pi) * np.cos(theta) #1,0
    #Yp = -0.5 * np.sqrt(3./(2.*np.pi)) * np.sin(theta) * np.exp(phi * 1j) #1,1
    #Yp = 0.5 * np.sqrt(3./(2.*np.pi)) * np.sin(theta) * np.exp(phi * -1j) #1,-1
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
