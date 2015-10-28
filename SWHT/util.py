"""
Utility functions
"""

import numpy as np
import datetime

def sph2cart(theta, phi, r=None):
    """Convert spherical coordinates to 3D cartesian
    theta, phi, and r must be the same size and shape, if no r is provided then unit sphere coordinates are assumed (r=1)
    theta: colatitude/elevation angle, 0(north pole) =< theta =< pi (south pole)
    phi: azimuthial angle, 0 <= phi <= 2pi
    r: radius, 0 =< r < inf
    returns X, Y, Z arrays of the same shape as theta, phi, r
    see: http://mathworld.wolfram.com/SphericalCoordinates.html
    """
    if r is None: r = np.ones_like(theta) #if no r, assume unit sphere radius

    #elevation is pi/2 - theta
    #azimuth is ranged (-pi, pi]
    X = np.cos((np.pi/2.)-theta) * np.cos(phi-np.pi)
    Y = np.cos((np.pi/2.)-theta) * np.sin(phi-np.pi)
    Z = np.sin((np.pi/2.)-theta)

    return X, Y, Z

def cart2sph(X, Y, Z):
    """Convert 3D cartesian coordinates to spherical coordinates
    X, Y, Z: arrays of the same shape and size
    returns r: radius, 0 =< r < inf
            phi: azimuthial angle, 0 <= phi <= 2pi
            theta: colatitude/elevation angle, 0(north pole) =< theta =< pi (south pole)
    see: http://mathworld.wolfram.com/SphericalCoordinates.html
    """
    r = np.sqrt(X**2. + Y**2. + Z**2.)
    phi = np.arctan2(Y, X) + np.pi #convert azimuth (-pi, pi] to phi (0, 2pi]
    theta = np.pi/2. - np.arctan2(Z, np.sqrt(X**2. + Y**2.)) #convert elevation [pi/2, -pi/2] to theta [0, pi]

    return r, phi, theta

def vectorize(mat):
    """Convert upper-left triangle of mat to rank 1 vector
    """
    idx = np.triu_indices(mat.shape[0])
    return mat[idx]

def vectorize3D(mat):
    """Convert upper-left triangle of a 3D array to rank 1 vector, assumes the first axis is 
    """
    idx = np.triu_indices(mat.shape[1])
    return mat[:, idx[0], idx[1]]

def convert_arg_range(arg):
    """Split apart command-line lists/ranges into a list of numbers."""
    if arg is None: return None

    arg = arg.split(',')
    outList = []
    for aa in arg:
        rr = map(int, aa.split('_'))
        if len(rr)==1: outList.append(rr[0])
        elif len(rr)==2:
            outList.extend(range(rr[0], rr[1]+1))
    return outList

def meanTimeDelta(l):
    """Return the mean of a list of datetime.timedelta objects"""
    tDelta = datetime.timedelta(0)
    for td in l:
        tDelta += td
    return tDelta/len(l)

#Functions taken from healpy/sphtfunc.py
class Alm(object):
    """This class provides some static methods for alm index computation.

    Methods
    -------
    getlm
    getidx
    getsize
    getlmax
    """
    def __init__(self):
        pass

    @staticmethod
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

    @staticmethod
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

    @staticmethod
    def getsize(lmax,mmax = None):
        """Returns the size of the array needed to store alm up to *lmax* and *mmax*

        Parameters
        ----------
        lmax : int
          The maximum l, defines the alm layout
        mmax : int, optional
          The maximum m, defines the alm layout. Default: lmax.

        Returns
        -------
        size : int
          The size of the array needed to store alm up to lmax, mmax.
        """
        if mmax is None or mmax < 0 or mmax > lmax:
            mmax = lmax
        return mmax * (2 * lmax + 1 - mmax) // 2 + lmax + 1

    @staticmethod
    def getlmax(s, mmax = None):
        """Returns the lmax corresponding to a given array size.
        
        Parameters
        ----------
        s : int
          Size of the array
        mmax : None or int, optional
          The maximum m, defines the alm layout. Default: lmax.

        Returns
        -------
        lmax : int
          The maximum l of the array, or -1 if it is not a valid size.
        """
        if mmax is not None and mmax >= 0:
            x = (2 * s + mmax ** 2 - mmax - 2) / (2 * mmax + 2)
        else:
            x = (-3 + np.sqrt(1 + 8 * s)) / 2
        if x != np.floor(x):
            return -1
        else:
            return int(x)

def almVec2array(vec, lmax):
    """Convert the vector output of healpy.map2alm into a 2-D array of the same format as used in the SWHT
    healpy.map2alm returns coefficients for 0=<l<=lmax and 0<=m<=l
    vec: output of healpy.map2alm
    lmax: maximum l number used in healpy.map2alm"""
    lmaxp1 = lmax + 1 #account for the 0 mode
    coeffs = np.zeros((lmaxp1, 2*lmaxp1-1), dtype='complex')

    for l in np.arange(lmaxp1):
        for m in np.arange(l+1):
            #These calls use the healpy Alm() calls which are reproduced in util.py
            #coeffs[l,l-m] = ((-1.)**m) * np.conj(vec[hp.Alm.getidx(lmax,l,m)]) #since the map is real, the a_l,-m = (-1)**m * a_l,m.conj
            #coeffs[l,l+m] = vec[hp.Alm.getidx(lmax,l,m)]
            coeffs[l,l-m] = ((-1.)**m) * np.conj(vec[Alm.getidx(lmax,l,m)]) #since the map is real, the a_l,-m = (-1)**m * a_l,m.conj
            coeffs[l,l+m] = vec[Alm.getidx(lmax,l,m)]

    return coeffs

def array2almVec(arr):
    """Convert a 2-D array of coefficients used in the SWHT into a vector that is the same as that of healpy.map2alm such that healpy.alm2map can be called with the output vector
    healpy.map2alm returns coefficients for 0=<l<=lmax and 0<=m<=l
    arr: 2-D array of SWHT coefficients [lmax + 1, 2*lmax + 1]
    """
    lmax = arr.shape[0] - 1
    ncoeffs = Alm.getsize(lmax)
    vec = np.zeros(ncoeffs, dtype='complex')

    for l in np.arange(lmax + 1):
        for m in np.arange(l + 1):
            vec[Alm.getidx(lmax,l,m)] = arr[l,l+m]

    return vec

import numpy as np

if __name__ == '__main__':
    print 'Running test cases'

    [theta, phi] = np.meshgrid(np.linspace(0, np.pi, num=128, endpoint=False), np.linspace(0, 2.*np.pi, num=128, endpoint=False))
    X, Y, Z = sph2cart(theta, phi)
    r0, phi0, theta0 = cart2sph(X, Y, Z)

    print np.allclose(theta, theta0)
    print np.allclose(phi, phi0)

    aa = np.arange(5*5).reshape(5,5)
    print aa
    print vectorize(aa)

    from matplotlib import pyplot as plt
    plt.subplot(221)
    plt.imshow(theta, interpolation='nearest')
    plt.colorbar()
    plt.subplot(222)
    plt.imshow(phi, interpolation='nearest')
    plt.colorbar()
    plt.subplot(223)
    plt.imshow(theta0, interpolation='nearest')
    plt.colorbar()
    plt.subplot(224)
    plt.imshow(phi0, interpolation='nearest')
    plt.colorbar()
    plt.show()

    print 'Made it through without any errors.'

