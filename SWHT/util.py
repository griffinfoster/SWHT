"""
Utility functions
"""

import numpy as np

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

