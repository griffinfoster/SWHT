"""Utility functions
"""

#TODO: make sure sph2cart and cart2sph are consistant with each other

def sph2cart(theta, phi, r=None):
    """Convert spherical coordinates to 3D cartesian
    theta, phi, and r must be the same size and shape, if no r is provided then unit sphere coordinates are assumed (r=1)
    theta: colatitude/elevation angle, 0(north pole) =< theta =< pi (south pole)
    phi: azimuthial angle, 0 <= phi <= 2pi
    r: radius, 0 =< r < inf
    returns X, Y, Z arrays of the same shape as theta, phi, r
    """
    if r is None: r = np.ones_like(theta) #if no r, assume unit sphere radius

    X = r * np.sin(theta) * np.cos(phi)
    Y = r * np.sin(theta) * np.sin(phi)
    Z = r * np.cos(theta)

    return X, Y, Z

def cart2sph(X, Y, Z):
    """Convert 3D cartesian coordinates to spherical coordinates
    X, Y, Z: arrays of the same shape and size
    returns r: radius, 0 =< r < inf
            phi: azimuthial angle, 0 <= phi <= 2pi
            theta: colatitude/elevation angle, 0(north pole) =< theta =< pi (south pole)
    """
    r = np.sqrt(X**2. + Y**2. + Z**2.)
    phi = np.arctan2(Y, X)
    theta = np.arctan2(Z, np.sqrt(X**2. + Y**2.))

    return r, phi, theta

import numpy as np

if __name__ == '__main__':
    print 'Running test cases'

    [theta, phi] = np.meshgrid(np.linspace(0, np.pi, num=128, endpoint=False), np.linspace(0, 2.*np.pi, num=128, endpoint=False))
    X, Y, Z = sph2cart(theta, phi)
    r0, theta0, phi0 = cart2sph(X, Y, Z)

    print 'Made it through without any errors.'
