"""
An implementation on spherical harmonics in python becasue scipy.special.sph_harm in scipy<=0.13 is very slow

Originally written by Jozef Vesely
https://github.com/scipy/scipy/issues/1280
"""

import numpy as np

def xfact(m):
    # computes (2m-1)!!/sqrt((2m)!)
    res = 1.
    for i in xrange(1, 2*m+1):
        if i % 2: res *= i # (2m-1)!!
        res /= np.sqrt(i) # sqrt((2m)!)
    return res

def lplm_n(l, m, x):
    # associated legendre polynomials normalized as in Ylm, from Numerical Recipes 6.7
    l,m = int(l),int(m)
    assert 0<=m<=l and np.all(np.abs(x)<=1.)

    norm = np.sqrt(2. * l + 1.) / np.sqrt(4. * np.pi)
    if m == 0:
        pmm = norm * np.ones_like(x)
    else:
        pmm = (-1.)**m * norm * xfact(m) * (1.-x**2.)**(m/2.)
    if l == m:
        return pmm
    pmmp1 = x * pmm * np.sqrt(2.*m+1.)
    if l == m+1:
        return pmmp1
    for ll in xrange(m+2, l+1):
        pll = (x*(2.*ll-1.)*pmmp1 - np.sqrt( (ll-1.)**2. - m**2.)*pmm)/np.sqrt(ll**2.-m**2.)
        pmm = pmmp1
        pmmp1 = pll
    return pll

def Ylm(l, m, phi, theta):
    # spherical harmonics
    # theta is from 0 to pi with pi/2 on equator
    l,m = int(l),int(m)
    assert 0 <= np.abs(m) <=l
    if m > 0:
        return lplm_n(l, m, np.cos(theta)) * np.exp(1J * m * phi)
    elif m < 0:
        return (-1.)**m * lplm_n(l, -m, np.cos(theta)) * np.exp(1J * m * phi)
    return lplm_n(l, m, np.cos(theta)) * np.ones_like(phi)

def Ylmr(l, m, phi, theta):
    # real spherical harmonics
    # theta is from 0 to pi with pi/2 on equator
    l,m = int(l),int(m)
    assert 0 <= np.abs(m) <=l
    if m > 0:
        return lplm_n(l, m, np.cos(theta)) * np.cos(m * phi) * np.sqrt(2.)
    elif m < 0:
        return (-1.)**m * lplm_n(l, -m, np.cos(theta)) * np.sin(-m * phi) * np.sqrt(2.)
    return lplm_n(l, m, np.cos(theta)) * np.ones_like(phi)

if __name__ == "__main__":
    from scipy.special import sph_harm
    from scipy.misc import factorial2, factorial
    from timeit import Timer

    def ref_xfact(m):
        return factorial2(2*m-1)/np.sqrt(factorial(2*m))
 
    print "Time: xfact(10)", Timer("xfact(10)",
        "from __main__ import xfact, ref_xfact").timeit(100)
    print "Time: ref_xfact(10)", Timer("ref_xfact(10)",
        "from __main__ import xfact, ref_xfact").timeit(100)
    print "Time: xfact(80)", Timer("xfact(80)",
        "from __main__ import xfact, ref_xfact").timeit(100)
    print "Time: ref_xfact(80)", Timer("ref_xfact(80)",
        "from __main__ import xfact, ref_xfact").timeit(100)
    
    print "m", "xfact", "ref_xfact" 
    for m in range(10) + range(80,90):
        a = xfact(m)
        b = ref_xfact(m)
        print m, a, b 

    phi, theta = np.ogrid[0:2*np.pi:10j,-np.pi/2:np.pi/2:10j]

    print "Time: Ylm(1,1,phi,theta)", Timer("Ylm(1,1,phi,theta)",
        "from __main__ import Ylm, sph_harm, phi, theta").timeit(10)
    print "Time: sph_harm(1,1,phi,theta)", Timer("sph_harm(1,1,phi,theta)",
        "from __main__ import Ylm, sph_harm, phi, theta").timeit(10)

    
    print "l", "m", "max|Ylm-sph_harm|" 
    for l in xrange(0,10):
        for m in xrange(-l,l+1):
            a = Ylm(l,m,phi,theta) 
            b = sph_harm(m,l,phi,theta)
            print l,m, np.amax(np.abs(a-b))

