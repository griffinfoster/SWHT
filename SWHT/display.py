"""
Functions to display images and coefficients
"""

import matplotlib.pyplot as plt
import numpy as np
import swht, util

def disp3D(img, phi, theta, dmode='abs', cmap='jet'):
    """Display 3D, equal in phi and theta (Driscoll and Healy mapping) image
    img: 2D array of image values
    phi: 2D array of phi values
    theta: 2D array of thetat values
        img, phi, theta are of the same shape, they are the output of swht.make3Dimage()
    dim: 2 element list of resolution in phi and theta [delta phi, delta theta]
    dmode: string, data mode (abs, real, imaginary, phase)
    cmap: string, matplotlib colormap name
    """
    if dmode.startswith('abs'): img = np.abs(img)
    elif dmode.startswith('real'): img = img.real
    elif dmode.startswith('imag'): img = img.imag
    elif dmode.startswith('phase'): img = np.angle(img)
    else:
        print 'WARNING: Unknown data mode, defaulting to absolute value'
        img = np.abs(img)

    #X = np.cos(theta-(np.pi/2.)) * np.cos(phi)
    #Y = np.cos(theta-(np.pi/2.)) * np.sin(phi)
    #Z = np.sin(theta-(np.pi/2.))
    X, Y, Z = util.sph2cart(theta, phi)

    #http://stackoverflow.com/questions/22175533/what-is-the-equivalent-of-matlabs-surfx-y-z-c-in-matplotlib
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    from matplotlib.colors import Normalize
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    imin = img.min()
    imax = img.max()
    scalarMap = cm.ScalarMappable(norm=Normalize(vmin=imin, vmax=imax), cmap=plt.get_cmap(cmap))
    #scalarMap = cm.ScalarMappable(norm=Normalize(vmin=imin, vmax=imax), cmap=cm.jet)
    #scalarMap = cm.ScalarMappable(norm=Normalize(vmin=imin, vmax=imax), cmap=cm.gist_earth_r)
    scalarMap.set_array(img)
    C = scalarMap.to_rgba(img)
    surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors=C, antialiased=True)
    fig.colorbar(scalarMap)
    return fig, ax

def dispCoeffs(imgCoeffs, zeroDC=True, vis=False):
    """Display SWHT image coefficient values in a 2 rows x 3 columns plot
    imgCoeffs: SWHT image coefficients array
    zeroDC: zero out DC coefficient values
    vis: convert image coefficients to visibility coefficients

    returns: matplotlib figure and subplots
    """
    fig, ax = plt.subplots(2, 3)

    if vis: imgCoeffs = swht.computeblm(imgCoeffs, reverse=True) #convert brightness coefficients to visibility coefficients
    
    if zeroDC: imgCoeffs[0,0] = 0 #zero out DC offset component

    # Top Left
    plt.axes(ax[0,0])
    plt.imshow(imgCoeffs.real, interpolation='nearest')
    plt.title('Real Components')
    plt.colorbar()

    # Top Centre
    plt.axes(ax[0,1])
    plt.title('Imaginary Components')
    plt.imshow(imgCoeffs.imag, interpolation='nearest')
    plt.colorbar()

    # Bottom Left
    plt.axes(ax[1,0])
    plt.title('Amplitude (dB)')
    plt.imshow(10.*np.log10(np.abs(imgCoeffs)), interpolation='nearest')
    plt.colorbar()

    # Bottom Centre
    plt.axes(ax[1,1])
    plt.title('Phase')
    plt.imshow(np.angle(imgCoeffs), interpolation='nearest')
    plt.colorbar()

    # Top Right
    plt.axes(ax[0,2])
    coeffsFlat = []
    mms = []
    lls = []
    for ll in np.arange(imgCoeffs.shape[0]):
        for mm in np.arange(-1*ll, ll+1):
            mms.append(mm)
            lls.append(ll)
            coeffsFlat.append(imgCoeffs[ll,ll+mm])
    coeffsFlat = np.array(coeffsFlat)
    plt.ylabel('Amplitude (dB)')
    plt.xlabel('l')
    plt.plot(lls, 10.*np.log10(np.abs(coeffsFlat)), '.')

    # Bottom Right
    plt.axes(ax[1,2])
    plt.ylabel('Amplitude (dB)')
    plt.xlabel('m')
    plt.plot(mms, 10.*np.log10(np.abs(coeffsFlat)), '.')

    return fig, ax

