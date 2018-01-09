"""
Functions to display images and coefficients
"""

import matplotlib.pyplot as plt
import matplotlib.patches
import numpy as np
import swht, util

def disp2D(img, dmode='dB', cmap='jet'):
    """Display 2D projected image
    img: 2D array of complex flux values
    dmode: string, data mode (abs, dB (absolute value in log units), real, imaginary, phase)
    cmap: string, matplotlib colormap name
    """
    if dmode.startswith('abs'): img = np.abs(img)
    elif dmode.startswith('dB'): img = 10. * np.log10(np.abs(img))
    elif dmode.startswith('real'): img = img.real
    elif dmode.startswith('imag'): img = img.imag
    elif dmode.startswith('phase'): img = np.angle(img)
    else:
        print 'WARNING: Unknown data mode, defaulting to absolute value'
        img = np.abs(img)

    img = np.fliplr(img)

    fig, ax = plt.subplots(1, 1)
    ax.yaxis.set_ticks([]) # hide tick marks and labels
    ax.xaxis.set_ticks([]) # hide tick marks and labels
    ax.patch.set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.imshow(img, interpolation='nearest', cmap=plt.get_cmap(cmap))

    # draw alt-az grid (assume square image)
    xc = img.shape[0]/2. - 0.5
    yc = img.shape[1]/2. - 0.5
    ax.add_patch( matplotlib.patches.Circle((xc, yc), img.shape[0]/2., fill=False)) # alt 0
    altLines = 5 # number of lines at constant altitude, including alt=0 and alt=90
    deltaAlt = img.shape[0] / (2.*altLines)
    for al in range(1, altLines): # skip alt=0 and alt=90
        ax.add_patch( matplotlib.patches.Circle((xc, yc), xc-deltaAlt*al, fill=False, linestyle='dotted', alpha=0.7))
    azLines = 6 # number of constant azimuth lines
    deltaAz = np.pi / azLines
    for az in range(azLines):
        plt.plot(np.array([np.cos(az*deltaAz), np.cos(az*deltaAz + np.pi)])/2. + 0.5, np.array([np.sin(az*deltaAz), np.sin(az*deltaAz + np.pi)])/2. + 0.5, 'k:', alpha=0.7, transform=ax.transAxes)

    # grid labels
    for az in range(2*azLines):
        #plt.text(img.shape[0] * (1.06 * np.cos(az*deltaAz)/2. + 0.5), img.shape[0] * (1.06 * np.sin(az*deltaAz)/2. + 0.5), '%.0f'%(az*deltaAz*180./np.pi), horizontalalignment='center')
        plt.text(img.shape[0] * (1.06 * np.sin(az*deltaAz - np.pi)/2. + 0.5), img.shape[0] * (1.06 * np.cos(az*deltaAz - np.pi)/2. + 0.5), '%.0f'%(az*deltaAz*180./np.pi), horizontalalignment='center')

    plt.colorbar()

    return fig, ax

def disp2DStokes(xx, xy, yx, yy, cmap='jet'):
    """Display 2D projected Stokes images
    xx, xy, yx, yy: 2D array of complex flux values
    cmap: string, matplotlib colormap name
    """
    iIm = (xx + yy).real
    qIm = (xx - yy).real
    uIm = (xy + yx).real
    vIm = (yx - xy).imag

    fig, ax = plt.subplots(2, 2, figsize=(10,8))

    # Top Left
    plt.axes(ax[0,0])
    plt.imshow(iIm)
    #plt.xlabel('Pixels (E-W)')
    plt.ylabel('Pixels (N-S)')
    plt.colorbar()
    plt.title('I')

    # Top Right
    plt.axes(ax[0,1])
    plt.imshow(qIm)
    #plt.xlabel('Pixels (E-W)')
    #plt.ylabel('Pixels (N-S)')
    plt.colorbar()
    plt.title('Q')

    # Bottom Left
    plt.axes(ax[1,0])
    plt.imshow(uIm)
    plt.xlabel('Pixels (E-W)')
    plt.ylabel('Pixels (N-S)')
    plt.colorbar()
    plt.title('U')

    # Bottom Right
    plt.axes(ax[1,1])
    plt.imshow(vIm)
    plt.xlabel('Pixels (E-W)')
    #plt.ylabel('Pixels (N-S)')
    plt.colorbar()
    plt.title('V')

    return fig, ax

# TODO: Add RA/Dec grid
def disp3D(img, phi, theta, dmode='abs', cmap='jet'):
    """Display 3D, equal in phi and theta (Driscoll and Healy mapping) image
    img: 2D array of complex flux values
    phi: 2D array of phi values
    theta: 2D array of theta values
        img, phi, theta are of the same shape, they are the output of swht.make3Dimage()
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

    # North-South Flip
    #Z *= -1.

    #http://stackoverflow.com/questions/22175533/what-is-the-equivalent-of-matlabs-surfx-y-z-c-in-matplotlib
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    from matplotlib.colors import Normalize

    fig = plt.figure(figsize=(10,8))
    ax = fig.gca(projection='3d')
    imin = img.min()
    imax = img.max()

    scalarMap = cm.ScalarMappable(norm=Normalize(vmin=imin, vmax=imax), cmap=plt.get_cmap(cmap))
    #scalarMap = cm.ScalarMappable(norm=Normalize(vmin=imin, vmax=imax), cmap=cm.jet)
    #scalarMap = cm.ScalarMappable(norm=Normalize(vmin=imin, vmax=imax), cmap=cm.gist_earth_r)
    scalarMap.set_array(img)
    C = scalarMap.to_rgba(img)
    surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors=C, antialiased=True, linewidth=1)
    #wire = ax.plot_wireframe(1.01*X, 1.01*Y, 1.01*Z, rstride=5, cstride=5, color='black') # RA/Dec grid

    fig.colorbar(scalarMap, shrink=0.7)
    ax.set_axis_off()
    ax.set_xlim(-0.75,0.75)
    ax.set_ylim(-0.75,0.75)
    ax.set_zlim(-0.75,0.75)
    ax.set_aspect('auto')

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

def dispVis3D(uvw):
    # TODO: sometimes wavelength, sometimes metres
    """3D plot of UVW coverage/sampling
    uvw: (N, 3) array of UVW coordinates

    returns: matplotlib figure and subplots
    """
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure(figsize=(8,6))
    ax = fig.gca(projection='3d')
    
    ax.scatter(uvw[:,0], uvw[:,1], uvw[:,2], edgecolors=None, alpha=0.5)
    ax.set_xlabel('U (m)')
    ax.set_ylabel('V (m)')
    ax.set_zlabel('W (m)')

    return fig, ax

def dispVis2D(uvw):
    # TODO: sometimes wavelength, sometimes metres
    """2D plot of UV coverage/sampling
    uvw: (N, 3) array of UVW coordinates

    returns: matplotlib figure and subplots
    """
    fig, ax = plt.subplots(2, 1, figsize=(6,9))
    ax[0].plot(uvw[:,0], uvw[:,1], '.', alpha=0.5)
    ax[0].set_xlabel('U ($\lambda$)')
    ax[0].set_ylabel('V ($\lambda$)')
    ax[0].set_aspect('equal', 'datalim')

    ax[1].plot(np.sqrt(uvw[:,0]**2. + uvw[:,1]**2.), uvw[:,2], '.', alpha=0.5)
    ax[1].set_xlabel('UVdist ($\lambda$)')
    ax[1].set_ylabel('W ($\lambda$)')

    return fig, ax

