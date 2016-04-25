"""
Functions to display images and coefficients
"""

import matplotlib.pyplot as plt
import numpy as np
import swht

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

