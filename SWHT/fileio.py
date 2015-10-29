"""
functions to read/write data for the imagers
"""

#TODO: hdf5 wrappers for visibilities and images

import cPickle as pkl
import numpy as np
import datetime

def parse(fn, fmt=None):
    """Parse an input visibility filename to determine meta data and type
    XST files are assumed to follow the SE607 format: <date>_<time>_rcu<id>_sb<subband>_int<integration length>_dur<duration of observation>[_<HBA config in hex>]_xst.dat
    fmt: if None then automatically determines format based on filename, else can be set to 'ms' (measurement set), 'acc' (LOFAR ACC), 'xst' (LOFAR XST)
    returns: dictionary"""
    fDict = {}
    fDict['fn'] = fn
    if fn.lower().endswith('.ms') or fmt=='ms':
        fDict['fmt'] = 'ms'
    elif fn.lower().endswith('.dat') or fn.lower().endswith('.dat.sim') or fmt=='acc' or fmt=='xst':
        metaData = fn.split('/')[-1].split('_')
        fDict['ts'] = datetime.datetime(year=int(metaData[0][:4]), month=int(metaData[0][4:6]), day=int(metaData[0][6:]), hour=int(metaData[1][:2]), minute=int(metaData[1][2:4]), second=int(metaData[1][4:]))
        if metaData[2].startswith('acc'): #the file is a LOFAR ACC file
            fDict['fmt'] = 'acc'
            fDict['shape'] = map(int, metaData[3].split('.')[0].split('x'))
        elif metaData[-1].startswith('xst.dat'): #the file is a SE607 format LOFAR XST file
            fDict['fmt'] = 'xst'
            fDict['rcu'] = int(metaData[2][3:])
            fDict['sb'] = np.array( [int(metaData[3][2:])] )
            fDict['int'] = float(metaData[4][3:])
            fDict['dur'] = float(metaData[5][3:])
            if len(metaData)==8: #HBA all-sky file, get element identifiers
                fDict['elem'] = metaData[6][2:]
    elif fn.lower().endswith('.pkl') or fmt=='pkl': #the file is a set of SWHT image coefficients
        fDict['fmt'] = 'pkl'
    else:
        #unknown data format, returns warning
        fDict['fmt'] = -1
    return fDict

def writeCoeffPkl(fn, coeffs, phs=[0., 0.], lst=0.):
    """Write SWHT image coefficients to file
    fn: str, pickle filename
    coeffs: 2D array of complex coefficients
    phs: [float, float], RA and Dec (radians) position at the center of the image
    lst: float, local sidereal time of snapshot
    """
    coeffDict = {
        'coeffs': coeffs,
        'phs': phs,
        'lst': lst
    }
    fh = open(fn, 'wb')
    pkl.dump(coeffDict, fh)
    fh.close()
    
def readCoeffPkl(fn):
    """Read SWHT image coefficients from a pickle file, see writeCoeffPkl() for contents"""
    fh = open(fn,'rb')
    coeffDict = pkl.load(fh)
    fh.close()
    return coeffDict

def writeImgPkl(fn, d, fDict, res=None, fttype=None, imtype=None):
    """Write an image cube to a pickle file
    fn: str, pickle filename
    d: numpy array, image data
    fDict: dict, meta data from original visibility file
    res: float, resolution at zenith (radians)
    fftype: str, dft or fft convolution function name
    imtype: str, complex or Stokes"""
    imgDict = {
        'meta': fDict,
        'res': res,
        'fttype': fttype,
        'imtype': imtype,
        'img': d}
    fh = open(fn, 'wb')
    pkl.dump(imgDict, fh)
    fh.close()

def readImgPkl(fn):
    """Read an image cube from a pickle file, see writeImgPkl() for contents"""
    fh = open(fn,'rb')
    imgDict = pkl.load(fh)
    fh.close()
    return imgDict

def writeSWHTImgPkl(fn, d, fDict, mode):
    """Write a SWHT image cube to a pickle file
    fn: str, pickle filename
    d: numpy array, image data
    fDict: dict, meta data from original visibility file
    """
    imgDict = {
        'meta': fDict,
        'mode': mode}
    if mode.startswith('3D'):
        imgDict['img'] = d[0]
        imgDict['phi'] = d[1]
        imgDict['theta'] = d[2]
    else:
        imgDict['img'] = d
    fh = open(fn, 'wb')
    pkl.dump(imgDict, fh)
    fh.close()

def readSWHTImgPkl(fn):
    """Read an image cube from a pickle file, see writeSWHTImgPkl() for contents"""
    fh = open(fn,'rb')
    imgDict = pkl.load(fh)
    fh.close()
    return imgDict

if __name__ == '__main__':
    print 'Running test cases...'

    fDict = parse('../examples/20150607_122433_acc_512x192x192.dat')
    print fDict

    fDict = parse('../examples/zen.2455819.69771.uvcRREM.MS')
    print fDict

    fDict = parse('../examples/20150915_191137_rcu5_sb60_int10_dur10_elf0f39fe2034ea85fc02b3cc1544863053b328fd83291e880cd0bf3c3d3a50a164a3f3e0c070c73d073f4e43849c0e93b_xst.dat')
    print fDict

    print '...Made it through without errors'

