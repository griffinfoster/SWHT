#!/usr/bin/env python
"""
Plot the UVW coverage/sampling 
"""

import numpy as np
import matplotlib.pyplot as plt
import datetime
import sys,os
import SWHT

#import scipy.constants
#cc = scipy.constants.c
cc = 299792458.0 #speed of light, m/s

# TODO: plot as a function of wavelength
# TODO: the rotation is off, need to sort this out

if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] ACC/XST/MS/PKL FILE')
    o.set_description(__doc__)
    o.add_option('--station', dest='station', default=None,
        help = 'LOFAR ONLY: station name, e.g. SE607, if this is used then the ant_field and ant_array options are not required, default: None')
    o.add_option('-F', '--ant_field', dest='ant_field', default=None,
        help = 'LOFAR ONLY: AntennaField.conf file for the LOFAR station of the ACC files, default: None')
    o.add_option('-A', '--ant_array', dest='ant_array', default=None,
        help = 'LOFAR ONLY(NOT REQUIRED): AntennaArray.conf file for the LOFAR station geographical coordinates, default: None')
    o.add_option('-D', '--deltas', dest='deltas', default=None,
        help = 'LOFAR ONLY: iHBADeltas.conf file, only required for HBA imaging, default: None')
    o.add_option('-r', '--rcumode', dest='rcumode', default=3, type='int',
        help = 'LOFAR ONLY: Station RCU Mode, usually 3,5,6,7, for XST it will override filename metadata default: 3(LBA High)')
    o.add_option('-s', '--subband', dest='subband', default='0',
        help = 'Select which subband(s) to image, for ACC and MS it will select, for multiple subbands use X,Y,Z and for range use X_Y notation, for XST it will override filename metadata, default:0')
    o.add_option('-C', '--cal', dest='calfile', default=None,
        help = 'LOFAR ONLY: Apply a calibration soultion file to the data.')
    o.add_option('-S', '--save', dest='savefig', default=None,
        help = 'Save the figure using this name, type is determined by extension')
    o.add_option('--nodisplay', dest='nodisplay', action='store_true',
        help = 'Do not display the generated image')
    o.add_option('-i', '--int', dest='int_time', default=1., type='float',
        help = 'LOFAR ONLY: Integration time, used for accurate zenith pointing, for XST it will override filename metadata, default: 1 second')
    o.add_option('-c', '--column', dest='column', default='CORRECTED_DATA', type='str',
        help = 'MS ONLY: select which data column to image, default: CORRECTED_DATA')
    o.add_option('--override', dest='override', action='store_true',
        help = 'LOFAR XST ONLY: override filename metadata for RCU, integration length, and subband')
    o.add_option('-m', '--mode', dest='mode', default='3D',
        help = 'UV coverage mode, 3D: UVW coverage, 2D: UV coverage, default: 3D')
    o.add_option('-t', '--times', dest='times', default='0',
        help = 'KAIRA ONLY: Select which integration(s) to image, can use a[seconds] to average, d[step size] to decimate, of a specific range of integrations similar to the subband selection option, default:0 (select the first integration of the file)')
    opts, args = o.parse_args(sys.argv[1:])

    # parse subbands
    sbs = np.array(SWHT.util.convert_arg_range(opts.subband))

    # setup variables for combined visibilities and uvw samples
    visComb = np.array([]).reshape(4, 0, len(sbs))
    uvwComb = np.array([]).reshape(0, 3, len(sbs))

    dataFmt = None
    if (not (opts.station is None)) or (not (opts.ant_field is None)): # If using LOFAR data, get station information
        lofarStation = SWHT.lofarConfig.getLofarStation(name=opts.station, affn=opts.ant_field, aafn=opts.ant_array, deltas=opts.deltas)
        antGains = None # Setup variable so that the gain table isn't re-read for every file if used
        if lofarStation.name=='KAIRA': dataFmt='KAIRA'

    ####################
    ## Read Visibilities
    ####################
    visFiles = args # filenames to image
    for vid,visFn in enumerate(visFiles):
        print 'Using %s (%i/%i)'%(visFn, vid+1, len(visFiles))
        fDict = SWHT.fileio.parse(visFn, fmt=dataFmt)

        # Pull out the visibility data in a (u,v,w) format
        if fDict['fmt']=='acc': # LOFAR station all subbands ACC file visibilities
            decomp = True
            fDict['rcu'] = opts.rcumode # add the RCU mode to the meta data of an ACC file
            fDict['sb'] = sbs # select subbands to use
            fDict['int'] = opts.int_time # set integration length (usually 1 second)

            vis, uvw, freqs, obsInfo = SWHT.fileio.readACC(visFn, fDict, lofarStation, sbs, calTable=opts.calfile)
            [obsLat, obsLong, LSTangle] = obsInfo

            # add visibilities to previously processed files
            visComb = np.concatenate((visComb, vis), axis=1)
            uvwComb = np.concatenate((uvwComb, uvw), axis=0)

        elif fDict['fmt']=='xst': # LOFAR XST format visibilities
            decomp = True
            if opts.override: # Override XST filename metadata
                fDict['rcu'] = opts.rcumode
                fDict['sb'] = sbs
                fDict['int'] = opts.int_time
            else:
                sbs = fDict['sb']

            vis, uvw, freqs, obsInfo = SWHT.fileio.readXST(visFn, fDict, lofarStation, sbs, calTable=opts.calfile)
            [obsLat, obsLong, LSTangle] = obsInfo

            # add visibilities to previously processed files
            visComb = np.concatenate((visComb, vis), axis=1)
            uvwComb = np.concatenate((uvwComb, uvw), axis=0)

        elif fDict['fmt']=='KAIRA': # KAIRA LOFAR XST format visibilities
            decomp = True
            if opts.override: # Override XST filename metadata
                fDict['rcu'] = opts.rcumode
                fDict['sb'] = sbs
                fDict['int'] = opts.int_time
            else:
                sbs = fDict['sb']

            vis, uvw, freqs, obsInfo = SWHT.fileio.readKAIRAXST(visFn, fDict, lofarStation, sbs, times=opts.times)
            [obsLat, obsLong, LSTangle] = obsInfo

            # add visibilities to previously processed files
            visComb = np.concatenate((visComb, vis), axis=1)
            uvwComb = np.concatenate((uvwComb, uvw), axis=0)

        elif fDict['fmt']=='ms': # MS-based visibilities
            decomp = True

            fDict['sb'] = sbs

            vis, uvw, freqs, obsInfo = SWHT.fileio.readMS(visFn, sbs, column=opts.column)
            [obsLat, obsLong, LSTangle] = obsInfo

            # add visibilities to previously processed files
            visComb = np.concatenate((visComb, vis), axis=1)
            uvwComb = np.concatenate((uvwComb, uvw), axis=0)

        elif fDict['fmt']=='pkl':
            print 'Loading Image Coefficients file:', visFn
            coeffDict = SWHT.fileio.readCoeffPkl(visFn)
            iImgCoeffs = coeffDict['coeffs']
            LSTangle = coeffDict['lst']
            obsLong = coeffDict['phs'][0]
            obsLat = coeffDict['phs'][1]
            decomp = False

        else:
            print 'ERROR: unknown data format, exiting'
            exit()

    if opts.mode.lower() == '3d':
        fig, ax = SWHT.display.dispVis3D(uvwComb)
    else:
        fig, ax = SWHT.display.dispVis2D(uvwComb)

    if not (opts.savefig is None): plt.savefig(opts.savefig)
    if not opts.nodisplay: plt.show()

