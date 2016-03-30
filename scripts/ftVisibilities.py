#!/usr/bin/env python
"""
Perform a Fourier Transform (standard or fast) on LOFAR ACC/XST data or widefield MS data (e.g. PAPER) to form a complex or Stokes dirty image dirty image, single files only
"""

import numpy as np
from matplotlib import pyplot as plt
import datetime
import ephem
import casacore.tables as tbls
import sys,os
import SWHT

#import scipy.constants
#cc = scipy.constants.c
cc = 299792458.0 #speed of light, m/s

if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] ACC/XST/MS FILE')
    o.set_description(__doc__)
    o.add_option('--station', dest='station', default=None,
        help = 'LOFAR ONLY: station name, e.g. SE607, if this is used then the ant_field and ant_array options are not required, default: None')
    o.add_option('-F', '--ant_field', dest='ant_field', default=None,
        help = 'LOFAR ONLY: AntennaField.conf file for the LOFAR station of the ACC files, default: None')
    o.add_option('-A', '--ant_array', dest='ant_array', default=None,
        help = 'LOFAR ONLY: AntennaArray.conf file for the LOFAR station geographical coordinates, default: None')
    o.add_option('-D', '--deltas', dest='deltas', default=None,
        help = 'LOFAR ONLY: iHBADeltas.conf file, only required for HBA imaging, default: None')
    o.add_option('-r', '--rcumode', dest='rcumode', default=3, type='int',
        help = 'LOFAR ONLY: Station RCU Mode, usually 3,5,6,7, for XST it will override filename metadata default: 3(LBA High)')
    o.add_option('-s', '--subband', dest='subband', default='0',
        help = 'Select which subband(s) to image, for ACC and MS it will select, for multiple subbands use X,Y,Z and for range use X_Y notation, for XST it will override filename metadata, default:0')
    o.add_option('-p', '--pixels', dest='pixels', default=64, type='int',
        help = 'Width of image in pixels, default: 64')
    o.add_option('-C', '--cal', dest='calfile', default=None,
        help = 'LOFAR ONLY: Apply a calibration soultion file to the data.')
    o.add_option('-S', '--save', dest='savefig', default=None,
        help = 'Save the figure using this name, type is determined by extension')
    o.add_option('--conv', dest='conv', default='fast',
        help = 'If using FFT, choose a convolution function: fast(nearest neighbor), rectangle, gaussian, prolate spheroid. default:fast')
    o.add_option('--dft', dest='dft', action='store_true',
        help = 'Form image with a direct FT instead of an FFT')
    o.add_option('--nodisplay', dest='nodisplay', action='store_true',
        help = 'Do not display the generated image')
    o.add_option('--pkl', dest='pkl', default=None,
        help = 'Save complex images in a numpy array in a pickle file using this name (include .pkl extention), default: tempImage.pkl')
    o.add_option('-i', '--int', dest='int_time', default=1., type='float',
        help = 'LOFAR ONLY: Integration time, used for accurate zenith pointing, for XST it will override filename metadata, default: 1 second')
    o.add_option('-c', '--column', dest='column', default='CORRECTED_DATA', type='str',
        help = 'MS ONLY: select which data column to image, default: CORRECTED_DATA')
    o.add_option('--override', dest='override', action='store_true',
        help = 'LOFAR XST ONLY: override filename metadata for RCU, integration length, and subband')
    o.add_option('--autos', dest='autos', action='store_true',
        help = 'Include the auto-correlation in the image, by default they are blanked')
    o.add_option('--weight', dest='weighting', default='natural',
        help = 'Weighting mode, natural (default), uniform')
    o.add_option('--fov', dest='fov', default=180., type='float',
        help = 'Field of View in degrees, default: 180 (all-sky)')
    opts, args = o.parse_args(sys.argv[1:])

    visFile = args[0]
    fDict = SWHT.fileio.parse(visFile)
    print fDict

    #parse subbands
    sbs = np.array(SWHT.util.convert_arg_range(opts.subband))
    
    #Pull out the visibility data in a (u,v,w) format
    if fDict['fmt']=='acc' or fDict['fmt']=='xst': #LOFAR visibilities
        if fDict['fmt']=='acc' or opts.override:
            fDict['rcu'] = opts.rcumode #add the RCU mode to the meta data of an ACC file, or override the XST metadat
            fDict['sb'] = sbs
            fDict['int'] = opts.int_time
        else:
            sbs = fDict['sb']

        lofarStation = SWHT.lofarConfig.getLofarStation(name=opts.station, affn=opts.ant_field, aafn=opts.ant_array, deltas=opts.deltas) #get station position information

        #longitude and latitude of array
        #lon, lat, elev = lofarStation.antArrays.location[SWHT.lofarConfig.rcuInfo[fDict['rcu']]['array_type']]
        arr_xyz = lofarStation.antField.location[SWHT.lofarConfig.rcuInfo[fDict['rcu']]['array_type']]
        lat, lon, elev = SWHT.ecef.ecef2geodetic(arr_xyz[0], arr_xyz[1], arr_xyz[2], degrees=True)
        print 'LON(deg):', lon, 'LAT(deg):', lat, 'ELEV(m):', elev

        #antenna positions
        ants = lofarStation.antField.antpos[SWHT.lofarConfig.rcuInfo[fDict['rcu']]['array_type']]
        if 'elem' in fDict: #update the antenna positions if there is an element string
            if lofarStation.deltas is None:
                print 'Warning: HBA element string found, but HBADeltas file is missing, your image is probably not going to make sense'
            else:
                print 'Updating antenna positions with HBA element deltas'
                for aid in np.arange(ants.shape[0]):
                    delta = lofarStation.deltas[int(fDict['elem'][aid], 16)]
                    delta = np.array([delta, delta])
                    ants[aid] += delta
        nants = ants.shape[0]
        print 'NANTENNAS:', nants

        #frequency information
        nchan = SWHT.lofarConfig.rcuInfo[fDict['rcu']]['nchan']
        bw = SWHT.lofarConfig.rcuInfo[fDict['rcu']]['bw']
        df = bw/nchan
        freqs = sbs*df + SWHT.lofarConfig.rcuInfo[fDict['rcu']]['offset'] + (df/2.) #df/2 to centre the band
        print 'SUBBANDS:', sbs, '(', freqs/1e6, 'MHz)'
        npols = 2
        
        #read LOFAR Calibration Table
        if not (opts.calfile is None):
            print 'Using CalTable:', opts.calfile
            antGains = SWHT.lofarConfig.readCalTable(opts.calfile, nants, nchan, npols)
        else: antGains = None

        #get correlation matrix for subbands selected
        nantpol = nants * npols
        print 'Reading in visibility data file ...',
        if fDict['fmt']=='acc':
            tDeltas = [] #subband timestamp deltas from the end of file
            corrMatrix = np.fromfile(visFile, dtype='complex').reshape(nchan, nantpol, nantpol) #read in the complete correlation matrix
            sbCorrMatrix = np.zeros((sbs.shape[0], nantpol, nantpol), dtype=complex)
            for sbIdx, sb in enumerate(sbs):
                if antGains is None:
                    sbCorrMatrix[sbIdx] = corrMatrix[sb, :, :] #select out a single subband, shape (nantpol, nantpol)
                else: #Apply Gains
                    sbAntGains = antGains[sb][np.newaxis].T
                    sbVisGains = np.conjugate(np.dot(sbAntGains, sbAntGains.T)) # from Tobia, visibility gains are computed as (G . G^T)*
                    sbCorrMatrix[sbIdx] = np.multiply(sbVisGains, corrMatrix[sb, :, :]) #select out a single subband, shape (nantpol, nantpol)

                #correct the time due to subband stepping
                tOffset = (nchan - sb) * fDict['int'] #the time stamp in the filename in for the last subband
                rem = tOffset - int(tOffset) #subsecond remainder
                tDeltas.append(datetime.timedelta(0, int(tOffset), rem*1e6))
            meants = fDict['ts'] - SWHT.util.meanTimeDelta(tDeltas) #if using multiple subbands, use the mean offset time
            
        elif fDict['fmt']=='xst':
            corrMatrix = np.fromfile(visFile, dtype='complex').reshape(1, nantpol, nantpol) #read in the correlation matrix
            if antGains is None:
                sbCorrMatrix = corrMatrix #shape (nantpol, nantpol)
            else: #Apply Gains
                sbAntGains = antGains[fDict['sb']][np.newaxis].T
                sbVisGains = np.conjugate(np.dot(sbAntGains, sbAntGains.T)) # from Tobia, visibility gains are computed as (G . G^T)*
                sbCorrMatrix = np.multiply(sbVisGains, corrMatrix) #shape (nantpol, nantpol)
            meants = fDict['ts']

        print 'done'
        print 'CORRELATION MATRIX SHAPE', corrMatrix.shape
        
        obs = ephem.Observer() #create an observer at the array location
        obs.long = lon * (np.pi/180.)
        obs.lat = lat * (np.pi/180.)
        obs.elevation = float(elev)
        obs.epoch = meants
        obs.date = meants
        print 'Observatory:', obs

        src = ephem.FixedBody() #create a source at zenith
        src._ra = obs.sidereal_time()
        src._dec = obs.lat
        src.compute(obs)
        
        #get antenna positions in ITRF (x,y,z) format and compute the (u,v,w) coordinates pointing at zenith
        xyz = []
        for a in ants: xyz.append([a[0,0]+arr_xyz[0], a[0,1]+arr_xyz[1], a[0,2]+arr_xyz[2]])
        xyz = np.array(xyz)
        uvw = SWHT.ft.xyz2uvw(xyz, src, obs, freqs)
        uu = SWHT.util.vectorize3D(uvw[:,:,:,0]).flatten()
        vv = SWHT.util.vectorize3D(uvw[:,:,:,1]).flatten()
        ww = SWHT.util.vectorize3D(uvw[:,:,:,2]).flatten()
        uvw = np.vstack((uu, vv, ww)).T

        #split up polarizations, vectorize the correlation matrix, and drop the lower triangle
        xxVis = SWHT.util.vectorize3D(sbCorrMatrix[:, 0::2, 0::2]).flatten()
        xyVis = SWHT.util.vectorize3D(sbCorrMatrix[:, 0::2, 1::2]).flatten()
        yxVis = SWHT.util.vectorize3D(sbCorrMatrix[:, 1::2, 0::2]).flatten()
        yyVis = SWHT.util.vectorize3D(sbCorrMatrix[:, 1::2, 1::2]).flatten()

        ##uv coverage plot
        #plt.plot(uvw[:,0], uvw[:,1], '.')
        #plt.show()

    elif fDict['fmt']=='ms': #MS-based visibilities

        fDict['sb'] = sbs

        MS = tbls.table(visFile, readonly=True)
        data_column = opts.column.upper()
        uvw = MS.col('UVW').getcol() # [vis id, (u,v,w)]
        vis = MS.col(data_column).getcol() #[vis id, freq id, stokes id]
        vis = vis[:,fDict['sb'],:] #select a single subband
        MS.close()

        #freq information, convert uvw coordinates
        SW = tbls.table(visFile + '/SPECTRAL_WINDOW')
        freqs = SW.col('CHAN_FREQ').getcol()[0, sbs][np.newaxis] # [1, nchan]
        #convert (u,v,w) from metres to wavelengths
        uu = np.dot(uvw[:,0][np.newaxis].T, freqs).flatten()
        vv = np.dot(uvw[:,1][np.newaxis].T, freqs).flatten()
        ww = np.dot(uvw[:,2][np.newaxis].T, freqs).flatten()
        uvw = np.vstack((uu, vv, ww)).T / cc
        print 'SUBBANDS:', sbs, '(', freqs/1e6, 'MHz)'
        SW.close()

        #split up polarizations
        print vis.shape
        xxVis = vis[:,:,0].flatten() 
        xyVis = vis[:,:,1].flatten()
        yxVis = vis[:,:,2].flatten()
        yyVis = vis[:,:,3].flatten()
        print xxVis.shape

        ##uv coverage plot
        #plt.plot(uvw[:,0], uvw[:,1], '.')
        #plt.show()

    else:
        print 'ERROR: unknown data format, exiting'
        exit()
    
    #remove auto-correlations
    print 'AUTO-CORRELATIONS:', opts.autos
    if not opts.autos:
        autoIdx = np.argwhere(uvw[:,0]**2. + uvw[:,1]**2. + uvw[:,2]**2. == 0.)
        xxVis[autoIdx] = 0.
        xyVis[autoIdx] = 0.
        yxVis[autoIdx] = 0.
        yyVis[autoIdx] = 0.

    #prepare for Fourier transform
    print 'Performing Fourier Transform'
    pixels = opts.pixels
    px = [pixels,pixels]
    fov = opts.fov * (np.pi/180.) #Field of View in radians
    res = fov / px[0] #pixel resolution
    print 'Resolution(deg):', res*180./np.pi

    #perform DFT or FFT
    if opts.dft:
        print 'DFT'
        xxIm = SWHT.ft.dftImage(xxVis, uvw, px, res, mask=False, rescale=False, stokes=False)
        xyIm = SWHT.ft.dftImage(xyVis, uvw, px, res, mask=False, rescale=False, stokes=False)
        yxIm = SWHT.ft.dftImage(yxVis, uvw, px, res, mask=False, rescale=False, stokes=False)
        yyIm = SWHT.ft.dftImage(yyVis, uvw, px, res, mask=False, rescale=False, stokes=False)
    else:
        print 'FFT'
        conv = opts.conv
        xxIm = SWHT.ft.fftImage(xxVis, uvw, px, res, mask=False, wgt=opts.weighting, conv=conv)
        xyIm = SWHT.ft.fftImage(xyVis, uvw, px, res, mask=False, wgt=opts.weighting, conv=conv)
        yxIm = SWHT.ft.fftImage(yxVis, uvw, px, res, mask=False, wgt=opts.weighting, conv=conv)
        yyIm = SWHT.ft.fftImage(yyVis, uvw, px, res, mask=False, wgt=opts.weighting, conv=conv)

    #save complex image to pickle file
    if opts.pkl is None: outPklFn = 'tempImage.pkl'
    else: outPklFn = opts.pkl
    if opts.dft: fttype = 'dft'
    else: fttype = opts.conv
    print 'Writing image to file %s ...'%outPklFn,
    SWHT.fileio.writeImgPkl(outPklFn, np.array([xxIm,xyIm,yxIm,yyIm]), fDict, res=res, fttype=fttype, imtype='complex')
    print 'done'
    
    #display stokes plots
    if not opts.nodisplay or not (opts.savefig is None):
        #generate stokes images
        iIm = (xxIm + yyIm).real
        qIm = (xxIm - yyIm).real
        uIm = (xyIm + yxIm).real
        vIm = (yxIm - xyIm).imag
    
        plt.subplot(2,2,1)
        plt.imshow(iIm)
        plt.xlabel('Pixels (E-W)')
        plt.ylabel('Pixels (N-S)')
        plt.title('I')
        plt.colorbar()
        plt.subplot(2,2,2)
        plt.imshow(qIm)
        plt.xlabel('Pixels (E-W)')
        plt.ylabel('Pixels (N-S)')
        plt.title('Q')
        plt.colorbar()
        plt.subplot(2,2,3)
        plt.imshow(uIm)
        plt.xlabel('Pixels (E-W)')
        plt.ylabel('Pixels (N-S)')
        plt.title('U')
        plt.colorbar()
        plt.subplot(2,2,4)
        plt.imshow(vIm)
        plt.xlabel('Pixels (E-W)')
        plt.ylabel('Pixels (N-S)')
        plt.title('V')
        plt.colorbar()
    if not (opts.savefig is None): plt.savefig(opts.savefig)
    if not opts.nodisplay: plt.show()

