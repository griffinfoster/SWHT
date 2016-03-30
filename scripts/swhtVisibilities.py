#!/usr/bin/env python
"""
Perform a Spherical Wave Harmonic Transform on LOFAR ACC/XST data or widefield MS data (e.g. PAPER) to form a complex or Stokes dirty image dirty image
"""

#TODO: 3D, HEALPix mask

import numpy as np
from matplotlib import pyplot as plt
import datetime
import ephem
import casacore.tables as tbls
import healpy as hp
import sys,os
import SWHT

#import scipy.constants
#cc = scipy.constants.c
cc = 299792458.0 #speed of light, m/s

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
    o.add_option('-p', '--pixels', dest='pixels', default=64, type='int',
        help = 'Width of 2D image in pixels, or number of steps in 3D image, or NSIDE in HEALPix image, default: 64')
    o.add_option('-C', '--cal', dest='calfile', default=None,
        help = 'LOFAR ONLY: Apply a calibration soultion file to the data.')
    o.add_option('-S', '--save', dest='savefig', default=None,
        help = 'Save the figure using this name, type is determined by extension')
    o.add_option('--nodisplay', dest='nodisplay', action='store_true',
        help = 'Do not display the generated image')
    o.add_option('--of', dest='of', default=None,
        help = 'Save complex images in a numpy array in a pickle file or HEALPix map using this name (include .pkl or .hpx extention), default: tempImage.pkl')
    o.add_option('-i', '--int', dest='int_time', default=1., type='float',
        help = 'LOFAR ONLY: Integration time, used for accurate zenith pointing, for XST it will override filename metadata, default: 1 second')
    o.add_option('-c', '--column', dest='column', default='CORRECTED_DATA', type='str',
        help = 'MS ONLY: select which data column to image, default: CORRECTED_DATA')
    o.add_option('--override', dest='override', action='store_true',
        help = 'LOFAR XST ONLY: override filename metadata for RCU, integration length, and subband')
    o.add_option('--autos', dest='autos', action='store_true',
        help = 'Include the auto-correlations in the image')
    o.add_option('--fov', dest='fov', default=180., type='float',
        help = '2D IMAGING MODE ONLY: Field of View in degrees, default: 180 (all-sky)')
    o.add_option('-l', '--lmax', dest='lmax', default=32, type='int',
        help = 'Maximum l spherical harmonic quantal number, rule-of-thumb: lmax ~ (pi/longest baseline resolution), default: 32')
    o.add_option('--lmin', dest='lmin', default=0, type='int',
        help = 'Minimum l spherical harmonic quantal number, usually left as 0, default: 0')
    o.add_option('--ocoeffs', dest='ocoeffs', default=None,
        help = 'Save output image coefficients to a pickle file using this name (include .pkl extention), default: tempCoeffs.pkl')
    o.add_option('-I', '--image', dest='imageMode', default='2D',
        help='Imaging mode: 2D (hemisphere flattened), 3D, healpix, coeff (coefficients) default: 2D')
    opts, args = o.parse_args(sys.argv[1:])

    #parse subbands
    sbs = np.array(SWHT.util.convert_arg_range(opts.subband))

    #setup variables for combined visibilities and uvw positions
    xxVisComb = np.array([]).reshape(0, len(sbs))
    xyVisComb = np.array([]).reshape(0, len(sbs))
    yxVisComb = np.array([]).reshape(0, len(sbs))
    yyVisComb = np.array([]).reshape(0, len(sbs))
    uvwComb = np.array([]).reshape(0, 3, len(sbs))

    #from mpl_toolkits.mplot3d import Axes3D
    #fig = plt.figure()
    #ax = fig.add_subplot(111, projection='3d')
    #clr = ['r', 'g', 'b', 'k', 'y', 'm']

    if (not (opts.station is None)) or (not (opts.ant_field is None)): #If using LOFAR data, get station information
        lofarStation = SWHT.lofarConfig.getLofarStation(name=opts.station, affn=opts.ant_field, aafn=opts.ant_array, deltas=opts.deltas)
        antGains = None #Setup variable so that the gain table isn't re-read for every file if used

    #get filenames to image
    visFiles = args
    for vid,visFn in enumerate(visFiles):
        print 'Using %s (%i/%i)'%(visFn, vid+1, len(visFiles))
        fDict = SWHT.fileio.parse(visFn)

        #Pull out the visibility data in a (u,v,w) format
        if fDict['fmt']=='acc' or fDict['fmt']=='xst': #LOFAR visibilities
            decomp = True
            if fDict['fmt']=='acc' or opts.override:
                fDict['rcu'] = opts.rcumode #add the RCU mode to the meta data of an ACC file, or override the XST metadat
                fDict['sb'] = sbs
                fDict['int'] = opts.int_time
            else:
                sbs = fDict['sb']

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
                if antGains is None: #read the Cal Table only once
                    print 'Using CalTable:', opts.calfile
                    antGains = SWHT.lofarConfig.readCalTable(opts.calfile, nants, nchan, npols)
            else: antGains = None

            #get correlation matrix for subbands selected
            nantpol = nants * npols
            print 'Reading in visibility data file ...',
            if fDict['fmt']=='acc':
                tDeltas = [] #subband timestamp deltas from the end of file
                corrMatrix = np.fromfile(visFn, dtype='complex').reshape(nchan, nantpol, nantpol) #read in the complete correlation matrix
                sbCorrMatrix = np.zeros((sbs.shape[0], nantpol, nantpol), dtype=complex)
                for sbIdx, sb in enumerate(sbs):
                    if antGains is None:
                        sbCorrMatrix[sbIdx] = corrMatrix[sb, :, :] #select out a single subband, shape (nantpol, nantpol)
                    else: #Apply Gains
                        sbAntGains = antGains[sb][np.newaxis].T
                        sbVisGains = np.conjugate(np.dot(sbAntGains, sbAntGains.T)) # from Tobia, visibility gains are computed as (G . G^T)*
                        sbCorrMatrix[sbIdx] = np.multiply(sbVisGains, corrMatrix[sb, :, :]) #select out a single subband, shape (nantpol, nantpol)

                    #correct the time due to subband stepping
                    tOffset = (nchan - sb) * fDict['int'] #the time stamp in the filename is for the last subband
                    rem = tOffset - int(tOffset) #subsecond remainder
                    tDeltas.append(datetime.timedelta(0, int(tOffset), rem*1e6))

            elif fDict['fmt']=='xst':
                corrMatrix = np.fromfile(visFn, dtype='complex').reshape(1, nantpol, nantpol) #read in the correlation matrix
                if antGains is None:
                    sbCorrMatrix = corrMatrix #shape (nantpol, nantpol)
                else: #Apply Gains
                    sbAntGains = antGains[fDict['sb']][np.newaxis].T
                    sbVisGains = np.conjugate(np.dot(sbAntGains, sbAntGains.T)) # from Tobia, visibility gains are computed as (G . G^T)*
                    sbCorrMatrix = np.multiply(sbVisGains, corrMatrix) #shape (nantpol, nantpol)
                tDeltas = [datetime.timedelta(0, 0)] #no time offset

            print 'done'
            print 'CORRELATION MATRIX SHAPE', corrMatrix.shape
            
            obs = ephem.Observer() #create an observer at the array location
            obs.long = lon * (np.pi/180.)
            obs.lat = lat * (np.pi/180.)
            obs.elevation = float(elev)
            obs.epoch = fDict['ts']
            obs.date = fDict['ts']
            obsLat = float(obs.lat) #radians
            obsLong = float(obs.long) #radians
            print 'Observatory:', obs

            #get the UVW and visibilities for the different subbands
            ncorrs = nants*(nants+1)/2
            uvw = np.zeros((ncorrs, 3, len(sbs)), dtype=float)
            xxVis = np.zeros((ncorrs, len(sbs)), dtype=complex)
            yxVis = np.zeros((ncorrs, len(sbs)), dtype=complex)
            xyVis = np.zeros((ncorrs, len(sbs)), dtype=complex)
            yyVis = np.zeros((ncorrs, len(sbs)), dtype=complex)
            for sbIdx, sb in enumerate(sbs):
                obs.epoch = fDict['ts'] - tDeltas[sbIdx]
                obs.date = fDict['ts'] - tDeltas[sbIdx]

                #in order to accommodate multiple observations/subbands at different times/sidereal times all the positions need to be rotated relative to sidereal time 0
                LSTangle = obs.sidereal_time() #radians
                print 'LST:',  LSTangle
                rotAngle = float(LSTangle) - float(obs.long) #adjust LST to that of the Observatory longitutude to make the LST that at Greenwich
                #to be honest, the next two lines change the LST to make the images come out but i haven't worked out the coordinate transforms, so for now these work without justification
                rotAngle += np.pi
                rotAngle *= -1
                #Rotation matrix for antenna positions
                rotMatrix = np.array([[np.cos(rotAngle), -1.*np.sin(rotAngle), 0.],
                                      [np.sin(rotAngle), np.cos(rotAngle),     0.],
                                      [0.,               0.,                   1.]]) #rotate about the z-axis

                #get antenna positions in ITRF (x,y,z) format and compute the (u,v,w) coordinates referenced to sidereal time 0, this works only for zenith snapshot xyz->uvw conversion
                xyz = np.dot(ants[:,0,:], rotMatrix)

                repxyz = np.repeat(xyz, nants, axis=0).reshape((nants, nants, 3))
                uu = SWHT.util.vectorize(repxyz[:,:,0] - repxyz[:,:,0].T)
                vv = SWHT.util.vectorize(repxyz[:,:,1] - repxyz[:,:,1].T)
                ww = SWHT.util.vectorize(repxyz[:,:,2] - repxyz[:,:,2].T)
                uvw[:, :, sbIdx] = np.vstack((uu, vv, ww)).T

                #split up polarizations, vectorize the correlation matrix, and drop the lower triangle
                xxVis[:, sbIdx] = SWHT.util.vectorize(sbCorrMatrix[sbIdx, 0::2, 0::2])
                yxVis[:, sbIdx] = SWHT.util.vectorize(sbCorrMatrix[sbIdx, 0::2, 1::2])
                xyVis[:, sbIdx] = SWHT.util.vectorize(sbCorrMatrix[sbIdx, 1::2, 0::2])
                yyVis[:, sbIdx] = SWHT.util.vectorize(sbCorrMatrix[sbIdx, 1::2, 1::2])

            #add visibilities to previously processed files
            xxVisComb = np.concatenate((xxVisComb, xxVis))
            xyVisComb = np.concatenate((xyVisComb, xyVis))
            yxVisComb = np.concatenate((yxVisComb, yxVis))
            yyVisComb = np.concatenate((yyVisComb, yyVis))
            uvwComb = np.concatenate((uvwComb, uvw))

            #ax.scatter(uvw[:,0], uvw[:,1], uvw[:,2], c=clr[vid])

            ##uv coverage plot
            #plt.plot(uvw[:,0,:], uvw[:,1,:], '.')
            #plt.show()

        elif fDict['fmt']=='ms': #MS-based visibilities
            decomp = True

            fDict['sb'] = sbs

            MS = tbls.table(visFn, readonly=True)
            data_column = opts.column.upper()
            uvw = MS.col('UVW').getcol() # [vis id, (u,v,w)]
            vis = MS.col(data_column).getcol() #[vis id, freq id, stokes id]
            vis = vis[:,sbs,:] #select subbands
            MS.close()

            #lat/long/lst information
            ANTS = tbls.table(visFn + '/ANTENNA')
            positions = ANTS.col('POSITION').getcol()
            ant0Lat, ant0Long, ant0hgt = SWHT.ecef.ecef2geodetic(positions[0,0], positions[0,1], positions[0,2], degrees=False) #use the first antenna in the table to get the array lat/long
            ANTS.close()
            SRC = tbls.table(visFn + '/SOURCE')
            direction = SRC.col('DIRECTION').getcol()
            obsLat = direction[0,1]
            obsLong = ant0Long
            LSTangle = direction[0,0]
            SRC.close()

            #freq information, convert uvw coordinates
            SW = tbls.table(visFn + '/SPECTRAL_WINDOW')
            freqs = SW.col('CHAN_FREQ').getcol()[0, sbs] # [nchan]
            print 'SUBBANDS:', sbs, '(', freqs/1e6, 'MHz)'
            SW.close()

            #in order to accommodate multiple observations at different times/sidereal times all the positions need to be rotated relative to sidereal time 0
            print 'LST:',  LSTangle
            rotAngle = float(LSTangle) - obsLong #adjust LST to that of the Observatory longitutude to make the LST that at Greenwich
            #to be honest, the next two lines change the LST to make the images come out but i haven't worked out the coordinate transforms, so for now these work without justification
            rotAngle += np.pi
            rotAngle *= -1
            #Rotation matrix for antenna positions
            rotMatrix = np.array([[np.cos(rotAngle), -1.*np.sin(rotAngle), 0.],
                                  [np.sin(rotAngle), np.cos(rotAngle),     0.],
                                  [0.,               0.,                   1.]]) #rotate about the z-axis
            uvwRot = np.dot(uvw, rotMatrix).reshape(uvw.shape[0], uvw.shape[1], 1)
            uvwRotRepeat = np.repeat(uvwRot, len(sbs), axis=2)

            #split up polarizations
            xxVis = vis[:,:,0] 
            xyVis = vis[:,:,1]
            yxVis = vis[:,:,2]
            yyVis = vis[:,:,3]

            #add visibilities to previously processed files
            xxVisComb = np.concatenate((xxVisComb, xxVis))
            xyVisComb = np.concatenate((xyVisComb, xyVis))
            yxVisComb = np.concatenate((yxVisComb, yxVis))
            yyVisComb = np.concatenate((yyVisComb, yyVis))
            uvwComb = np.concatenate((uvwComb, uvwRotRepeat))

            ##uv coverage plot
            #plt.plot(uvw[:,0], uvw[:,1], '.')
            #plt.show()

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

    #compute the ideal l_max given the average solid angle angular resolution of an l-mode is Omega ~ 4pi / 2l steradian, and if the PSF is circular theta ~ pi / l radians
    blLen = np.sqrt(uvwComb[:,0,:]**2. + uvwComb[:,1,:]**2. + uvwComb[:,2,:]**2.) #compute the baseline lengths (in meters)
    maxBl = np.max(blLen) #maximum baseline length (in meters)
    meanWl = cc / np.mean(freqs) #mean observing wavelength
    maxRes = 1.22 * meanWl / maxBl
    print 'MAXIMUM RES: %f (radians) %f (deg)'%(maxRes, maxRes * (180. / np.pi))
    idealLmax = int(np.pi / maxRes)
    print 'SUGGESTED L_MAX: %i, %i (oversample 3), %i (oversample 5)'%(idealLmax, idealLmax*3, idealLmax*5)

    #from mpl_toolkits.mplot3d import Axes3D
    #fig = plt.figure()
    #ax = fig.add_subplot(111, projection='3d')
    #ax.scatter(uvwComb[:,0], uvwComb[:,1], uvwComb[:,2])
    #plt.show()
    #exit()

    #decompose the input visibilities into spherical harmonics visibility coefficeints
    if decomp:
        #remove auto-correlations
        print 'AUTO-CORRELATIONS:', opts.autos
        if not opts.autos:
            autoIdx = np.argwhere(uvwComb[:,0]**2. + uvwComb[:,1]**2. + uvwComb[:,2]**2. == 0.)
            xxVisComb[autoIdx] = 0.
            xyVisComb[autoIdx] = 0.
            yxVisComb[autoIdx] = 0.
            yyVisComb[autoIdx] = 0.

        #prepare for SWHT
        print 'Performing Spherical Wave Harmonic Transform'
        print 'LMAX:', opts.lmax
        #TODO: only doing total intensity right now
        iImgCoeffs = SWHT.swht.swhtImageCoeffs(xxVisComb+yyVisComb, uvwComb, freqs, lmax=opts.lmax, lmin=opts.lmin)

        #save image coefficients to file
        if opts.ocoeffs is None: outCoeffPklFn = 'tempCoeffs.pkl'
        else: outCoeffPklFn = opts.pkl
        SWHT.fileio.writeCoeffPkl(outCoeffPklFn, iImgCoeffs, [float(obsLong), float(obsLat)], float(LSTangle))

    #Imaging
    if opts.of is None:
        if opts.imageMode.startswith('heal'): outFn = 'tempImage.hpx'
        else: outFn = 'tempImage.pkl'
    else: outFn = opts.of

    if opts.imageMode.startswith('2'): #Make a 2D hemispheric image
        fov = opts.fov * (np.pi/180.) #Field of View in radians
        px = [opts.pixels, opts.pixels]
        res = fov/px[0] #pixel resolution
        print 'Generating 2D Hemisphere Image of size (%i, %i)'%(px[0], px[1])
        print 'Resolution(deg):', res*180./np.pi
        img = SWHT.swht.make2Dimage(iImgCoeffs, res, px, phs=[0., float(obsLat)]) #0 because the positions have already been rotated to the zenith RA of the first snapshot, if multiple snaphsots this needs to be reconsidered
        img = np.fliplr(img)
        plt.imshow(np.abs(img), interpolation='nearest')
        plt.colorbar()

        #save complex image to pickle file
        print 'Writing image to file %s ...'%outFn,
        SWHT.fileio.writeSWHTImgPkl(outFn, img, fDict, mode='2D')
        print 'done'

    elif opts.imageMode.startswith('3'): #Make a 3D equal stepped image
        print 'Generating 3D Image with %i steps in theta and %i steps in phi'%(opts.pixels, opts.pixels)
        img, phi, theta = SWHT.swht.make3Dimage(iImgCoeffs, dim=[opts.pixels, opts.pixels])
        img = np.abs(img)
        #X = np.cos(theta-(np.pi/2.)) * np.cos(phi)
        #Y = np.cos(theta-(np.pi/2.)) * np.sin(phi)
        #Z = np.sin(theta-(np.pi/2.))
        X, Y, Z = SWHT.util.sph2cart(theta, phi)

        #http://stackoverflow.com/questions/22175533/what-is-the-equivalent-of-matlabs-surfx-y-z-c-in-matplotlib
        from mpl_toolkits.mplot3d import Axes3D
        from matplotlib import cm
        from matplotlib.colors import Normalize
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        imin = img.min()
        imax = img.max()
        scalarMap = cm.ScalarMappable(norm=Normalize(vmin=imin, vmax=imax), cmap=cm.jet)
        scalarMap.set_array(img)
        C = scalarMap.to_rgba(img)
        surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors=C, antialiased=True)
        fig.colorbar(scalarMap)

        #save complex image to pickle file
        print 'Writing image to file %s ...'%outFn,
        SWHT.fileio.writeSWHTImgPkl(outFn, [img, phi, theta], fDict, mode='3D')
        print 'done'

    elif opts.imageMode.startswith('heal'): #plot healpix and save healpix file using the opts.pkl name
        print 'Generating HEALPix Image with %i NSIDE'%(opts.pixels)
        #use the healpy.alm2map function as it is much faster, there is a ~1% difference between the 2 functions, this is probably due to the inner workings of healpy
        #m = SWHT.swht.makeHEALPix(iImgCoeffs, nside=opts.pixels)
        m = hp.alm2map(SWHT.util.array2almVec(iImgCoeffs), opts.pixels)

        #save complex image to HEALPix file
        print 'Writing image to file %s ...'%outFn,
        hp.write_map(outFn, np.abs(m), coord='C')
        print 'done'
    
    elif opts.imageMode.startswith('coeff'): #plot the complex coefficients
        iImgCoeffs[0,0] = 0 #zero out DC offset component

        plt.subplot(221)
        plt.title('Real Components')
        plt.imshow(iImgCoeffs.real, interpolation='nearest')
        plt.colorbar()

        plt.subplot(222)
        plt.title('Imaginary Components')
        plt.imshow(iImgCoeffs.real, interpolation='nearest')
        plt.imshow(iImgCoeffs.imag, interpolation='nearest')
        plt.colorbar()

        plt.subplot(223)
        plt.title('Amplitude (dB)')
        plt.imshow(iImgCoeffs.real, interpolation='nearest')
        plt.imshow(10.*np.log10(np.abs(iImgCoeffs)), interpolation='nearest')
        plt.colorbar()

        plt.subplot(224)
        plt.title('Phase')
        plt.imshow(iImgCoeffs.real, interpolation='nearest')
        plt.imshow(np.angle(iImgCoeffs), interpolation='nearest')
        plt.colorbar()

        plt.subplot(233)
        coeffsFlat = []
        mms = []
        lls = []
        for ll in np.arange(iImgCoeffs.shape[0]):
            for mm in np.arange(-1*ll, ll+1):
                mms.append(mm)
                lls.append(ll)
                coeffsFlat.append(iImgCoeffs[ll,ll+mm])
        coeffsFlat = np.array(coeffsFlat)
        plt.ylabel('Amplitude (dB)')
        plt.xlabel('l')
        plt.plot(lls, 10.*np.log10(np.abs(coeffsFlat)), '.')

        plt.subplot(236)
        plt.ylabel('Amplitude (dB)')
        plt.xlabel('m')
        plt.plot(mms, 10.*np.log10(np.abs(coeffsFlat)), '.')

    if not (opts.savefig is None): plt.savefig(opts.savefig)
    if not opts.nodisplay:
        if opts.imageMode.startswith('heal'): hp.mollview(np.abs(m))
        plt.show()
    
