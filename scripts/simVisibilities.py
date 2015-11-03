#!/usr/bin/env python
"""
Simulate visibilties based on a LOFAR station or Measurement Set from a HEALPIX map or set of Spherical Harmonics coefficients
"""

import sys,os
import numpy as np
import healpy as hp
import datetime
import ephem
import SWHT

#TODO: include beam
#TODO: mask out dec min/max
#import matplotlib.pyplot as plt
#TODO: test sim, ACC, XST, MS

if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] -i HEALPIX_MAP/COEFF_PKL MS/ACC/XST FILES\n' \
        '    LOFAR XST, HEALPIX map/coefficient pickle, station\n' \
        '    LOFAR ACC, HEALPIX map/coefficient pickle, station, subbands (optional)\n' \
        '    HEALPIX map/coefficient pickle, station, subbands, rcumode, timestamp\n' \
        '    Measurement Set, HEALPIX map/coefficient pickle, subbands (optional), column, mode\n')
    o.set_description(__doc__)
    o.add_option('-i', '--imap', dest='imap', default=None,
        help='REQUIRED: Input HEALPIX map or spherical harmonics coefficient file')
    o.add_option('--station', dest='station', default=None,
        help = 'LOFAR ONLY: station name, e.g. SE607, if this is used then the ant_field option is not required, default: None')
    o.add_option('-F', '--ant_field', dest='ant_field', default=None,
        help = 'LOFAR ONLY: AntennaField.conf file for the LOFAR station of the ACC files, default: None')
    o.add_option('-D', '--deltas', dest='deltas', default=None,
        help = 'LOFAR ONLY: iHBADeltas.conf file, only required for HBA imaging, default: None')
    o.add_option('--override', dest='override', action='store_true',
        help = 'LOFAR XST ONLY: override filename metadata for RCU, integration length, and subband')
    o.add_option('-r', '--rcumode', dest='rcumode', default=3, type='int',
        help = 'LOFAR ONLY: Station RCU Mode for simulation, usually 3,5,6,7, for XST it will override filename metadata default: 3(LBA High)')
    o.add_option('-t', '--ts', dest='ts', default=None,
        help = 'LOFAR ONLY: UTC timestamp for simulation, if using multiple subbands this timestamp is for the highest subband and nt_time option is used, format YYYYMMDD_HHMMSS') 
    o.add_option('-I', '--int', dest='int_time', default=1., type='float',
        help = 'LOFAR ONLY: Integration time, used for accurate zenith pointing, for XST it will override filename metadata, default: 1 second')
    o.add_option('-l', '--lmax', dest='lmax', default=32, type='int',
        help = 'HEALPIX ONLY: Maximum l spherical harmonic quantal number, default: 32')
    o.add_option('-c', '--column', dest='column', default='CORRECTED_DATA', type='str',
        help = 'MS ONLY: select which data column write visibilities to, default: CORRECTED_DATA')
    o.add_option('-m', '--mode', dest='mode', default='replace', type='str',
        help = 'MS ONLY: mode to write simulated visibilities to MS: replace, add, subtract, default: replace')
    o.add_option('-s', '--subband', dest='subband', default=None,
        help = 'Select which subband(s) to image, for ACC and MS it will select, for multiple subbands use X,Y,Z and for range use X_Y notation, for XST it will override filename metadata, default: all')
    opts, args = o.parse_args(sys.argv[1:])

    if opts.imap is None:
        print "ERROR: no input HEALPIX map or coefficient file set with -i/--imap option"
        exit()
    elif opts.imap.endswith('.hpx'): #input HEALPIX map, decompose into spherical harmonic coefficients
        #Get HEALPIX map
        m = None
        w = None
        print 'Opening:', opts.imap
        hpMap = hp.read_map(opts.imap, field=None, h=True)
        if len(hpMap)==2: #no weight map
            m, hdr = hpMap
        elif len(hpMap)==3: #weight map
            m, w, hdr = hpMap

        if w is not None: m /= w #divide by the pixel weights
        print 'Map :: min=%f :: max=%f'%(np.nanmin(m), np.nanmax(m))

        #Convert HEALPIX map into Alm spherical harmonics coefficients
        print 'Generating Spherical Harmonic Coefficients from map...',
        alms = hp.sphtfunc.map2alm(m, lmax=opts.lmax, mmax=opts.lmax)
        blm = SWHT.util.almVec2array(alms, opts.lmax)
        print 'done'

    elif opts.imap.endswith('.pkl'): #SWHT coefficient pickle
        print 'Loading Image Coefficients file:', opts.imap
        coeffDict = SWHT.fileio.readCoeffPkl(opts.imap)
        blm = coeffDict['coeffs']

    #parse subbands
    sbs = np.array(SWHT.util.convert_arg_range(opts.subband))
    
    #if simulating LOFAR visibilities, get the station information
    lofarStation = None
    if (not (opts.station is None)) or (not (opts.ant_field is None)): #If using LOFAR data, get station information
        lofarStation = SWHT.lofarConfig.getLofarStation(name=opts.station, affn=opts.ant_field, aafn=None, deltas=opts.deltas)

    #TODO: test
    if args==[]: #if not using measurement set or LOFAR data set as a template then try to do a simulation from the input options
        print 'No input datasets as a template, assuming LOFAR simulation with given input options'

        if opts.ts is None:
            print 'Error: missing a timestamp'
            exit()
        if lofarStation is None:
            print 'Error: missing a station'
            exit()
        if opts.subband is None: #use all subbands for LOFAR, 512 channels
            sbs = np.array(SWHT.util.convert_arg_range('0_511'))

        dd,tt = opts.ts.split('_')
        ts = datetime.datetime(year=int(dd[:4]), month=int(dd[4:6]), day=int(dd[6:]), hour=int(tt[:2]), minute=int(tt[2:4]), second=int(tt[4:]))
        rcumode = opts.rcumode

        print 'Using Station: %s, RCU Mode: %i, Time: %s, Subbands: '%(lofarStation.name, opts.rcumode, opts.ts), sbs

        #longitude and latitude of array
        arr_xyz = lofarStation.antField.location[SWHT.lofarConfig.rcuInfo[rcumode]['array_type']]
        lat, lon, elev = SWHT.ecef.ecef2geodetic(arr_xyz[0], arr_xyz[1], arr_xyz[2], degrees=True)
        print 'LON(deg):', lon, 'LAT(deg):', lat, 'ELEV(m):', elev

        #antenna positions
        ants = lofarStation.antField.antpos[SWHT.lofarConfig.rcuInfo[rcumode]['array_type']]
        nants = ants.shape[0]
        print 'NANTENNAS:', nants
        npols = 2
        nantpol = nants * npols

        #frequency information
        nchan = SWHT.lofarConfig.rcuInfo[rcumode]['nchan']
        bw = SWHT.lofarConfig.rcuInfo[rcumode]['bw']
        df = bw/nchan
        freqs = sbs*df + SWHT.lofarConfig.rcuInfo[rcumode]['offset'] + (df/2.) #df/2 to centre the band
        print 'SUBBANDS:', sbs, '(', freqs/1e6, 'MHz)'
        npols = 2

        obs = ephem.Observer() #create an observer at the array location
        obs.long = lon * (np.pi/180.)
        obs.lat = lat * (np.pi/180.)
        obs.elevation = float(elev)
        obs.epoch = ts
        obs.date = ts
        obsLat = float(obs.lat) #radians
        obsLong = float(obs.long) #radians
        print 'Observatory:', obs

        #get the UVW and visibilities for the different subbands
        ncorrs = nants*(nants+1)/2
        maxChan = np.max(sbs)
        for sbIdx, sb in enumerate(sbs):
            #correct the time due to subband stepping
            tOffset = (maxChan - sb) * opts.int_time
            rem = tOffset - int(tOffset) #subsecond remainder
            tDelta = datetime.timedelta(0, int(tOffset), rem*1e6)
            obs.epoch = ts - tDelta
            obs.date = ts - tDelta

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
            uvw = np.vstack((uu, vv, ww)).T

            #Compute visibilities from brightness coefficients
            vis = SWHT.swht.iswhtVisibilities(blm, uvw, np.array([freqs[sbIdx]]))

            #remove auto-correlations
            #autoIdx = np.argwhere(uvw[:,0]**2. + uvw[:,1]**2. + uvw[:,2]**2. == 0.)
            #vis[autoIdx] = 0.
            iImgCoeffs = SWHT.swht.swhtImageCoeffs(vis, uvw, np.array([freqs[sbIdx]]), lmax=blm.shape[0]-1)
            
            from matplotlib import pyplot as plt
            plt.subplot(131)
            plt.imshow(10.*np.log10(np.abs(blm)), interpolation='nearest')
            plt.colorbar()
            plt.subplot(132)
            plt.imshow(10.*np.log10(np.abs(iImgCoeffs)), interpolation='nearest')
            plt.colorbar()
            plt.subplot(133)
            plt.imshow(10.*np.log10(np.abs(iImgCoeffs-blm)), interpolation='nearest')
            plt.colorbar()
            plt.show()

            SWHT.fileio.writeCoeffPkl('reverseTestCoeffs.pkl', iImgCoeffs, [0., 0.], 0.)
            exit()

            #Build a correlation matrix for a single polarization
            corrMatrix = np.zeros((nants, nants), dtype=complex)
            triu_indices = np.triu_indices(nants)
            tril_indices = np.tril_indices(nants)
            corrMatrix[tril_indices[0], tril_indices[1]] = np.conjugate(vis[:,0]/2.) #lower half is redundant, a conjugate of the upper half, 1/2 power for one polarization
            corrMatrix[triu_indices[0], triu_indices[1]] = vis[:,0]/2.

            #Create a LOFAR XST correlation matrix, XY and YX cross-pols are set to 0.
            fullCorrMatrix = np.zeros((nantpol, nantpol), dtype=complex)
            fullCorrMatrix[::2, ::2] = corrMatrix #XX
            fullCorrMatrix[1::2, 1::2] = corrMatrix #YY

            #Save simulated visibilities to XST file
            #20150915_191137_rcu5_sb60_int10_dur10_xst.dat
            xst_fn = dd + '_' + tt + '_rcu%i'%rcumode + '_sb%i'%sb + '_int%i'%int(opts.int_time) + '_dur%i'%int(opts.int_time) + '_xst.dat.sim'
            print 'Saving simulated visibilities to', xst_fn
            fullCorrMatrix.tofile(xst_fn)

    else: #using input files as templates
        visFiles = args
        for vid,visFn in enumerate(visFiles):
            print 'Simulating visibilities for %s (%i/%i)'%(visFn, vid+1, len(visFiles))
            fDict = SWHT.fileio.parse(visFn)

            #TODO: test
            if fDict['fmt']=='acc' or fDict['fmt']=='xst': #LOFAR visibilities
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
                for sbIdx, sb in enumerate(sbs):
                    #correct the time due to subband stepping
                    tOffset = (nchan - sb) * fDict['int'] #the time stamp in the filename is for the last subband
                    rem = tOffset - int(tOffset) #subsecond remainder
                    tDelta = datetime.timedelta(0, int(tOffset), rem*1e6)
                    obs.epoch = fDict['ts'] - tDelta
                    obs.date = fDict['ts'] - tDelta

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
                    uvw = np.vstack((uu, vv, ww)).T

                    #Compute visibilities from brightness coefficients
                    vis = SWHT.swht.iswhtVisibilities(blm, uvw, np.array([freqs[sbIdx]]))

                    #remove auto-correlations
                    #autoIdx = np.argwhere(uvw[:,0]**2. + uvw[:,1]**2. + uvw[:,2]**2. == 0.)
                    #vis[autoIdx] = 0.
                    iImgCoeffs = SWHT.swht.swhtImageCoeffs(vis, uvw, np.array([freqs[sbIdx]]), lmax=blm.shape[0]-1)
                    
                    from matplotlib import pyplot as plt
                    plt.subplot(131)
                    plt.imshow(10.*np.log10(np.abs(blm)), interpolation='nearest')
                    plt.colorbar()
                    plt.subplot(132)
                    plt.imshow(10.*np.log10(np.abs(iImgCoeffs)), interpolation='nearest')
                    plt.colorbar()
                    plt.subplot(133)
                    plt.imshow(10.*np.log10(np.abs(iImgCoeffs-blm)), interpolation='nearest')
                    plt.colorbar()
                    plt.show()

                    SWHT.fileio.writeCoeffPkl('reverseTestCoeffs.pkl', iImgCoeffs, [0., 0.], 0.)
                    exit()

                    #Build a correlation matrix for a single polarization
                    corrMatrix = np.zeros((nants, nants), dtype=complex)
                    triu_indices = np.triu_indices(nants)
                    tril_indices = np.tril_indices(nants)
                    corrMatrix[tril_indices[0], tril_indices[1]] = np.conjugate(vis[:,0]/2.) #lower half is redundant, a conjugate of the upper half, 1/2 power for one polarization
                    corrMatrix[triu_indices[0], triu_indices[1]] = vis[:,0]/2.

                    #Create a LOFAR XST correlation matrix, XY and YX cross-pols are set to 0.
                    fullCorrMatrix = np.zeros((nantpol, nantpol), dtype=complex)
                    fullCorrMatrix[::2, ::2] = corrMatrix #XX
                    fullCorrMatrix[1::2, 1::2] = corrMatrix #YY

                    #Save simulated visibilities to XST file
                    #20150915_191137_rcu5_sb60_int10_dur10_xst.dat
                    xst_fn = dd + '_' + tt + '_rcu%i'%rcumode + '_sb%i'%sb + '_int%i'%int(opts.int_time) + '_dur%i'%int(opts.int_time) + '_xst.dat.sim'
                    print 'Saving simulated visibilities to', xst_fn
                    fullCorrMatrix.tofile(xst_fn)

            #TODO: test
            elif fDict['fmt']=='ms': #MS-based visibilities
                print visFn
                #TODO
                #open MS, get XYZ positions, convert to r,theta,phi
                #Compute visibilities from brightness coefficients
                #vis = iswhtVisibilities(blm, uvw, freqs)
                #save visibilities to column

            else:
                print 'ERROR: unknown data format, exiting'
                exit()

