"""
Functions and classes to read and parse LOFAR station configuration files
"""

import numpy as np
import glob
import os
import struct

#HARDCODED path, need to make this configurable to the module, something in setup.py?
#StaticMetaData = '/home/griffin/sw/SWHT/SWHT/data/LOFAR/StaticMetaData/'
StaticMetaData = __file__.split('lofarConfig.py')[0] + 'data/LOFAR/StaticMetaData/'

rcuInfo = [ {'mode':'OFF', 'rcuID':0, 'array_type':'LBA', 'bw':100000000., 'offset':0., 'nchan':512},                    #0
            {'mode':'LBL_HPF10MHZ', 'rcuID':1, 'array_type':'LBA', 'bw':100000000., 'offset':0., 'nchan':512},           #1
            {'mode':'LBL_HPF30MHZ', 'rcuID':2, 'array_type':'LBA', 'bw':100000000., 'offset':0., 'nchan':512},           #2
            {'mode':'LBH_HPF10MHZ', 'rcuID':3, 'array_type':'LBA', 'bw':100000000., 'offset':0., 'nchan':512},           #3
            {'mdde':'LBH_HPF30MHZ', 'rcuID':4, 'array_type':'LBA', 'bw':100000000., 'offset':0., 'nchan':512},           #4
            {'mode':'HBA_110_190MHZ', 'rcuID':5, 'array_type':'HBA', 'bw':100000000., 'offset':100000000., 'nchan':512}, #5
            {'mode':'HBA_170_230MHZ', 'rcuID':6, 'array_type':'HBA', 'bw': 80000000., 'offset':160000000., 'nchan':512}, #6
            {'mode':'HBA_210_290MHZ', 'rcuID':7, 'array_type':'HBA', 'bw':100000000., 'offset':200000000., 'nchan':512}] #7

def getLofarStation(name=None, affn=None, aafn=None, deltas=None, noarrays=True):
    """Create an instance of the lofarStation class based on a station name and the configuration repo, or from the antenna array and field files
    name: station name (required if no filenames used)
    affn, aafn: AntennaArray filename, AntennaField filename (required if station name is not used)
    deltas: HBA tile deltas filename, optional, only used in HBA imaging
    noarrays: do not require a *-AntennaArrays.conf file
    """
    nameValid = False #used to check if the input station name is valid
    confValid = False #used to check if the input ant_array and ant_field pair is valid
    dfn = None #HBA deltas filename to pass on, favour input filename over repository file
    if name is not None: #check if the station name is valid
        repoaafn = glob.glob(StaticMetaData+name+'-AntennaArrays.conf')
        repoaffn = glob.glob(StaticMetaData+name+'-AntennaField.conf')
        repodfn = glob.glob(StaticMetaData+name+'-iHBADeltas.conf')
        if len(repoaafn)==1 and noarrays: #case: only require AntennaField.conf file
            repoaffn = repoaffn[0]
            repoaafn = None
            nameValid = True
        elif len(repoaafn)==1 and len(repoaffn)==1: #case: using bother AntennaArrays.conf and AntennaField.conf files
            repoaafn = repoaafn[0]
            repoaffn = repoaffn[0]
            nameValid = True
        if len(repodfn)==1: dfn = repodfn[0] #found a valid HBADeltas file in the repo

    if affn is not None and aafn is not None: #check if ant_array and ant_field is valid
        confValid = True

    if nameValid is False and confValid is False: #no valid information for forming station
        print 'ERROR: input station name or configuration files are no good, exiting'
        exit()

    #favour the input config files over the station name and repo config files
    if deltas is not None: dfn = deltas
    if confValid:
        stationName = affn.split('/')[-1].split('-')[0]
        print 'Station Name:', stationName
        print 'AntennaArray:', aafn
        print 'AntennaField:', affn
        print 'iHBADeltas:', dfn
        return lofarStation(stationName, affn, aafn, deltas=dfn)
    elif nameValid:
        print 'Station Name:', name
        print 'AntennaArray:', repoaafn
        print 'AntennaField:', repoaffn
        print 'iHBADeltas:', dfn
        return lofarStation(name, repoaffn, repoaafn, deltas=dfn)

class lofarStation():
    def __init__(self, name, affn, aafn=None, deltas=None):
        """deltas: optional HBA deltas file
        """
        self.name = name
        self.antField = antennaField(name, affn)
        if aafn is not None: self.antArrays = antennaArrays(name, aafn)

        if deltas is not None: self.deltas = getHBADeltas(deltas)
        else: self.deltas = None
        
        if name.lower().startswith('cs'): self.stype = 'core'
        elif name.lower().startswith('rs'): self.stype = 'remote'
        else: self.stype = 'international'

def getHBADeltas(fn):
    """Interface to a 16x3 delta offset XYZ  position for the elements of an HBA tile relative to the position in antennaField
    """
    fh=open(fn)
    lines = []
    lines = fh.read().split('\n')
    fh.close()
    cleanStr = ''
    for l in lines:
        if l=='' or l.startswith('#'): continue
        cleanStr += " ".join(l.split())+" "
    cleanStr = " ".join(cleanStr.split())
    return np.array(map(float, cleanStr.split(' ')[5:-1])).reshape((16, 3))

class antennaField():
    def __init__(self, name, fn):
        self.name = name
        self.rotMatrix = {} #roation matrix to transform local horizon positions to relative ITRF positions, needs to be inverted to convert relative ITRF positions to local horizon positions
        self.normVec = {}
        self.antpos = {} #ITRF positions relative to station ITRF in self.location
        self.localAntPos = {} #local horizon plane positions obtained from applying station rotation matrix to self.antpos
        self.location = {} #station ITRF position
        fh=open(fn)
        lines = []
        lines = fh.read().split('\n')
        fh.close()
        cleanStr = ''
        for l in lines:
            if l=='' or l.startswith('#'): continue
            cleanStr += " ".join(l.split())+" "
        cleanStr = " ".join(cleanStr.split())
        lastMode = None
        for l in cleanStr.split(']'):
            if l=='': continue
            infoStr,dataStr = l.split('[')
            infoStr = infoStr.strip()
            if infoStr.lower().startswith('normal'):
                name,mode,length = infoStr.split(' ')
                self.normVec[mode] = np.array(map(float, dataStr.strip().split(' ')))
            elif infoStr.lower().startswith('rotation'):
                name,mode,l0,fill,l1 = infoStr.split(' ')
                self.rotMatrix[mode] = np.array(map(float, dataStr.strip().split(' '))).reshape((int(l0), int(l1)))
            elif infoStr.lower().startswith('lba') or infoStr.lower().startswith('hba'):
                mode,length = infoStr.split(' ')
                self.location[mode] = np.array(map(float, dataStr.strip().split(' ')))
                lastMode = mode
            else: #antenna positions
                l0,f0,l1,f1,l2 = infoStr.split(' ')
                self.antpos[lastMode] = np.array(map(float, dataStr.strip().split(' '))).reshape((int(l0), int(l1), int(l2)))
        #convert antenna positions to local horizon coordinate system
        for mode in self.antpos:
            #a bit hacky, but some stations use HBA0 and HBA1 for the rotation matrix and HBA for the antenna postions
            for rkey in self.rotMatrix:
                if mode.startswith(rkey):
                    rotMode = mode
                    continue
            self.localAntPos[mode] = np.zeros_like(self.antpos[mode])
            self.localAntPos[mode][:,0,:] = np.linalg.lstsq(self.rotMatrix[rotMode], self.antpos[mode][:,0,:].T)[0].T
            self.localAntPos[mode][:,1,:] = np.linalg.lstsq(self.rotMatrix[rotMode], self.antpos[mode][:,1,:].T)[0].T

class antennaArrays():
    def __init__(self, name, fn):
        """Parse the AntenneArrays file, most of the informationis redundant to the AntennaField file, but contains the (lat, long, height)
        DEPRECIATED: these files were used to get the station (lat, lon, h) only, but that is now computed with the array X,Y,Z positions in antennaField() using ecef.py
        """
        self.name = name
        self.antpos = {}
        self.location = {}
        fh = open(fn)
        lines = []
        lines = fh.read().split('\n')
        fh.close()
        cleanStr = ''
        for l in lines:
            if l=='' or l.startswith('#'): continue
            cleanStr += " ".join(l.split())+" "
        cleanStr = " ".join(cleanStr.split())
        lastMode = None
        for l in cleanStr.split(']'):
            if l=='': continue
            infoStr,dataStr = l.split('[')
            infoStr = infoStr.strip()
            if infoStr.lower().startswith('lba') or infoStr.lower().startswith('hba'):
                mode,length = infoStr.split(' ')
                self.location[mode] = np.array(map(float, dataStr.strip().split(' ')))
                lastMode = mode
            else: #antenna positions
                l0,f0,l1,f1,l2 = infoStr.split(' ')
                self.antpos[lastMode] = np.array(map(float, dataStr.strip().split(' '))).reshape((int(l0), int(l1), int(l2)))
        #print self.antpos.keys(),self.location.keys()

def rotationMatrix(alpha, beta, gamma):
    """Generic rotation matrix to apply to an XYZ co-ordinate"""
    return np.matrix([[np.cos(beta)*np.cos(gamma), np.cos(gamma)*np.sin(alpha)*np.sin(beta)-np.cos(alpha)*np.sin(gamma), np.cos(alpha)*np.cos(gamma)*np.sin(beta)+np.sin(alpha)*np.sin(gamma)],
                      [np.cos(beta)*np.sin(gamma), np.cos(alpha)*np.cos(gamma)+np.sin(alpha)*np.sin(beta)*np.sin(gamma), -1.*np.cos(gamma)*np.sin(alpha)+np.cos(alpha)*np.sin(beta)*np.sin(gamma)],
                      [-1.*np.sin(beta),           np.cos(beta)*np.sin(alpha),                                           np.cos(alpha)*np.cos(beta)]])

def rotMatrixfromXYZ(station, mode='LBA'):
    """Return a rotation matrix which will rotate a station to (0,0,1)"""
    loc = station.antField.location[mode]
    longRotMat = rotationMatrix(0., 0., -1.*np.arctan(loc[1]/loc[0]))
    loc0 = np.dot(longRotMat, loc)
    latRotMat = rotationMatrix(0., np.arctan(loc0[0,2]/loc0[0,0]), 0.)
    return np.dot(latRotMat, longRotMat)

def applyRotMatrix(station, rm, mode='LBA'):
    """Apply a rotation matrix to a station location, returns new location after apply rotation matrix"""
    loc = station.antField.location[mode]
    return np.dot(rm,loc)

def relativeStationOffset(s0, s1, mode='LBA'):
    """Return the relative offset of station s1 from station s0"""
    rotMat = rotMatrixfromXYZ(s0, 'LBA')
    s0loc = applyRotMatrix(s0, rotMat, mode)
    s1loc = applyRotMatrix(s1, rotMat, mode)
    return np.array(s1loc-s0loc)[0][::-1]

def readCalTable(fn, nants=96, nsbs=512, npols=2):
    """Parse a LOFAR Calibration Table
    return: [NSBS, NANTS * NPOLS] complex array, the X and Y pols are interlaced, not serial, i.e. xGains = antGains[:, 0::2] and yGains = antGains[:, 1::2]
    """
    fh = open(fn)
    # Test for header record above raw data - present in newer caltables (starting 2012)
    line = fh.readline()
    if 'HeaderStart' in line:
        while not 'HeaderStop' in line:
            line = fh.readline()
    else:  # no header present, seek to starting position
        file.seek(0)

    fmt = str(nants * npols * nsbs * 2) + 'd'
    sz = struct.calcsize(fmt)
    antGains = np.array( struct.unpack(fmt, fh.read(sz)) )
    fh.close()

    antGains = antGains[0::2] + 1j * antGains[1::2] #make the array complex
    antGains.resize(nsbs, nants * npols)

    return antGains

if __name__ == '__main__':
    print 'Running test cases'

    print 'Using data from this directory:', StaticMetaData

    deltas = getHBADeltas('data/LOFAR/StaticMetaData/SE607-iHBADeltas.conf')
    print deltas

    #antfield=antennaField('CS013','data/LOFAR/StaticMetaData/CS013-AntennaField.conf')
    #antfield=antennaField('RS208','data/LOFAR/StaticMetaData/RS208-AntennaField.conf')
    #antfield=antennaField('UK608','data/LOFAR/StaticMetaData/UK608-AntennaField.conf')

    #antArrys=antennaArrays('CS013','data/LOFAR/StaticMetaData/CS013-AntennaArrays.conf')
    #antArrys=antennaArrays('RS208','data/LOFAR/StaticMetaData/RS208-AntennaArrays.conf')
    #antArrys=antennaArrays('UK608','data/LOFAR/StaticMetaData/UK608-AntennaArrays.conf')
    CS013 = lofarStation('CS013','data/LOFAR/StaticMetaData/CS013-AntennaField.conf', 'data/LOFAR/StaticMetaData/CS013-AntennaArrays.conf')
    print CS013.name
    RS208 = lofarStation('RS208','data/LOFAR/StaticMetaData/RS208-AntennaField.conf', 'data/LOFAR/StaticMetaData/RS208-AntennaArrays.conf')
    print RS208.name
    UK608 = lofarStation('UK608','data/LOFAR/StaticMetaData/UK608-AntennaField.conf', 'data/LOFAR/StaticMetaData/UK608-AntennaArrays.conf')
    print UK608.name
    SE607 = lofarStation('SE607','data/LOFAR/StaticMetaData/SE607-AntennaField.conf')
    print SE607.name

    CS002 = lofarStation('CS002','data/LOFAR/StaticMetaData/CS002-AntennaField.conf', 'data/LOFAR/StaticMetaData/CS002-AntennaArrays.conf')
    CS003 = lofarStation('CS003','data/LOFAR/StaticMetaData/CS003-AntennaField.conf', 'data/LOFAR/StaticMetaData/CS003-AntennaArrays.conf')
    rotMat = rotMatrixfromXYZ(CS002,'LBA')
    #print applyRotMatrix(CS002,rotMat,'LBA')
    #print applyRotMatrix(CS003,rotMat,'LBA')
    #print applyRotMatrix(CS013,rotMat,'LBA')
    #print applyRotMatrix(RS208,rotMat,'LBA')
    print relativeStationOffset(CS002,CS003)
    print relativeStationOffset(CS002,CS013)
    print relativeStationOffset(CS002,RS208)
    print relativeStationOffset(CS002,UK608)

    #print '\n'
    #aa=CS002.antField.location['LBA']
    #print aa
    #print n.sqrt(aa[0]**2.+aa[1]**2.+aa[2]**2.)
    #print CS002.antField.rotMatrix['LBA']
    #ab=n.dot(n.linalg.inv(CS002.antField.rotMatrix['LBA']),CS002.antField.location['LBA'])
    #print ab
    #print n.sqrt(ab[0]**2.+ab[1]**2.+ab[2]**2.)

    getLofarStation(name='SE607')
    getLofarStation(affn='data/LOFAR/StaticMetaData/SE607-AntennaField.conf', aafn='data/LOFAR/StaticMetaData/SE607-AntennaArrays.conf')
    getLofarStation(affn='data/LOFAR/StaticMetaData/SE607-AntennaField.conf', aafn='data/LOFAR/StaticMetaData/SE607-AntennaArrays.conf', deltas='data/LOFAR/StaticMetaData/SE607-iHBADeltas.conf')

    print 'Made it through without any errors.'

