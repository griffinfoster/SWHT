### LOFAR Station Configuration Files

Most files taken from ASTRON SVN (Revision 38931): https://svn.astron.nl/LOFAR/trunk/MAC/Deployment/data/StaticMetaData/

* HBAtile.txt : the local tangent plane position (in meters) of the 16 HBA elements in a tile. I can not confirm the x,y,z reference directions.
* *-AntennaField.conf: station position in ITRF X,Y,Z; antenna positions in relative ITRF X,Y,Z coordinates
* *-AntennaArrays.conf: station (lat, lon, h) and antenna positions in local tangent plane positions (depreciated)
* *-iHBADeltas.conf: relative ITRF X,Y,Z coordinates of the 16 HBA elements
* *-HBADeltas.conf: not sure what coordinate system these are in

