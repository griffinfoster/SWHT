### LOFAR Station Configuration Files

This is a copy of the LOFAR station configuration files directory from the UK608 station (from 2013 or so, it would be good to get updated files).

* HBAtile.txt : the local tangent plane position (in meters) of the 16 HBA elements in a tile. I can not confirm the x,y,z reference directions.
* *-AntennaField.conf: station position in ITRF X,Y,Z; antenna positions in relative ITRF X,Y,Z coordinates
* *-AntennaArrays.conf: station (lat, lon, h) and antenna positions in local tangent plane positions (depreciated)
* *-iHBADeltas.conf: relative ITRF X,Y,Z coordinates of the 16 HBA elements
* *-HBADeltas.conf: not sure what coordinate system these are in

#### Of note:

* AntennaArrays.conf files are largely depreciated, only AntennaField.conf files are used
* SE607-AntennaArrays.conf is a copy of UK608-AntennaArrays.conf
* DE602-AntennaArrays.conf, DE603-AntennaArrays.conf and DE605-AntennaArrays.conf are the same
* There is no FI609-AntennaArrays.conf file
