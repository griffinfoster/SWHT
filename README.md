SWHT
===
  
Contact: griffin.foster@gmail.com  

A python package for generating radio interferometry images from LOFAR station ACC and XST files, and from widefield, low-frequency measurement sets (e.g. PAPER) using a Spherical Wave Harmonic Transform ([Imaging on a Sphere with Interferometers: the Spherical Wave Harmonic Transform](http://arxiv.org/abs/1504.04485)) and a standard 2D Fourier Transform.

#### Required Python Modules

* matplotlib
* numpy
* scipy [special functions]
* [ephem](http://rhodesmill.org/pyephem/) [observatories]

#### Optional Python Modules

* [healpy](https://healpy.readthedocs.org/en/latest/) [HEALPIX interface]
* [python-casacore](https://github.com/casacore/python-casacore)

#### Install

To install the current stable version (0.1.2) use pip:

```
pip install SWHT
```

While developing it is useful to do a developer install:

```
sudo python setup.py develop
```

Otherwise, the standard install will install the package:

```
sudo python setup.py install  
```

#### Scripts

* ftVisibilities.py: 2D Fourier Transform of LOFAR ACC, XST files and Measurement Sets  
* gsm2healpix.py: convert the output of the [GSM (Global Sky Model)](http://space.mit.edu/~angelica/gsm/index.html) to a HEALPIX map
* imageSWHTcoeffs.py: generate images and HEALPIX maps from pre-computed SWHT image coefficients  
* plotHealpix.py: general HEALPIX plotting script
* simVisibilities.py: simulate visibilities from Spherical Harmonics coefficients or a HEALPIX map
* swhtVisibilities.py: Spherical Wave Harmonic Transform of LOFAR ACC, XST files and Measurement Sets  

#### Examples

Example LBA and HBA correlation files is available at:  

* LBA (rcumode 1) :
	* UK608: https://zenodo.org/record/840405/files/20150607_122433_acc_512x192x192.dat
	* SE607: https://zenodo.org/record/840405/files/20120513_052251_acc_512x192x192.dat
* HBA (rcumode 5) :
	* SE607:  https://zenodo.org/record/840405/files/20150915_191137_rcu5_sb60_int10_dur10_elf0f39fe2034ea85fc02b3cc1544863053b328fd83291e880cd0bf3c3d3a50a164a3f3e0c070c73d073f4e43849c0e93b_xst.dat

For any script, use the '-h' argument to print out help on available input options.

```
ftVisibilities.py ../examples/20150915_191137_rcu5_sb60_int10_dur10_elf0f39fe2034ea85fc02b3cc1544863053b328fd83291e880cd0bf3c3d3a50a164a3f3e0c070c73d073f4e43849c0e93b_xst.dat --station=SE607 -p 128 --conv=prolate --autos
ftVisibilities.py ../examples/20150607_122433_acc_512x192x192.dat -s 300 --station=SE607 --conv=gauss -p 64
ftVisibilities.py ../examples/zen.2455819.69771.uvcRREM.MS -s 40 --conv=fast -p 256

swhtVisibilities.py --station=UK608 ../examples/20120513_052251_acc_512x192x192.dat -s 299 -l 24
swhtVisibilities.py --station=SE607 ../examples/20150607_122433_acc_512x192x192.dat -s 299 -l 24
swhtVisibilities.py --station=SE607 ../examples/20150915_191137_rcu5_sb60_int10_dur10_elf0f39fe2034ea85fc02b3cc1544863053b328fd83291e880cd0bf3c3d3a50a164a3f3e0c070c73d073f4e43849c0e93b_xst.dat -s 100 -l 32
swhtVisibilities.py --station=UK608 tempCoeffs.pkl -I coeff
```
