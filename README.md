SWHT
===

Created: 30.09.15 
Last Modified: 09.10.15  
Contact: griffin.foster@gmail.com  

A python package for generating radio interferometry images from LOFAR station ACC and XST files, and from widefield, low-frequency measurement sets (e.g. PAPER) using a standard Fourier Transform and with a Spherical Wave Harmonic Transform ([Imaging on a Sphere with Interferometers: the Spherical Wave Harmonic Transform](http://arxiv.org/abs/1504.04485)).

#### Required Python Modules

* matplotlib 
* numpy 
* pyrap 
* ephem 

#### Install

While developing it is useful to do a developer install:

```
sudo python setup.py develop
```

Otherwise, the standard install will work:

```
sudo python setup.py install  
```

#### Examples

```
./ftVisibilities.py ../examples/20150915_191137_rcu5_sb60_int10_dur10_elf0f39fe2034ea85fc02b3cc1544863053b328fd83291e880cd0bf3c3d3a50a164a3f3e0c070c73d073f4e43849c0e93b_xst.dat --station=SE607 -p 128 --conv=prolate --autos
./ftVisibilities.py ../examples/20150607_122433_acc_512x192x192.dat -s 300 --station=SE607 --conv=gauss -p 64
./ftVisibilities.py ../examples/zen.2455819.69771.uvcRREM.MS -s 40 --conv=fast -p 256
./swhtVisibilities.py --station=UK608 ../examples/20120513_052251_acc_512x192x192.dat -s 299 -l 24 --autos
./swhtVisibilities.py --station=SE607 ../examples/20150607_122433_acc_512x192x192.dat -s 299 -l 24 --autos
./swhtVisibilities.py --station=SE607 ../examples/20150915_191137_rcu5_sb60_int10_dur10_elf0f39fe2034ea85fc02b3cc1544863053b328fd83291e880cd0bf3c3d3a50a164a3f3e0c070c73d073f4e43849c0e93b_xst.dat -s 100 -l 32 --autos

```

