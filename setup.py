from setuptools import setup, find_packages
import numpy as np
import os, sys, glob

__version__ = '0.1.1' #this needs to be kept up to date with SWHT/__init__.py

setup(name = 'SWHT',
    version = __version__,
    description = 'Spherical Wave Harmonic Transform for radio astronomy visibility data',
    long_description = "Spherical Wave Harmonic Transform for radio astronomy visibility data as described by T. Carozzi in 'Imaging on a Sphere with Interferometers: the Spherical Wave Harmonic Transform'",
    author = 'Griffin Foster',
    author_email = 'griffin.foster@gmail.com',
    url='https://github.com/griffinfoster/SWHT',
    platforms = ['*nix'],
    license = 'GPL',
    requires = ['distutils','numpy','scipy','matplotlib','healpy','ephem'],
    provides = ['SWHT'],
    packages = ['SWHT'],
    include_dirs = [np.get_include()],
    package_data = {'SWHT': ['data/LOFAR/StaticMetaData/*.conf', 'data/LOFAR/*.txt']},
    #scripts = glob.glob('scripts/*.py'),
    scripts = ['scripts/ftVisibilities.py', 'scripts/gsm2healpix.py', 'scripts/imageSWHTcoeffs.py', 'scripts/plotHealpix.py', 'scripts/simVisibilities.py', 'scripts/swhtVisibilities.py'],
    classifiers = [
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Natural Language :: English',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 2.7',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Topic :: Scientific/Engineering :: Astronomy',
    ],
)

