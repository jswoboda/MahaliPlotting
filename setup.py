#!/usr/bin/env python
from setuptools import setup

# apt install libgeos-dev libgeos++-dev
# Note basemap 1.1 should be back on Pypi.
req = ['nose','numpy','matplotlib','seaborn','pytables','h5py','xarray','pandas', 'astropy','scipy','pyqt']
pipreq=['basemap']

try:
    import conda.cli
    conda.cli.main('install',*req)
except Exception as e:
    import pip
    pip.main(['install'] + req)
#%% install
from setuptools import setup

setup(name='MahaliPlotting',
      packages=['MahaliPlotting'],
      url='https://github.com/jswoboda/MahaliPlotting',
      classifiers=[
      'Intended Audience :: Science/Research',
      'Development Status :: 3 - Alpha',
      'Topic :: Scientific/Engineering :: Atmospheric Science',
      'Programming Language :: Python',
      ],
      install_requires=req+pipreq,
      dependency_links=['https://downloads.sourceforge.net/project/matplotlib/matplotlib-toolkits/basemap-1.0.7/basemap-1.0.7.tar.gz'],
	  )
