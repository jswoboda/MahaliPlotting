#!/usr/bin/env python
from setuptools import setup

# Note basemap 1.1 should be back on Pypi.
req = ['nose','numpy','matplotlib','seaborn',
       'pyqt',
       'basemap',
       ]

#%% install
setup(name='MahaliPlotting',
      packages=['MahaliPlotting'],
      url='https://github.com/scivision/weaksig-plot',
      classifiers=[
      'Intended Audience :: Science/Research',
      'Development Status :: 3 - Alpha',
      'Topic :: Scientific/Engineering :: Atmospheric Science',
      'Programming Language :: Python',
      ],
      install_requires=req,
      dependency_links=['https://downloads.sourceforge.net/project/matplotlib/matplotlib-toolkits/basemap-1.0.7/basemap-1.0.7.tar.gz'],
	  )