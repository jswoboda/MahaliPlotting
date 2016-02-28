#!/usr/bin/env python
"""
Created on Sun Feb 28 15:56:10 2016

@author: John Swoboda
"""

#!/usr/bin/env python
"""

"""
import os, glob,getopt,sys
import scipy as sp
import ConfigParser

import matplotlib
matplotlib.use('Agg') # for use where you're running on a command line
import pdb
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pytz
from datetime import datetime
from dateutil import parser
from mpl_toolkits.basemap import Basemap
from GeoData.plotting import scatterGD, slice2DGD,insertinfo
from GeoData.GeoData import GeoData
from GeoData.utilityfuncs import readIonofiles, readAllskyFITS,readSRI_h5

INIOPTIONS = ['latbounds','lonbounds','timebounds','timewin','date','asgamma','aslim','gpslim','isrparams','paramlim','isrheight']
class PlotClass(object):
    def __init__(self,inifile,GPSloc=None,ASloc=None,ISRloc=None):
        self.inifile = inifile
        self.params = readini(inifile)
        self.GeoData={'ISR':None,'GPS':None,'AS':None}
        
        self.GPSRead(GPSloc)
        self.ASRead(ASloc)
        self.ISRRead(ISRloc)
        
        self.RegisterData()
        
    def RegisterData():
        
        
        
    def GPSRead(self,GPSloc):
        if GPSloc is None:
            return
    def ASRead(self,ASloc):
        if ASloc is None:
            return  
    def ISRRead(self,ISRloc):
        if ISRloc is None:
            return
            
            
def readini(inifile):
    