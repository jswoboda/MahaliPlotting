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

INIOPTIONS = ['latbounds','lonbounds','timebounds','timewin','date','asgamma','aslim','gpslim','isrparams','paramlim','isrheight','reinterp','paramheight']
class PlotClass(object):
    def __init__(self,inifile,GPSloc=None,ASloc=None,ISRloc=None):
        self.inifile = inifile
        self.params = readini(inifile)
        self.GDISR = None
        self.GDGPS = None
        self.GDAS = None
        
        #GeoData objects
        self.GPSRead(GPSloc)
        self.ASRead(ASloc)
        self.ISRRead(ISRloc)
        
        self.RegisterData()
        
    
        
    def GPSRead(self,GPSloc):
        if GPSloc is None:
            return
        
        timelim=self.params['timebounds']
        TEClist = []
        TECfiles = glob.glob(os.path.join(GPSloc,'*.iono'))
        TECtime = [sp.Inf,-sp.Inf];

        for ifile in TECfiles:
            TECGD = GeoData(readIonofiles,(ifile,))
            if timelim is not None:
                TECGD.timereduce(timelim)
                
            if len(TECGD.times)==0:
                continue
            TEClist.append(TECGD)
            TECtime[0] = min(min(TECGD.times[:,0]),TECtime[0])
            TECtime[1] = max(max(TECGD.times[:,0]),TECtime[1])
            
            self.GDGPS = TECGD
        
    def ASRead(self,ASloc):
        if ASloc is None:
            return  
            
        wlstr ='*_0'+wl+'_*.FITS'
        interpsavedfile = os.path.join(ASloc,'interp'+wl+'.h5')
        reinterp=self.params['reinterp']
        timelim=self.params['timebounds']
        if reinterp or (not os.path.isfile(interpsavedfile)):
            pfalla = sp.array([65.136667,-147.447222,689.])
    
            flist558 = glob.glob(os.path.join(ASloc,wlstr))
            allsky_data = GeoData(readAllskyFITS,(flist558,'PKR_20111006_AZ_10deg.FITS','PKR_20111006_EL_10deg.FITS',150.,pfalla))
            if timelim is not None:
                allsky_data.timereduce(timelim)


            xcoords = allsky_data.__changecoords__('WGS84')
            latlim=[xcoords[:,0].min(),xcoords[:,0].max()]
            lonlim=[xcoords[:,1].min(),xcoords[:,1].max()]
            nlat = 256
            nlon = 256
    
            latvec = sp.linspace(latlim[0],latlim[1],nlat)
            lonvec = sp.linspace(lonlim[0],lonlim[1],nlon)
            [LATM,LONM] = sp.meshgrid(latvec,lonvec)
    
            newcoords = sp.column_stack((LATM.flatten(),LONM.flatten(),150.*sp.ones(LONM.size)))
            allsky_data.interpolate(newcoords,'WGS84',method='linear',twodinterp=True)
            allsky_data.write_h5(interpsavedfile)
        else:
            allsky_data = GeoData.read_h5(interpsavedfile)
            if timelim is not None:
                allsky_data.timereduce(timelim)
            
        self.GDAS = allsky_data
    def ISRRead(self,ISRloc):
        
        if ISRloc is None:
            return
         
        pnheights= self.params['paramheight']
        paramstr = self.params['isrparams']
        SRIh5 = GeoData(readSRI_h5,(ISRloc,paramstr))
        dt1ts,dt2ts = self.params['timebounds']
    
        timelist = sp.where((SRIh5.times[:,0]>=dt1ts)&(SRIh5.times[:,0]<=dt2ts))[0]
    
        if len(timelist)==0:
            return
    
        SRIh5 = SRIh5.timeslice(timelist)
    
        hset = sp.array([i[1] for i in pnheights])
        uh,uhs =sp.unique(hset,return_inverse=True)
    
        newcoordname = 'WGS84'
        
        changed_coords = SRIh5.__changecoords__(newcoordname)
        
        latmin,latmax = [changed_coords[:,0].min(),changed_coords[:,0].max()]
        lonmin,lonmax = [changed_coords[:,1].min(),changed_coords[:,1].max()]
        latvec = sp.linspace(latmin,latmax,self.params['ISRLatnum'])
        lonvec = sp.linspace(lonmin,lonmax,self.params['ISRLonnum'])
        
        LON,LAT = sp.meshgrid(lonvec,latvec)
        xycoords = [LAT.flatten(),LON.flatten()]
            # interpolation
        ncoords = xycoords.shape[0]
        uhall = sp.repeat(uh,ncoords)
    
        coords = sp.tile(xycoords,(len(uh),1))
        coords = sp.column_stack((coords,uhall))
    
    
        SRIh5.interpolate(coords,newcoordname,method='linear')
        self.GDISR = SRIh5
    def RegisterData():
        """ """
        
            
def readini(inifile):
    