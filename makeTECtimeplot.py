#!/usr/bin/env python
"""
Created on Thu Mar  3 14:45:35 2016

@author: swoboj
"""

import os, glob,getopt,sys
import scipy as sp
import matplotlib
matplotlib.use('Agg') # for use where you're running on a command line
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from GeoData.plotting import scatterGD, slice2DGD,insertinfo
from GeoData.GeoData import GeoData
from GeoData.utilityfuncs import readIonofiles, readAllskyFITS,readSRI_h5
from plotdata import plottecvstime


if __name__== '__main__':

    datadir = 'TempGeo'
    flist1 = glob.glob(os.path.join(datadir,'*.h5'))
    TEClist1= map(GeoData.read_h5,flist1)
    flist=[]
    TEClist = []
    satnum=23
        
    for i,j in enumerate(TEClist1):
        if sp.any(j.data['satnum']==satnum):
            TEClist.append(j)
            flist.append(flist1[i])
        
    col = 2.
    numr = sp.ceil(len(flist)/col)
    


    dxs = 4.0
    dys = 2.
    fig, axmat = plt.subplots(int(numr),int(col),dpi=300,sharey=True,figsize=(dxs*col,dys*(numr+1)))
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    axvec= axmat.flatten()
    
    dnames = [os.path.splitext(os.path.split(i)[-1])[0] for i in flist]
    for i,iGD in enumerate(TEClist):
        lines = plottecvstime(iGD,satnum,fig,axvec[i])
        axvec[i].set_title(dnames[i])
    plt.suptitle('Data from Sat: {0:d}'.format(satnum))
    plt.subplots_adjust(top=0.85)
    plt.savefig('TECMaps')