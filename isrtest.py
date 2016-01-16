#!/usr/bin/env python
"""
Created on Fri Jan 15 12:48:45 2016

@author: John Swoboda
"""

import os
import numpy as np
from plotdata import getSRIhdf5



if __name__== '__main__':
    isrdir= os.path.expanduser('~/DATA/Mahali/ISRdata/2015_10_11')
    filename= os.path.join(isrdir,'20151011.004_lp_3min-cal.h5')
    xvec = np.linspace(-150.,150.,128)
    yvec = np.linspace(-50.,200.,128)
    Xmat,Ymat=np.meshgrid(xvec,yvec)
    xycoords = np.column_stack((Xmat.flatten(),Ymat.flatten()))
    newcordname='Cartesian'
    vbounds=[[1e10,5e11],[1000,3000],[1000,3000]]
    times = ['2015-10-11 10:00:00','2015-10-11 12:00:00']
    pnheights=[('Ne',150.),('Ti',250.),('Te',250.)]
    getSRIhdf5(filename,times,pnheights,xycoords,newcordname,vbounds,'/Users/Bodangles/DATA/Mahali/ISRdata/2015_10_11/Plots')