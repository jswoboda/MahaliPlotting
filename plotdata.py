#!/usr/bin/env python
"""  

"""
import os, glob
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib
import pdb
from mpl_toolkits.basemap import Basemap, cm
from GeoData.plotting import scatterGD, slice2DGD
from GeoData.GeoData import GeoData 
from GeoData.utilityfuncs import readIonofiles, readAllskyFITS


def main(allskydir,ionofdir,plotdir,wl = str(558),tint=5,reinterp=False):
    
    wlstr ='*_0'+wl+'_*.FITS'
    
    
    
    TEClist = []
    TECfiles = glob.glob(os.path.join(ionofdir,'*.iono'))
    TECtime = [sp.Inf,-sp.Inf];
    Geolatlim = [sp.Inf,-sp.Inf];
    Geolonlim = [sp.Inf,-sp.Inf];
    for ifile in TECfiles:
        TECGD = GeoData(readIonofiles,(ifile,))
        TEClist.append(TECGD)
        TECtime[0] = min(min(TECGD.times[:,0]),TECtime[0])
        TECtime[1] = max(max(TECGD.times[:,0]),TECtime[1])
        Geolatlim[0] = min(min(TECGD.dataloc[:,0]),Geolatlim[0])
        Geolatlim[1] = max(max(TECGD.dataloc[:,0]),Geolatlim[1])
    
        Geolonlim[0] = min(min(TECGD.dataloc[:,1]),Geolonlim[0])
        Geolonlim[1] = max(max(TECGD.dataloc[:,1]),Geolonlim[1])
    
    latlim2 = [45.,75.]
    lonlim2 = [-185.,-125.]
    interpsavedfile = os.path.join(allskydir,'interp'+wl+'.h5')
    if reinterp or (not os.path.isfile(interpsavedfile)):
        pfalla = sp.array([65.136667,-147.447222,689.])
    
        flist558 = glob.glob(os.path.join(allskydir,wlstr))
        allsky_data = GeoData(readAllskyFITS,(flist558,'PKR_DASC_20110112_AZ_10deg.FITS','PKR_DASC_20110112_EL_10deg.FITS',150.,pfalla))
        
            # reduce the size of the allskydata
        allskytime = allsky_data.times[:,0]
        allsky_data=allsky_data.timeslice(sp.where(sp.logical_and(allskytime>=TECtime[0],allskytime<TECtime[1] ))[0])
        allskytime=allsky_data.times[:,0]
        latlim2 = [45.,75.]
        lonlim2 = [-185.,-125.]
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
        allskytime=allsky_data.times[:,0]
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_axes([0.1,0.1,0.8,0.8])
    # create polar stereographic Basemap instance.
    m = Basemap(projection='merc',lon_0=sp.mean(lonlim2),lat_0=sp.mean(latlim2),\
        lat_ts=sp.mean(latlim2),llcrnrlat=latlim2[0],urcrnrlat=latlim2[1],\
        llcrnrlon=lonlim2[0],urcrnrlon=lonlim2[1],\
        rsphere=6371200.,resolution='i')
    # draw coastlines, state and country boundaries, edge of map.
    #m.drawcoastlines()
#    m.drawstates()
#    m.drawcountries()
    shp_info = m.readshapefile('st99_d00','states',drawbounds=True)
    meridians=sp.arange(lonlim2[0],lonlim2[1],10.)
    parallels = sp.arange(latlim2[0],latlim2[1],5.)
    parhand=m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
    mrdhand = m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
    
    tectime = sp.arange(TECtime[0],TECtime[1],60.*tint)
    nptimes= len(tectime)    
    
     
    allskylist = []
    gpslist = []
    tlistkeep = sp.zeros(nptimes-1,dtype=bool)
    for itasn in range(nptimes-1):
        techlist = []
        tlistemp=True
        itback=tectime[itasn]
        itfor = tectime[itasn+1]
        itas = sp.where(sp.logical_and(allskytime>=itback, allskytime<itfor))[0]
        if len(itas)==0:
            continue
        
        for k in range(len(TEClist)):
            Geoone=TEClist[k];
            timevec = Geoone.times[:,0];
            
            itgps = sp.where(sp.logical_and(timevec>=itback, timevec<itfor))[0]
            if len(itgps)>0:
                tlistemp=False
            techlist.append(itgps)
        if tlistemp:
            continue
        allskylist.append(itas) 
        gpslist.append(techlist)
        tlistkeep[itasn]=True
    plotgpsnoptics(allsky_data,TEClist,allskylist,gpslist,plotdir,m,ax,fig)
    plt.close(fig)
        
def plotgpsnoptics(optic_class,gps_class_list,opticallist,GPSlist,plotdir,m,ax,fig):
    
    maxplot = len(optic_class.times)
    strlen = sp.ceil(sp.log10(maxplot))+1
    fmstr = '{0:0>'+str(strlen)+'}_'
    plotnum=0
    for (optic_times,gps_cur)in zip(opticallist,GPSlist):
        gpshands = []        
        for igpsn, (igps,igpslist) in enumerate(zip(gps_class_list,gps_cur)):
            print('Plotting GPS data from rec {0} of {1}'.format(igpsn,len(gps_cur)))
            # check if there's anything to plot
            if len(igpslist)==0:
                continue
            
            (sctter,scatercb) = scatterGD(igps,'alt',3.5e5,vbounds=[0,15],time = igpslist,gkey = 'vTEC',cmap='jet',fig=fig,
                  ax=ax,title='',cbar=True,err=.1,m=m)
            gpshands.append(sctter)
            
        #change he z order
        minz = gpshands[0]
        for i in reversed(gpshands):
            i.set_zorder(i.get_zorder()+1)
        
        for iop in optic_times:
        
            (slice3,cbar3) = slice2DGD(optic_class,'alt',150,[100,800],title='',
                                time = iop,cmap='Greys',gkey = 'image',fig=fig,ax=ax,cbar=False,m=m)
            slice3.set_zorder(minz)
            print('Ploting {0} of {1} plots'.format(plotnum,maxplot))
            plt.savefig(os.path.join(plotdir,fmstr.format(plotnum)+'ASwGPS.png'))
            plotnum+=1
            slice3.remove()
        for i in reversed(gpshands):
            i.remove()
    
if __name__== '__main__':
    plotdir = os.path.expanduser('~/Documents/python/mahali/plots10172015')
    if not os.path.isdir(plotdir):
        os.mkdir(plotdir)
    allskydir = os.path.expanduser('~/DATA/Mahali/allsky')
    ionofdir = os.path.expanduser('~/DATA/Mahali/GPS_Processed_IonoFiles/site_280_2015')
    cmdrun=True
    if cmdrun:
        matplotlib.use('Agg') # for use where you're running on a command line

    main(allskydir,ionofdir,plotdir,tint=2.)
