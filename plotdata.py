#!/usr/bin/env python
"""

"""
import os, glob,getopt,sys
import scipy as sp
import matplotlib
matplotlib.use('Agg') # for use where you're running on a command line
import pdb
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter,MinuteLocator, HourLocator
import pytz
from datetime import datetime
from dateutil import parser
from mpl_toolkits.basemap import Basemap
from GeoData.plotting import scatterGD, slice2DGD,insertinfo
from GeoData.GeoData import GeoData
from GeoData.utilityfuncs import readIonofiles, readAllskyFITS,readSRI_h5


def main(allskydir,ionofdir,plotdir,latlim2,lonlim2,wl = str(558),tint=5,reinterp=False,timelim=None):
    """ This is the main function for plot data. This function will determine what is to be plotted
    and call all of the spatial and time regestration programs from GeoData.
    Inputs
        allskydir - The directory that holds the FITS files for the allsky data.
            If a None is passed then a the allsky is not plotted.
        ionofdir - The directory that holds all of the ionofiles. If a None is 
            passed then a the allsky is not plotted.
        plotdir - The directory where the plots are stored.
        wl - The wavelength of the allsky light in a string.
        tint - The number of minutes the GPS data is plotted over. Default is 5.
        reinterp - A bool that determines if the GeoData files for the optical data
            should be remade.
        timelim - A list that shows the time boundries in posix."""
#    latlim2 = [45.,75.]
#    lonlim2 = [-185.,-125.]
    
    
        #%% Make map
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_axes([0.1,0.1,0.8,0.8])
    # create polar stereographic Basemap instance.
    m = Basemap(projection='merc',lon_0=sp.mean(lonlim2),lat_0=sp.mean(latlim2),\
        lat_ts=sp.mean(latlim2),llcrnrlat=latlim2[0],urcrnrlat=latlim2[1],\
        llcrnrlon=lonlim2[0],urcrnrlon=lonlim2[1],\
        rsphere=6371200.,resolution='i',ax=ax)
    # draw coastlines, state and country boundaries, edge of map.
    #m.drawcoastlines()
#    m.drawstates()
#    m.drawcountries()
    shp_info = m.readshapefile('st99_d00','states',drawbounds=True)
    
    merstep = sp.round_((lonlim2[1]-lonlim2[0])/5.)
    parstep = sp.round_((latlim2[1]-latlim2[0])/5.)
    meridians=sp.arange(lonlim2[0],lonlim2[1],merstep)
    parallels = sp.arange(latlim2[0],latlim2[1],parstep)
    parhand=m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
    mrdhand = m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
    plt.hold(True)
    isallsky=False
    isgps = False
    if ionofdir is not  None:
        isgps=True
        TEClist = []
        TECfiles = glob.glob(os.path.join(ionofdir,'*.iono'))
        TECtime = [sp.Inf,-sp.Inf];
    #    Geolatlim = [sp.Inf,-sp.Inf];
    #    Geolonlim = [sp.Inf,-sp.Inf];
        for ifile in TECfiles:
            TECGD = GeoData(readIonofiles,(ifile,))
            if timelim is not None:
                TECGD.timereduce(timelim)
                
            if len(TECGD.times)==0:
                continue
            TEClist.append(TECGD)
            TECtime[0] = min(min(TECGD.times[:,0]),TECtime[0])
            TECtime[1] = max(max(TECGD.times[:,0]),TECtime[1])

    if allskydir is not None:
        isallsky=True
   
        wlstr ='*_0'+wl+'_*.FITS'
        interpsavedfile = os.path.join(allskydir,'interp'+wl+'.h5')
        if reinterp or (not os.path.isfile(interpsavedfile)):
            pfalla = sp.array([65.136667,-147.447222,689.])
    
            flist558 = glob.glob(os.path.join(allskydir,wlstr))
            allsky_data = GeoData(readAllskyFITS,(flist558,'PKR_20111006_AZ_10deg.FITS','PKR_20111006_EL_10deg.FITS',150.,pfalla))
            if timelim is not None:
                allsky_data.timereduce(timelim)
                # reduce the size of the allskydata
            allskytime = allsky_data.times[:,0]
            allsky_data=allsky_data.timeslice(sp.where(sp.logical_and(allskytime>=TECtime[0],allskytime<TECtime[1] ))[0])
            allskytime=allsky_data.times[:,0]

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
            allskytime=allsky_data.times[:,0]

    
    #%% make lists for plotting
    tectime = sp.arange(TECtime[0],TECtime[1],60.*tint)
    nptimes= len(tectime)

    if isallsky and isgps:
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
                itas = sp.where(allskytime<=itback)[0]
                if len(itas)==0:
                    continue
                itas = [itas[-1]]
    
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

        plotgpswoptics(allsky_data,TEClist,allskylist,gpslist,plotdir,m,ax,fig,latlim2,lonlim2)
    elif isgps:
        gpslist = []
        tlistkeep = sp.zeros(nptimes-1,dtype=bool)
        for itasn in range(nptimes-1):
            techlist = []
            tlistemp=True
            itback=tectime[itasn]
            itfor = tectime[itasn+1]
                
            for k in range(len(TEClist)):
                Geoone=TEClist[k];
                timevec = Geoone.times[:,0];
    
                itgps = sp.where(sp.logical_and(timevec>=itback, timevec<itfor))[0]
                if len(itgps)>0:
                    tlistemp=False
                techlist.append(itgps)
            if tlistemp:
                continue
            gpslist.append(techlist)
            tlistkeep[itasn]=True
        plotgpsonly(TEClist,gpslist,plotdir,m,ax,fig,latlim2,lonlim2)
    elif isallsky:            
        plotopticsonly(allsky_data,plotdir,m,ax,fig,latlim2,lonlim2)
    
    plt.close(fig)
    
def plotgpsonly(TEClist,gpslist,plotdir,m,ax,fig,latlim,lonlim):
    """ Makes a set of plots when only gps data is avalible."""
    maxplot = len(gpslist)
    strlen = int(sp.ceil(sp.log10(maxplot))+1)
    fmstr = '{0:0>'+str(strlen)+'}_'
    plotnum=0
    for gps_cur in gpslist:
        gpshands = []
        gpsmin = sp.inf
        gpsmax = -sp.inf
        for igpsn, (igps,igpslist) in enumerate(zip(TEClist,gps_cur)):
            print('Plotting GPS data from rec {0} of {1}'.format(igpsn,len(gps_cur)))
            # check if there's anything to plot
            if len(igpslist)==0:
                continue

            (sctter,scatercb) = scatterGD(igps,'alt',3.5e5,vbounds=[0,20],time = igpslist,gkey = 'vTEC',cmap='plasma',fig=fig,
                  ax=ax,title='',cbar=True,err=.1,m=m)
            gpsmin = sp.minimum(igps.times[igpslist,0].min(),gpsmin)
            gpsmax = sp.maximum(igps.times[igpslist,0].max(),gpsmax)
            
            gpshands.append(sctter)
        scatercb.set_label('vTEC in TECu')
        #change he z order
        
        print('Ploting {0} of {1} plots'.format(plotnum,maxplot))
        plt.savefig(os.path.join(plotdir,fmstr.format(plotnum)+'GPSonly.png'))
        plotnum+=1
        for i in reversed(gpshands):
            i.set_zorder(i.get_zorder()+1)
def plotopticsonly(allsky_data,plotdir,m,ax,fig,latlim,lonlim):
    """ Make a set of pots when only all sky is avalible."""
    maxplot = len(allsky_data.times)
    strlen = int(sp.ceil(sp.log10(maxplot))+1)
    fmstr = '{0:0>'+str(strlen)+'}_'
    optictimes = allsky_data.times
    plotnum=0
    firstbar = True
    optbnds = [300,1100]
    for iop in range(len(optictimes)):
        (slice3,cbar3) = slice2DGD(allsky_data,'alt',150,optbnds,title='',
                            time = iop,cmap='gray',gkey = 'image',fig=fig,ax=ax,cbar=True,m=m)
                            
        slice3.set_norm(colors.PowerNorm(gamma=0.6,vmin=optbnds[0],vmax=optbnds[1]))
        if firstbar:
            firstbar=False
            cbaras = plt.colorbar(slice3,ax=ax,orientation='horizontal')
            cbaras.set_label('All Sky Scale')
        
        plt.title(insertinfo('All Sky $tmdy $thmsehms',posix=allsky_data.times[iop,0],posixend=allsky_data.times[iop,1]))
        print('Ploting {0} of {1} plots'.format(plotnum,maxplot))
        plt.savefig(os.path.join(plotdir,fmstr.format(plotnum)+'ASonly.png'))

        plotnum+=1
        slice3.remove()
    
def plotgpswoptics(allsky_data,TEClist,allskylist,gpslist,plotdir,m,ax,fig,latlim,lonlim):
    """ Make a set of plots when given both all sky ad GPS are given.
        Inputs
            allsky_data - The all sky data as a GeoData object.
            TEClist - The of GeoData objects derived from the ionofiles.
            allskylist - A list of list which determines which allsky times are used."""
    maxplot = len(allsky_data.times)
    maxplot = sp.array([len(i) for i in allskylist]).sum()
    strlen = int(sp.ceil(sp.log10(maxplot))+1)
    fmstr = '{0:0>'+str(strlen)+'}_'
    plotnum=0
    firstbar = True
    optbnds = [300,1100]
    for (optic_times,gps_cur)in zip(allskylist,gpslist):
        gpshands = []
        gpsmin = sp.inf
        gpsmax = -sp.inf
        for igpsn, (igps,igpslist) in enumerate(zip(TEClist,gps_cur)):
            print('Plotting GPS data from rec {0} of {1}'.format(igpsn,len(gps_cur)))
            # check if there's anything to plot
            if len(igpslist)==0:
                continue

            (sctter,scatercb) = scatterGD(igps,'alt',3.5e5,vbounds=[0,20],time = igpslist,gkey = 'vTEC',cmap='plasma',fig=fig,
                  ax=ax,title='',cbar=True,err=.1,m=m)
            gpsmin = sp.minimum(igps.times[igpslist,0].min(),gpsmin)
            gpsmax = sp.maximum(igps.times[igpslist,0].max(),gpsmax)
            gpshands.append(sctter)
        scatercb.set_label('vTEC in TECu')
        #change he z order
        minz = gpshands[0].get_zorder()
        for i in reversed(gpshands):
            i.set_zorder(i.get_zorder()+1)

        for iop in optic_times:
            (slice3,cbar3) = slice2DGD(allsky_data,'alt',150,optbnds,title='',
                                time = iop,cmap='gray',gkey = 'image',fig=fig,ax=ax,cbar=False,m=m)
            slice3.set_norm(colors.PowerNorm(gamma=0.6,vmin=optbnds[0],vmax=optbnds[1]))
           
            if firstbar:
                firstbar=False
                cbaras = plt.colorbar(slice3,ax=ax,orientation='horizontal')
                cbaras.set_label('All Sky Scale')
            slice3.set_zorder(minz)
            plt.title(insertinfo('GPS $tmdy $thmsehms',posix=gpsmin,posixend=gpsmax)+'\n'+insertinfo('All Sky $tmdy $thmsehms',posix=allsky_data.times[iop,0],posixend=allsky_data.times[iop,1]))
            print('Ploting {0} of {1} plots'.format(plotnum,maxplot))
            plt.savefig(os.path.join(plotdir,fmstr.format(plotnum)+'ASwGPS.png'))

            plotnum+=1
            slice3.remove()
        for i in reversed(gpshands):
            i.remove()


def plottecvstime(TECGD,satnum,fig,ax):
    
    keep = TECGD.data['satnum']==satnum
    times = TECGD.times[:,0][keep]
    vtec = TECGD.data['vTEC'][keep]
    dts = map(datetime.fromtimestamp, times)
    dtfmt = DateFormatter('%H:%M:%S')
    
    lines = ax.plot(dts,vtec)
    
    ax.xaxis.set_major_locator(HourLocator())
    ax.xaxis.set_major_formatter(dtfmt)
    ax.set_xlabel('Time')
    ax.set_ylabel('vTEC')
    ax.set_title('Data From Sat {0:d}'.format(satnum))
    return lines
    
def getSRIhdf5(filename,times,pnheights,xycoords,newcordname,vbounds,pltdir =None):
    """ Plots a set of ISR data in SRI's data format."""
    paramstr = ['Ne','Ti','Te']
    SRIh5 = GeoData(readSRI_h5,(filename,paramstr))
    (dt1,dt2) = parser.parse(times[0]),parser.parse(times[1])
    dt1 =dt1.replace(tzinfo=pytz.utc)
    dt2 = dt2.replace(tzinfo=pytz.utc)
    dt1ts = (dt1 -datetime(1970,1,1,0,0,0,tzinfo=pytz.utc)).total_seconds()
    dt2ts = (dt2 -datetime(1970,1,1,0,0,0,tzinfo=pytz.utc)).total_seconds()

    timelist = sp.where((SRIh5.times[:,0]>=dt1ts)&(SRIh5.times[:,0]<=dt2ts))[0]

    if len(timelist)==0:
        return

    SRIh5 = SRIh5.timeslice(timelist)

    hset = sp.array([i[1] for i in pnheights])
    uh,uhs =sp.unique(hset,return_inverse=True)

        # interpolation
    ncoords = xycoords.shape[0]
    uhall = sp.repeat(uh,ncoords)

    coords = sp.tile(xycoords,(len(uh),1))
    coords = sp.column_stack((coords,uhall))


    SRIh5.interpolate(coords,newcordname,method='linear')


    maxplot = len(timelist)
    strlen = int(sp.ceil(sp.log10(maxplot))+1)
    fmstr = '{0:0>'+str(strlen)+'}_'
    for itn in range(len(timelist)):
        fig, axmat = plt.subplots(nrows=len(pnheights),ncols=1)
        axvec = axmat.flatten()


        for icase,(iparam,iheight) in enumerate(pnheights):
           (plth,cbh) =  slice2DGD(SRIh5,'z',uhs[icase],vbounds=vbounds[icase],time = itn,gkey = iparam,cmap='jet',fig=fig,
                  ax=axvec[icase],title=iparam + ' at {0} km'.format(iheight),cbar=True)
           if iparam.lower()!='ne':
               ntics = sp.linspace(vbounds[icase][0],vbounds[icase][1],5)
               cbh.set_ticks(ntics)
               cbh.formatter.fmt = '%d'
               cbh.update_ticks()
           else:
               ntics = sp.linspace(vbounds[icase][0],vbounds[icase][1],5)
               cbh.set_ticks(ntics)
               cbh.formatter.fmt = '%.1e'
               cbh.update_ticks()
        outstr = insertinfo('ISR Data at $tmdy $thmsehms',posix=SRIh5.times[itn,0],posixend = SRIh5.times[itn,1])
        plt.suptitle(outstr)
        fname = 'SRIData'+fmstr.format(itn)+'.png'
        plt.tight_layout()
        plt.subplots_adjust(top=0.85)
        if not pltdir is None:
            fname=os.path.join(pltdir,fname)
        print 'Plotting '+fname
        plt.savefig(fname)
        plt.close(fig)



if __name__== '__main__':
    argv = sys.argv[1:]
    outstr = ''' 
             Usage: plotdata.py -a <all skypath> -w <wavelength>, -i <ionofile dir>, -t <time interval>, -d <date>, -b <begining time>, -e <endtime>, -p <plotdirectory> -r <type y to reinterpolate all sky data> -s <SRI File>

             or 
             
             python plotdata.py -h
             
             This script will run Mahali Plotting software. The user needs to 
             specify the locations of the different types of data and the time 
             limits if they don't want all of the data processed. 
                
            Optional arguments
            -a The allsky data directory.
            -w The wavelength of the allsky data.
            -i The directory that holds all of the TEC data in ionofile formats.
            -t The time interval that the TEC data will be plotted over.
            -d The date of the data.
            -b The beginning of the time window that the data will be plotted over.
            -e The ending of the time window that the data will be plotted over.
            -l Latitude bounds, must put both in.
            -o Longitude bound, must be both.
            -r If a y follows this then the raw radar data will be remade. If
                this is not used the radar data will only be made if it does
                not exist in the file first.
            
             Example:
             python runsim.py -f radardata -f fitting -i ~/DATA/ExampleLongPulse -c ~/DATA/Example -r y'''


    try:
        opts, args = getopt.gnu_getopt(argv,"ha:w:i:t:r:p:d:b:e:l:o:")
    except getopt.GetoptError:
        print(outstr)
        sys.exit(2)

    remakealldata = False
    allskydir=None
    ionofdir=None
   
    wl='558'
    timelist=[None]*3
    
    latlim2 = [45.,75.]
    lonlim2 = [-175.,-125.]    
    
    latlist=[]
    lonlist=[]
    for opt, arg in opts:
        if opt == '-h':
            print(outstr)
            sys.exit()
        elif opt in ("-i", "--ifile"):
            ionofdir = os.path.expanduser(arg)

        elif opt in ("-w", "--wlen"):
            wl = arg
        elif opt in ("-a", "--asky"):
            allskydir=os.path.expanduser(arg)
        elif opt in ("-t", "--tin"):
            tint=float(arg)
        elif opt in ("-p", "--pdir"):
            plotdir=os.path.expanduser(arg)
        elif opt in ("-d","--date"):
            timelist[0]=arg
        elif opt in ("-b","--begtime"):
            timelist[1]=arg
        elif opt in ("-e","--endtime"):
            timelist[2]=arg
        elif opt in ('-r', "--re"):
            if arg.lower() == 'y':
                remakealldata = True
        elif opt in ('-l', "--lat"):
            latlist.append(float(arg))
        elif opt in ('-o',"--lon"):
            lonlist.append(float(arg))
            
    if len(latlist)>1:
        minlat=min(latlist)
        maxlat=max(latlist)
        latlist=[sp.maximum(minlat,min(latlim2)),sp.minimum(maxlat,max(latlim2))]
        
    if len(lonlist)>1:
        minlon=min(lonlist)
        maxlon=max(lonlist)
        lonlist=[sp.maximum(minlon,min(lonlim2)),sp.minimum(maxlon,max(lonlim2))]
    if None in timelist:
        timelim=None
    else:
        (dt1,dt2) = parser.parse(timelist[0]+ ' '+timelist[1]),parser.parse(timelist[0]+ ' '+timelist[2])
        dt1 =dt1.replace(tzinfo=pytz.utc)
        dt2 = dt2.replace(tzinfo=pytz.utc)
        dt1ts = (dt1 -datetime(1970,1,1,0,0,0,tzinfo=pytz.utc)).total_seconds()
        dt2ts = (dt2 -datetime(1970,1,1,0,0,0,tzinfo=pytz.utc)).total_seconds()
        timelim=[dt1ts,dt2ts]
#    plotdir = os.path.expanduser('~/Documents/python/mahali/plots10172015')
    if not os.path.isdir(plotdir):
        os.mkdir(plotdir)
#    allskydir = os.path.expanduser('~/DATA/Mahali/allsky')
#    ionofdir = os.path.expanduser('~/DATA/Mahali/GPS_Processed_IonoFiles/site_280_2015')
    cmdrun=True
    if cmdrun:
        matplotlib.use('Agg') # for use where you're running on a command line
    
    main(allskydir,ionofdir,plotdir,latlist,lonlist,wl = wl,tint=tint,reinterp=remakealldata,timelim=timelim)
