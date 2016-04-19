#!/usr/bin/env python
"""
Created on Sun Feb 28 15:56:10 2016

@author: John Swoboda
"""

import os, glob, sys, getopt
import scipy as sp
import ConfigParser

#import matplotlib
#matplotlib.use('Agg') # for use where you're running on a command line
import pdb
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pytz
from datetime import datetime
from dateutil import parser
from mpl_toolkits.basemap import Basemap
from GeoData.plotting import scatterGD, slice2DGD,insertinfo, contourGD
from GeoData.GeoData import GeoData
from GeoData.utilityfuncs import readIonofiles, readAllskyFITS,readSRI_h5
from copy import copy
INIOPTIONS = ['latbounds','lonbounds','timebounds','timewin','asgamma','aslim','gpslim','paramlim','reinterp','paramheight','ISRLatnum','ISRLonnum','wl','TextList']

class PlotClass(object):
    """ This class will handle all of the reading , registration and plotting of data
    related to the mahali project. This is ment to be both the back end to a gui and 
    the basic code used to run a whole set of data. The class will be given an 
    ini file and location of data sets. The ini file will hold information regaurding
    the plotting and data registration.
    Variables:
        inifile - The name of the ini file used for the the parameters.
        params - A dixtionary that holds all of the parameters. The list INIOPTIONS
            holds all of the keys for this dicitonary.
        GDISR - A GeodData object for the ISR data. The default value is None.
        GDGPS - A GeodData object for the GPS data. The default value is None.
        GDAS - A GeodData object for the AllSky data. The default value is None.
        numGD - The number of Sensor classes avalible.
        GPSNames - The names of the GPS receivers.
        Regdict - A dictionary used to order all of the data. The keys are 
            TEC, AS, Time ISR. """
    def __init__(self,inifile,GPSloc=None,ASloc=None,ISRloc=None):
        """ This will create an instance of a PlotClass object.
            Inputs
            inifile - The name of the ini file used for the the parameters.
            GPSloc - The directory that holds all of the iono files. 
            ASloc - This can be either a directory holding the FITS files or a 
                h5 file thats been pre interpolated.
            ISRloc - This can be either a file from SRI or an h5 file
                thats been pre interpolated."""
        self.inifile = inifile
        self.params = readini(inifile)
        # GeoData objects
        self.GDISR = None
        self.GDGPS = None
        self.GPSNames = None
        self.GDAS = None
        self.numGD = 0
        # Read in GeoData objects
        self.GPSRead(GPSloc)
        self.ASRead(ASloc)
        self.ISRRead(ISRloc)
        # Time Register everyting
        self.RegisterData()
    
  #%% Read in data      
    def GPSRead(self,GPSloc):
        """ This function will read in the GPS data from ionofiles. It will assign 
            the class variable GDGPS to the resultant GeoData object.
            Input
                GPSloc - The directory that holds all of the iono files. """
        if GPSloc is None:
            return
        
        if not os.path.isdir(os.path.expanduser(GPSloc)):
            print('GPS path is not a directory')
            return
        self.numGD+=1
        print('Reading in GPS Data')
        timelim=self.params['timebounds']
        TEClist = []
        TECfiles = glob.glob(os.path.join(GPSloc,'*.iono'))

        for ifile in TECfiles:
            TECGD = GeoData(readIonofiles,(ifile,))
            if timelim is not None:
                TECGD.timereduce(timelim)
                
            if len(TECGD.times)==0:
                continue
            TEClist.append(TECGD)
            
        # Determine the reciver names
        self.GPSNames = [os.path.splitext(os.path.split(i)[-1])[0].split('-')[0] for i in TECfiles]
        self.GDGPS = TEClist
        print('Finished Reading in GPS Data')
        
    def ASRead(self,ASloc):
        """ This function will read in the All sky data from FITS files or structured
            h5 files for GeoData.
            Input
                ASloc - This can be either a directory holding the FITS files or a 
                h5 file thats been pre interpolated. It will assign the class variable GDAS to the
            resultant GeoData object."""
        
        if ASloc is None:
            return  
        if not os.path.isdir(os.path.expanduser(ASloc)) and not os.path.isfile(os.path.expanduser(ASloc)):
            print('All Sky Data cannot be read')
            return
        self.numGD+=1
        print('Reading in All sky data')
        wl = str(int(self.params['wl']))
        wlstr ='*_0'+wl+'_*.FITS'
        interpsavedfile = os.path.join(ASloc,'interp'+wl+'.h5')
        reinterp=self.params['reinterp']
        timelim=self.params['timebounds']
        if reinterp or (not os.path.isfile(ASloc)):
            pfalla = sp.array([65.136667,-147.447222,689.])
    
            filestr = os.path.join(ASloc,wlstr)
            flist558 = glob.glob(filestr)
            if len(flist558)==0:
                return
            allsky_data = GeoData(readAllskyFITS,(flist558,('PKR_20111006_AZ_10deg.FITS','PKR_20111006_EL_10deg.FITS'),150.,timelim))
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
            allsky_data = GeoData.read_h5(ASloc)
            if timelim is not None:
                allsky_data.timereduce(timelim)
            
        self.GDAS = allsky_data
        print('Finished Reading in Allsky Data')
    def ISRRead(self,ISRloc):
        """ This function will read in the ISR data from SRI files or structured
            h5 files for GeoData. It will assign the class variable GDISR to the
            resultant GeoData object.
            Input
                ISRloc - This can be either a file from SRI or an h5 file
                thats been pre interpolated. """
        if ISRloc is None:
            return
        if not os.path.isfile(os.path.expanduser(ISRloc)):
            print('ISR Data is not a file')
            return
        self.numGD+=1
        print('Reading in ISR Data') 
        pnheights= self.params['paramheight']
        
        paramstr = list(set([i[0] for i in pnheights]))
        try:
            SRIh5 = GeoData(readSRI_h5,(ISRloc,paramstr))
        except:
            try:
                SRIh5 = GeoData.read_h5(ISRloc)
            except:
                return                    
        dt1ts,dt2ts = self.params['timebounds']
    
        timelist = sp.where((SRIh5.times[:,1]>=dt1ts)&(SRIh5.times[:,0]<=dt2ts))[0]
    
        if len(timelist)==0:
            return
    
        SRIh5 = SRIh5.timeslice(timelist)
    
        hset = sp.array([i[1] for i in pnheights])
        uh,uhs =sp.unique(hset,return_inverse=True)
        uh=uh*1e3
        newcoordname = 'WGS84'
        
        if not (SRIh5.coordnames.lower() == newcoordname.lower()):
            changed_coords = SRIh5.__changecoords__(newcoordname)
            latmin,latmax = [changed_coords[:,0].min(),changed_coords[:,0].max()]
            lonmin,lonmax = [changed_coords[:,1].min(),changed_coords[:,1].max()]
            latvec = sp.linspace(latmin,latmax,self.params['ISRLatnum'])
            lonvec = sp.linspace(lonmin,lonmax,self.params['ISRLonnum'])
            
            LON,LAT = sp.meshgrid(lonvec,latvec)
            xycoords = sp.column_stack([LAT.flatten(),LON.flatten()])
                # interpolation
            ncoords = xycoords.shape[0]
            uhall = sp.repeat(uh,ncoords)
        
            coords = sp.tile(xycoords,(len(uh),1))
            coords = sp.column_stack((coords,uhall))
        
        
            SRIh5.interpolate(coords,newcoordname,method='linear')
        self.GDISR = SRIh5
        print('Finished Reading in ISR Data')
        #%% Registration
    def RegisterData(self):
        """ This function will register data and return a dictionary that will be used
            keep track of the time associations between all of the different data sets.
            """
        tbounds = self.params['timebounds']
        tint=self.params['timewin']
        # make lists for plotting
        tectime = sp.arange(tbounds[0],tbounds[1],60.*tint)
        
        nptimes= len(tectime)
        timelists = [[tectime[i],tectime[i+1]] for i in range(nptimes-1)]
        # teclist is a Ntime length list
        teclist = [[]]*(nptimes-1) 
        regdict = {'TEC':teclist,'AS':[None]*(nptimes-1),'ISR':[None]*(nptimes-1),'Time':timelists}
        
        if not self.GDGPS is None:
           for itasn in range(len(tectime)-1):        
                itback=tectime[itasn]
                itfor = tectime[itasn+1 ]
                templist = [[]]*len(self.GDGPS)
                for k in range(len(self.GDGPS)):
                    Geoone=self.GDGPS[k]
                    timevec = Geoone.times[:,0]
        
                    templist[k] = sp.where(sp.logical_and(timevec>=itback, timevec<itfor))[0]
                teclist[itasn]  =templist      
            
        if (not self.GDAS is None):
            GPS2AS=[sp.array([])]*(nptimes-1)
            GPS2ASlen = [0]*(nptimes-1)
            allskytime=self.GDAS.times[:,0]
            
            
            for itasn in range(len(tectime)-1):    
                itback=tectime[itasn]
                itfor = tectime[itasn+1]
                itas = sp.where(sp.logical_and(allskytime>=itback, allskytime<itfor))[0]
                if len(itas)==0:
                    itas = sp.where(allskytime<=itback)[0]
                    if len(itas)==0:
                        continue
                    itas = sp.array([itas[-1]])
                GPS2AS[itasn] = itas
                GPS2ASlen[itasn]=len(itas)
                
            GPS2ASsingle = []
            teclist2 =[]
            timelist2 = []
            #repeat for all sky values
            for i1,ilen in enumerate(GPS2ASlen):
                GPS2ASsingle=GPS2ASsingle+GPS2AS[i1].tolist()
                teclist2 =teclist2+[ teclist[i1]]*ilen
                timelist2 =timelist2+[ timelists[i1]]*ilen
            regdict['TEC']=teclist2
            regdict['AS']=GPS2ASsingle
            regdict['Time']=timelist2
            regdict['ISR'] = [None]*len(GPS2ASsingle)
            if (not self.GDISR is None):
                self.GDAS=self.GDAS.timeslice(sp.unique(GPS2ASsingle))
                as2radar =self.GDAS.timeregister(self.GDISR)
                
                teclist3=[]
                timelist3=[]
                GPS2ASsingle2=[]
                AS2ISRsingle=[]
                for j1,jasval in enumerate(as2radar):
                    jlen=len(jasval)
                    AS2ISRsingle =AS2ISRsingle +jasval.tolist()
                    teclist3=teclist3+[teclist2[j1]]*jlen
                    timelist3=timelist3+[timelist2[j1]]*jlen
                    GPS2ASsingle2 = GPS2ASsingle2+[GPS2ASsingle[j1]]*jlen
                
                regdict['TEC']=teclist3
                regdict['AS']=GPS2ASsingle2
                regdict['Time']=timelist3
                regdict['ISR'] = AS2ISRsingle
                
        elif (not self.GDISR is None):
            GPS2ISR=[sp.array([])]*(nptimes-1)
            GPS2ISRlen = [0]*(nptimes-1)
            isrtime=self.GDISR.times
            
           
            for itasn in range(len(tectime)-1):          
                itback=tectime[itasn]
                itfor = tectime[itasn+1]
                #need to fix this
                itas = sp.where(sp.logical_and(isrtime[:,1]>=itback, isrtime[:,0]<itfor))[0]
                if len(itas)==0:
                    itas = sp.where(isrtime<=itback)[0]
                    if len(itas)==0:
                        continue
                    itas = sp.array([itas[-1]])
                GPS2ISR[itasn] = itas
                GPS2ISRlen[itasn]=len(itas)
            
            GPS2ISRsingle = []
            teclist2 =[]
            timelist2 = []
            
            #repeat for all sky values
            for i1,ilen in enumerate(GPS2ISRlen):
                GPS2ISRsingle=GPS2ASsingle+GPS2ISR[i1].tolist()
                teclist2 =teclist2+[ teclist[i1]]*ilen
                timelist2 =timelist2+[ timelists[i1]]*ilen
            regdict['TEC']=teclist2
            regdict['ISR']=GPS2ISRsingle
            regdict['Time']=timelist2
            regdict['AS']=[None]*len(timelist2)
        
        self.Regdict=regdict
#%% Plotting
    def plotmap(self,fig,ax):
        """ This function will plot the map of Alaska. The data will be plotted
            over it and will use the basemap class to position everything.
            Input
                fig - The figure handle for the plots.
                ax - The axes handle that the map will be plotted over.
            Output
                m - This is the handle for the basemap object.
        """
        latlim2 = self.params['latbounds']
        lonlim2 = self.params['lonbounds']
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
        return m
        
    def plotalldata(self,plotdir,m,ax,fig):
        """ This method will plot all of the images from the time period specified.
            The png files will be numbered and have the ISR parameter name in 
            the title. The user can easily manipulate the files using comand line.
            Inputs 
                plotdir - This is the directory that all of the plots will be saved 
                    in.
                m - The map object that will be used to orient the data.
                fig - The figure handle for the plots.
                ax - The axes handle that the map will be plotted over.
        """
        
        timelist = self.Regdict['Time']
        Nt = len(timelist)
        strlen = int(sp.ceil(sp.log10(Nt))+1)
        fmstr = '{0:0>'+str(strlen)+'}_'
        plotnum=0
        cbarax = []
        if self.params['paramheight'] is None or self.GDISR is None:
            Ncase = 1
            paramstrs = ['']
        else:
            Ncase = len(self.params['paramheight'])
            paramstrs = [i[0]+str(int(i[1])) for i in self.params['paramheight']]
        Nplot = Ncase*Nt
        for itime in range(Nt):
            for icase in range(Ncase):
                hands,cbarax = self.plotsingle(m,ax,fig,itime,icase,cbarax)
                print('Ploting {0} of {1} plots'.format(plotnum,Nplot))
                plt.savefig(os.path.join(plotdir,fmstr.format(plotnum)+paramstrs[icase]+'ASwGPS.png'))
                
                if itime in self.params['TextList']:
                    hands = self.plotgpsnames(m,ax,fig,hands,timenum=itime)
                    plt.savefig(os.path.join(plotdir,'wloctext' + fmstr.format(plotnum)+paramstrs[icase]+'ASwGPS.png'))
                
                for ihand in hands:
                    if hasattr(ihand, "__len__"):
                        for ihand2 in ihand:
                            ihand2.remove()
                    elif hasattr(ihand,'collections'):
                        for ihand2 in ihand.collections:
                            ihand2.remove()
                    else:
                        ihand.remove()
                plotnum+=1
        # write out ini file to record plot parameters
        ininame = os.path.split(self.inifile)[-1]
        writeini(self.params,os.path.join(plotdir,ininame))
            
    def plotsingle(self,m,ax,fig,timenum=0,icase=0,cbarax=[]):
        """ Make single plot given the desired time and number associated with
            the desired ISR param.
            Inputs
                m - The map handle that is used to plot everything.
                fig - The figure handle for the plots.
                ax - The axes handle that the map will be plotted over.
                timenum - The of GeoData objects derived from the ionofiles.
                icase - A list of list which determines which allsky times are used.
                cbarax - The list of color bar axes.
            Outputs
                allhands - The list handles of the plotted data. 
                cbarax - A list of 
        """
        

        optbnds = self.params['aslim']
        gam=self.params['asgamma']
        
        curwin=self.Regdict['Time'][timenum]
        
        allhands = [[]]
        
        titlelist = []
        if len(cbarax)==0:
            wid = .3/self.numGD
            cbarax=[fig.add_axes([.7+i*wid,.3,wid/2.,.4]) for i in range(self.numGD)]
            fig.tight_layout(rect=[0,.05,.7,.95])
        cbcur=0
        if not self.GDGPS is None:
            
           gpshands = []
           gpsbounds = self.params['gpslim']
           for igps,igpslist in zip(self.GDGPS,self.Regdict['TEC'][timenum]):
               # check if there's anything to plot
               if len(igpslist)==0:
                   continue
               (sctter,scatercb) = scatterGD(igps,'alt',3.5e5,vbounds=gpsbounds,time = igpslist,gkey = 'vTEC',cmap='plasma',fig=fig, ax=ax,title='',cbar=False,err=.1,m=m)
                    
               gpshands.append(sctter)
           
           
           
             
           # If no gps data plots dont try to plot the color bar
           if len(gpshands)>0:
               scatercb = plt.colorbar(sctter,cax=cbarax[cbcur])
               
               scatercb.set_label('vTEC in TECu')
           cbcur+=1 
           allhands[0]=gpshands
           titlelist.append( insertinfo('GPS $tmdy $thmsehms',posix=curwin[0],posixend=curwin[1]))
           #change he z order
           
           allhands[0]=gpshands
        
        if not self.GDAS is None:
            
            iop = self.Regdict['AS'][timenum]
            (slice3,cbar3) = slice2DGD(self.GDAS,'alt',150,optbnds,title='',
                                time = iop,cmap='gray',gkey = 'image',fig=fig,ax=ax,cbar=False,m=m)
            slice3.set_norm(colors.PowerNorm(gamma=gam,vmin=optbnds[0],vmax=optbnds[1]))
            titlelist.append(insertinfo('All Sky $tmdy $thmsehms',posix=self.GDAS.times[iop,0],posixend=self.GDAS.times[iop,1]))
            cbaras = plt.colorbar(slice3,cax=cbarax[cbcur])
            cbcur+=1
            cbaras.set_label('All Sky Scale')
            minz=slice3.get_zorder()
            for i in reversed(allhands[0]):
                minz=sp.minimum(minz,i.get_zorder())
                i.set_zorder(i.get_zorder()+1)
            slice3.set_zorder(minz)
            allhands.append(slice3)
        if not self.GDISR is None:
            itn=self.Regdict['ISR'][timenum]
            curph=self.params['paramheight'][icase]
            vbounds = self.params['paramlim']
            iparam=curph[0]
            if iparam.lower()!='ne':
                levels = sp.logspace(sp.log10(vbounds[icase][0]),sp.log10(vbounds[icase][1]),5)
            else:
                levels = sp.linspace(vbounds[icase][0],vbounds[icase][1],5)
            (plth,cbh) =  contourGD(self.GDISR,'alt',curph[1],vbounds=vbounds[icase],time = itn,gkey = iparam,cmap='jet',fig=fig,
                  ax=ax,cbar=False,m=m,levels = levels)
            cbh = plt.colorbar(plth,cax=cbarax[cbcur])
            cbh.set_label(iparam)
            if iparam.lower()!='ne':
               
                fmt= '%d'
                
                
            else:

                plth.set_norm(colors.LogNorm(vmin=vbounds[icase][0],vmax=vbounds[icase][1]))
                fmt = '%.1e'
            titlelist.append( insertinfo('ISR Data at $tmdy $thmsehms',posix=self.GDISR.times[itn,0],posixend = self.GDISR.times[itn,1]))
#            minz=plth.get_zorder()
#            for i in reversed(allhands[0]):
#                minz=sp.minimum(minz,i.get_zorder)
#                i.set_zorder(i.get_zorder()+1)
#            plth.set_zorder(minz)
            cbh = plt.colorbar(plth,cax=cbarax[cbcur],format=fmt)
            cbh.set_label(iparam)
            allhands.append(plth)
        ax.set_title('\n'.join(titlelist) )
        return allhands,cbarax
    def plotgpsnames(self,m,ax,fig,allhands,timenum=0):
        """ This will plot the names of the GPS recievers next to the pierce points.
            Inputs
                m - The map handle that is used to plot everything.
                fig - The figure handle for the plots.
                ax - The axes handle that the map will be plotted over.
                allhands - The list handles of the plotted data. 
                timenum - The of GeoData objects derived from the ionofiles.                
            Outputs
                allhands - The list handles of the plotted data. 
            """
        if self.GDGPS is None:
            return allhands
        xoffset = 0.022*(m.xmax-m.xmin)
        yoffset = 0.022*(m.ymax-m.ymin)
        namehands = []
        for igpsnum, (igps,igpslist) in enumerate(zip(self.GDGPS,self.Regdict['TEC'][timenum])):
           # check if there's anything to plot
           if len(igpslist)==0:
               continue
           locs = igps.dataloc[igpslist]
           xwin,ywin = m(self.params['lonbounds'],self.params['latbounds'])
           x, y = m(locs[:,1], locs[:,0])
           keepboth =  (x>xwin[0])&(x<xwin[1])&(y>ywin[0])&(y<ywin[1])
           if keepboth.sum()==0:
               continue
           [xloc,yloc] = [x[keepboth][0],y[keepboth][0]]
           namehands.append(plt.text(xloc+xoffset,yloc+yoffset,self.GPSNames[igpsnum],fontsize=14,fontweight='bold',
                    color='r'))
        allhands.append(namehands)
        return allhands
        
        
    def writeiniclass(self,fname):
        """ The method for the class that calls the writeini function.
            Inputs
                fname - The name of the file it will be written out to.
        """
        writeini(self.params,fname)
#%% Write out file
def writeini(params,fname):
    """ This will write out a structured ini file that can be used by the PlotClass
        to fill its dictionary of parameters. 
        Inputs
            params - A dictionary that holds the parameters for the PlotClass.
            fname - The name of the file it will be written out to.
    """
    cfgfile = open(fname,'w')
    config = ConfigParser.ConfigParser(allow_no_value = True)
    
    config.add_section('params')
    config.add_section('paramsnames')
    for ip in INIOPTIONS:
        
        if not ip in params.keys():
            continue
        elif ip=='timebounds':
            config.set('params',ip,' '.join(posix2str(params[ip])))
        elif ip=='paramheight':
            temp= [item for sublist in params[ip] for item in sublist]
            data = ""
            for a in temp:
                data += str(a)
                data += " "
            config.set('params',ip,data)
        elif ip=='paramlim':
            temp= [item for sublist in params[ip] for item in sublist]
            data = ""
            for a in temp:
                data += str(a)
                data += " "
            config.set('params',ip,data)
        elif ip=='reinterp':
            if params[ip]:
                data='Yes'
            else:
                data='No'
            config.set('params',ip,data)
        elif type(params[ip]) in (sp.ndarray,list):
            data = ""
            for a in params[ip]:
                data += str(a)
                data += " "
            config.set('params',ip,data)
        else:
            config.set('params',ip,str(params[ip]))
        config.set('paramsnames',ip,ip)
    config.write(cfgfile)
    cfgfile.close()
    #%% Read in file        
def readini(inifile):
    """ This function will read in data from a configureation file.
        Inputs
            inifile- The name of the configuration file.
        Outputs 
            params - A dictionary with keys from INIOPTIONS that holds all of
                the plotting parameters.
    """
    config = ConfigParser.ConfigParser()
    config.read(inifile)
    params={i:None for i in INIOPTIONS}
    # Read in data from ini file
    for ip in config.options('params'):
        # get the original param name
        rname  = config.get('paramsnames',ip)
        # get the parameter and split it up
        params[rname] = config.get('params',ip)
        params[rname]=params[rname].split(" ")
        # If its a single object try to 
        if len(params[rname])==1:
            params[rname]=params[rname][0]
            try:
                params[rname]=float(params[rname])
            except:
                pass
        else:
            for a in range(len(params[rname])):
                try:
                    params[rname][a]=float(params[rname][a])
                except:
                    pass
    
    # turn the time bounds to time stamps
    if not params['timebounds']is None:
        timelist = params['timebounds']
        params['timebounds']=str2posix(timelist)
    # which times will have names
    if params['TextList'] is None:
        params['TextList']=[]
    
        
        
    # change param height to a list of lists 
    if not params['paramheight'] is None:
        l1 = params['paramheight'][::2]
        l2 = params['paramheight'][1::2]
        params['paramheight']=[[i,j] for i,j in zip(l1,l2)]
    if not params['paramlim'] is None:
        l1 = params['paramlim'][::2]
        l2 = params['paramlim'][1::2]
        params['paramlim']=[[i,j] for i,j in zip(l1,l2)]
    # Default for reinterp is false
    if params['reinterp']is None:
        params['reinterp']=False
    else:
        params['reinterp'] = params['reinterp'].lower()=='yes'
    return params

def str2posix(timelist):
    """ This will take a list of strings with the date along with a start and
        end time and make a list with the posix times.
        Inputs
            timelist - A list of strings with the data followed by two times. 
            The date for the second time can also be used, it will be at index
            2 and the second time will be at index 3.
        Outputs
            dtts - A list of posix times from the original inputs"""
    if len(timelist)==3:
        timelist.insert(2,timelist[0])
        
    (dt1,dt2) = parser.parse(timelist[0]+ ' '+timelist[1]),parser.parse(timelist[2]+ ' '+timelist[3])
    dt1 =dt1.replace(tzinfo=pytz.utc)
    dt2 = dt2.replace(tzinfo=pytz.utc)
    dt1ts = (dt1 -datetime(1970,1,1,0,0,0,tzinfo=pytz.utc)).total_seconds()
    dt2ts = (dt2 -datetime(1970,1,1,0,0,0,tzinfo=pytz.utc)).total_seconds()
    return [dt1ts,dt2ts]
    
def posix2str(posixlist):
    """ This function will chnage a list of posix times to a string in 
        %m/%d/%Y %H:%M:%S year format.
        Inputs
            posixlist - Length 2 list of posix times.
        Outputs
            data - A string of times.  
    """
    dts = map(datetime.utcfromtimestamp, posixlist)
    data = [datetime.strftime(dts[0],'%m/%d/%Y %H:%M:%S'), datetime.strftime(dts[1],'%m/%d/%Y %H:%M:%S')]
    return data

def runPlotClass(inifile,plotdir, gpsloc=None,ASloc=None,ISRloc=None):
    """ This will create a figure and axis and use the PlotClass to make images 
        of the data.
        Inputs
         GPSloc - The directory that holds all of the iono files. 
        ASloc - This can be either a directory holding the FITS files or a 
            h5 file thats been pre interpolated.
        ISRloc - This can be either a file from SRI or an h5 file
            thats been pre interpolated.
    """
    (fig,axmat) = plt.subplots(1,1,figsize=(16,12),facecolor='w')
    PC = PlotClass(inifile,GPSloc=gpsloc,ASloc=ASloc,ISRloc=ISRloc)
    m=PC.plotmap(fig,axmat)
    PC.plotalldata(plotdir,m,axmat,fig)
    plt.close(fig)
if __name__== '__main__':
    
    from argparse import ArgumentParser
    descr = '''
             This script will run Mahali Plotting software. The user needs to
             specify the locations of the different types of data and the time
             limits if they don't want all of the data processed. This specific
             module runs everything in a class struture so it can be used with
             other modules easier.
            '''
    p = ArgumentParser(description=descr)
    p.add_argument('-i','--ifile',help='The directory that holds all of the TEC data in ionofile formats.',default=None)
    p.add_argument("-a", "--asky",help='The allsky data directory.',default=None)
    p.add_argument("-p", "--pdir",help='plot output directory',default=os.getcwd())

    p.add_argument('-r', "--radar",help='Radar hdf5 file',default=None)#action='store_true')
 
    p.add_argument('-c',"--config",help='Config file name.')
    p = p.parse_args()
    
    if p.ifile is None:
        gpsloc=p.ifile
    else:
        gpsloc=os.path.expanduser(p.ifile)
    if p.asky is None:
        ASloc=p.asky
    else:
        ASloc=os.path.expanduser(p.asky)
    if p.radar is None:
        ISRloc=p.radar
    else:
        ISRloc=os.path.expanduser(p.radar)
   
    plotdir=os.path.expanduser(p.pdir)
    inifile = os.path.expanduser(p.config)
    
    runPlotClass(inifile,plotdir,gpsloc,ASloc,ISRloc)