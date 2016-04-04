#!/usr/bin/env python
import os, glob, datetime,sys,getopt
from GeoData.GeoData import GeoData
import matplotlib.pyplot as plt
import numpy as np
from GeoData.utilityfuncs import readIonofiles, readAllskyFITS,readSRI_h5
from PlottingClass import str2posix
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter,MinuteLocator, HourLocator

def PlotTECdiff(gpsloc,timelist,satnum,sublist,pname='mahdiff.png',tectype = 'TEC'):
    """ This will plot the differences between two gps sites. The code filters
        first by satilite number. The TEC is interpolated over the time period
        and then subtracted.
        Inputs
        gpsloc - The directory holding the ionofiles.
        timelist - A list of strings holding the time the data will be interpolated over
            [date1,time1,date2,time1]
            
        satnum - The satilite number as an int or float.
        sublist - This is the name of the recievers that will have their values compared
        pname - Name of the plot
        tectype - Either vTEC or TEC, determines which type will be compared."""
        
    # Figure out location or receivers
    flist1 = glob.glob(os.path.join(gpsloc,'*.iono'))
    
    fnames = np.array([os.path.splitext(os.path.split(i)[1])[0].lower().split('-')[0] for i in flist1])
    
    [f0,f1] = [np.argwhere(i.lower()==fnames)[0][0] for i in sublist]
    
    # read in the receivers and filter out by sat number
    mah0str = flist1[f0]
    mah1str = flist1[f1]
    
    mah0 = GeoData(readIonofiles,(mah0str,))
    mah1 = GeoData(readIonofiles,(mah1str,))
    
    sat230 = mah0.data['satnum']==satnum
    sat231 = mah1.data['satnum']==satnum

    timemah0 = mah0.times[sat230]
    timemah1 = mah1.times[sat231]

    TEC0 = mah0.data[tectype][sat230]
    TEC1 = mah1.data[tectype][sat231]
    
    # Interpolation
    xends = str2posix(timelist)
    xint = np.linspace(xends[0],xends[1],180)
    mah0int = np.interp(xint,timemah0[:,0],TEC0)
    mah1int = np.interp(xint,timemah1[:,0],TEC1)
    
    mahdif = mah0int-mah1int
    
    #plotting
    dts = map(datetime.datetime.utcfromtimestamp, xint)
    fig, axmat = plt.subplots(1,1,dpi=300)
    
    lines = axmat.plot(dts,mahdif)
    axmat.set_title(sublist[0]+' - '+sublist[1] +' '+ timelist[0])
    dtfmt = DateFormatter('%H:%M:%S')
    axmat.xaxis.set_major_locator(HourLocator())
    axmat.xaxis.set_major_formatter(dtfmt)
    axmat.set_xlabel('Time UT')
    axmat.set_ylabel(tectype)
    plt.savefig(pname)
    

if __name__== '__main__':
    
    argv = sys.argv[1:]
    outstr = ''' 
             Usage: mahlidiff.py  -i <ionofile dir>, -s <sat number>, -b <begining time>, -e <endtime>, -p <plotname> -m <first receiver> -n <second receiver>

             or 
             
             python plotdata.py -h
             
             This script will run comparethe TEC in two mahali recievers 
                
            Optional arguments

            -i The directory that holds all of the TEC data in ionofile formats.
            -s Satlilite number
            -b begining time.
            -e endtime
            -p plotdir
            -m first receiver
            -n seocnd receiver
            
            '''


    try:
        opts, args = getopt.gnu_getopt(argv,"hs:i:b:e:p:m:n:")
    except getopt.GetoptError:
        print(outstr)
        sys.exit(2)    
    
    gpsloc=''
    ASloc=''
    ISRloc=''
    
    inifile=None
    sublist = ['','']
    for opt, arg in opts:
        if opt == '-h':
            print(outstr)
            sys.exit()
        elif opt in ("-i", "--ifile"):
            gpsloc = os.path.expanduser(arg)   
        elif opt in ("-s", "--sat"):
            satnum=int(float(arg))
        elif opt in ("-b", "--btime"):
            btime=os.path.expanduser(arg)
        elif opt in ("-e", "--etime"):
            endtime = os.path.expanduser(arg) 
        elif opt in ('-p','--pdir'):
            plotdir=os.path.expanduser(arg)
        elif opt in ('-m','--minu'):
            sublist[0] = arg
        elif opt in ('-n','--ninu'):
            sublist[1] = arg
            
    timelist = btime.split(' ') +endtime.split(' ')
    PlotTECdiff(gpsloc,timelist,satnum,sublist)