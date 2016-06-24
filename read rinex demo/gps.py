"""
Greg Starr
this is based mostly on Bill Rideout's tec.py
scripts, I'm not using classes and I am using
numpy because I like it better, the rinex
reading was made by Michael Hirsch and Greg Starr
"""
from __future__ import division,absolute_import
import numpy as np
from datetime import datetime, timedelta
from pandas import DataFrame,Panel,Series,Panel4D
from pandas.io.pytables import read_hdf
from os.path import splitext,expanduser,getsize
from io import BytesIO
import time

f1 = 1575.42E6 #MHz
f2 = 1227.6E6  #MHz
def getRanges(data,svn,maxgap=3,maxjump=2.0): 
    if c2p2(data,svn):
        nans = np.logical_or.reduce((np.isnan(data[:,svn,'L1','data']),
                                     np.isnan(data[:,svn,'L2','data']),
                                     np.isnan(data[:,svn,'C1','data']),
                                     np.isnan(data[:,svn,'C2','data'])))
    else:
        nans = np.logical_or.reduce((np.isnan(data[:,svn,'L1','data']),
                                     np.isnan(data[:,svn,'L2','data']),
                                     np.isnan(data[:,svn,'C1','data']),
                                     np.isnan(data[:,svn,'P2','data'])))
    inarc=False
    start=[]
    end=[]
    phase=2.85E9*(data[:,svn,'L1','data']/f1-data[:,svn,'L2','data']/f2)
    lgi=0

    for i in range(len(nans)):
        if inarc:
            if nans[i]:
                if i-lgi>maxjump:
                    inarc=False
                    end.append(lgi)
            else:
                if abs(phase[i]-phase[lgi])>maxjump:
                    end.append(lgi)
                    start.append(i)
                lgi=i

        else:
            if not nans[i]:
                inarc=True
                start.append(i)
                lgi=i
                
    ranges = [(data.labels[a],data.labels[b]) for a,b in zip(start,end)]

    return ranges

def getTec(data,svn,drange,satbias=None):
    if c2p2(data,svn,drange):
        diffRange = (2.85E9/3.0E8)*(
            data[drange[0]:drange[1],svn,'C2','data']
            -data[drange[0]:drange[1],svn,'C1','data'])
    else:
        diffRange = (2.85E9/3.0E8)*(
            data[drange[0]:drange[1],svn,'P2','data']
            -data[drange[0]:drange[1],svn,'C1','data'])
        
    phase=2.85E9*(data[drange[0]:drange[1],svn,'L1','data']/f1
                  -data[drange[0]:drange[1],svn,'L2','data']/f2)        
        
    diffList = sorted(phase-diffRange)
    medianDiff = diffList[int(len(diffList)/2)]
    distWidth = diffList[int(len(diffList)*.75)]-diffList[int(len(diffList)*.25)]
    medianErr = distWidth/np.sqrt(len(diffList))
    tec = phase - medianDiff
    if satbias!=None: # FIX FOR SAT BIAS
        tec-=satbias.dict[svn]
    return tec,medianErr

def c2p2(data,svn,drange=(None,None)):
    return (np.sum(~np.isnan(
        data[drange[0]:drange[1],svn,'C2','data']))>
            np.sum(~np.isnan(
                data[drange[0]:drange[1],svn,'P2','data'])))

def c2p2alt(data,site,drange=(None,None)):
    return (np.sum(~np.isnan(
        data[site,drange[0]:drange[1],'C2','data']))>
            np.sum(~np.isnan(
                data[site,drange[0]:drange[1],'P2','data'])))



def getTecalt(data,site,drange,satbias=None):
    if c2p2alt(data,site,drange):
        diffRange = (2.85E9/3.0E8)*(
            data[site,drange[0]:drange[1],'C2','data']
            -data[site,drange[0]:drange[1],'C1','data'])
    else:
        diffRange = (2.85E9/3.0E8)*(
            data[site,drange[0]:drange[1],'P2','data']
            -data[site,drange[0]:drange[1],'C1','data'])
        
    phase=2.85E9*(data[site,drange[0]:drange[1],'L1','data']/f1
                  -data[site,drange[0]:drange[1],'L2','data']/f2)        
        
    diffList = sorted(phase-diffRange)
    medianDiff = diffList[int(len(diffList)/2)]
    distWidth = diffList[int(len(diffList)*.75)]-diffList[int(len(diffList)*.25)]
    medianErr = distWidth/np.sqrt(len(diffList))
    tec = phase - medianDiff
    if satbias!=None:
        tec-=satbias
    return tec,medianErr

def getRangesalt(data,site,maxgap=3,maxjump=2.0): 
    if c2p2alt(data,site):
        nans = np.logical_or.reduce((np.isnan(data[site,:,'L1','data']),
                                     np.isnan(data[site,:,'L2','data']),
                                     np.isnan(data[site,:,'C1','data']),
                                     np.isnan(data[site,:,'C2','data'])))
    else:
        nans = np.logical_or.reduce((np.isnan(data[site,:,'L1','data']),
                                     np.isnan(data[site,:,'L2','data']),
                                     np.isnan(data[site,:,'C1','data']),
                                     np.isnan(data[site,:,'P2','data'])))
    inarc=False
    start=[]
    end=[]
    phase=2.85E9*(data[site,:,'L1','data']/f1-data[site,:,'L2','data']/f2)
    lgi=0

    for i in range(len(nans)):
        if inarc:
            if nans[i]:
                if i-lgi>maxjump:
                    inarc=False
                    end.append(lgi)
            else:
                if abs(phase[i]-phase[lgi])>maxjump:
                    end.append(lgi)
                    start.append(i)
                lgi=i

        else:
            if not nans[i]:
                inarc=True
                start.append(i)
                lgi=i
                
    ranges = [(data.items[a],data.items[b]) for a,b in zip(start,end)]

    return ranges
            
class satelliteBias:
    """satelliteBias is a class to get satellite biases in tec units
    Once biases are loaded, get them using the dictionary attribute
    dict.  Key is tuple of prn (integer) and biasType (integer).  If
    TEC is calculated using C1, set biasType to 1.  If
    TEC is calculated using P1, set biasType to 0. If TEC is calculated
    using C1 and C2, set biasType to 2.
    """

    def __init__(self, satFile, C1BiasFile, L2C2BiasFile):
        """__init__ sets up the dictionary self.dict
        satFile - the ionex file with satellite biases as produced by
            JPL in ftp://cddis.gsfc.nasa.gov/pub/gps/products/ionex/
        C1BiasFile - the P1C1 bias file (may be None for verification only)
        L2C2BiasFile - the P2C2 bias file (may be None for verification only)
        """
        self.dict = {}
        self.__parseSatBiasFile(satFile)
        self.__parseC1BiasFile(C1BiasFile)
        self._parseC2BiasFile(L2C2BiasFile)


    def __parseSatBiasFile(self, satFile):
        """__parseSatBiasFile parses satellite bias file, and adds data
        to self.dict
        """
        indicatorStr = 'DIFFERENTIAL CODE BIASES'
        # conversionFactor in TECu
        conversionFactor = -0.463*6.158 # diff ns -> meters -> tec
        f = open(satFile)
        lines = f.readlines()
        f.close()
        lineFound = 0 # indicates right line found
        dataFound = 0 # indicates at least one line of data found
        for line in lines:
            if line[0:len(indicatorStr)] == indicatorStr:
                lineFound = 1
                continue
            if lineFound:
                items = line.split()
                # see if we're done
                try:
                    try:
                        sv = int(items[0])
                    except:
                        # see if last two characters are ints and first is G
                        if items[0][0] == 'G':
                            sv = int(items[0][-2:])
                        else:
                            raise(IOError, '')
                    bias = float(items[1])*conversionFactor
                    dataFound = 1
                    # add this data to dict
                    self.dict[sv] = bias
                except:
                    if dataFound == 0:
                        # no valid lines found
                        raise(IOError,
                              'No valid data found after indicator in %s'
                              % (satFile))
                    else:
                        return
        # if we got here, the indicator wasn't found, or the data was the last line in the file
        if dataFound == 1:
            return
        else:
            raise(IOError,
                  'No indicator string found in %s' % (satFile))

    def __parseC1BiasFile(self, C1BiasFile):
        """__parseC1BiasFile parses p1c1 bias file, and adds data
        to self.dict
        Bias is added to existing biases
        """
        conversionFactor = -0.463*6.158 # diff ns -> meters -> tec
        # allow no C1BiasFile for case where normal bias just being verified
        if (C1BiasFile == None):
            return
        f = open(C1BiasFile)
        lines = f.readlines()
        f.close()
        # print warning if no data found
        dataFound = False
        for line in lines:
            try:
                items = line.split()
                if items[0][0] in ('G', 'g'):
                    prn = int(items[0][1:])
                else:
                    prn = int(items[0])
                addBias =  float(items[1])* conversionFactor
                self.dict[(prn, 1)] =  self.dict[(prn, 0)] - addBias
                dataFound = True
            except:
                continue
        if not dataFound:
            print('WARNING: No valid data found in %s' % (C1BiasFile))

    def _parseC2BiasFile(self, L2C2BiasFile):
        """__parseC2BiasFile parses p2c2 bias file, and adds data
        to self.dict
        Bias is added to existing biases
        """
        conversionFactor = -0.463*6.158 # diff ns -> meters -> tec
        # allow no C1BiasFile for case where normal bias just being verified
        if (L2C2BiasFile == None):
            return
        f = open(L2C2BiasFile)
        lines = f.readlines()
        f.close()
        # print warning if no data found
        dataFound = False
        for line in lines:
            try:
                items = line.split()
                if items[0][0] in ('G', 'g'):
                    prn = int(items[0][1:])
                else:
                    prn = int(items[0])
                addBias =  float(items[1])* conversionFactor
                self.dict[(prn, 2)] =  self.dict[(prn, 1)] + addBias
                dataFound = True
            except:
                continue
        if not dataFound:
            print('WARNING: No valid data found in %s' % (L2C2BiasFile))




def rinexobs(obsfn,writeh5=None,maxtimes=None):
    stem,ext = splitext(expanduser(obsfn))
    if ext[-1].lower() == 'o': #raw text file
        with open(obsfn,'r') as f:
            t=time.time()
            lines = f.read().splitlines(True)
            lines.append('')
            header,version,headlines,obstimes,sats,svset = scan(lines)
            print('{} is a RINEX {} file, {} kB.'.format(obsfn,version,getsize(obsfn)/1000.0))
            data = processBlocks(lines,header,obstimes,svset,headlines,sats)
            print("finished in {0:.2f} seconds".format(time.time()-t))
    #%% save to disk (optional)
        if writeh5:
            h5fn = stem + '.h5'
            print('saving OBS data to {}'.format(h5fn))
            data.to_hdf(h5fn,key='OBS',mode='a',complevel=6,append=False)
    elif ext.lower() == '.h5':
        data = read_hdf(obsfn,key='OBS')
        print('loaded OBS data from {} to {}'.format(blocks.items[0],blocks.items[-1]))
    return data


# this will scan the document for the header info and for the line on
# which each block starts
def scan(lines):
    header={}        
    eoh=0
    for i,line in enumerate(lines):
        if "END OF HEADER" in line:
            eoh=i
            break
        if line[60:].strip() not in header:
            header[line[60:].strip()] = line[:60].strip()
        else:
            header[line[60:].strip()] += " "+line[:60].strip()
    verRinex = float(header['RINEX VERSION / TYPE'].split()[0])
    header['APPROX POSITION XYZ'] = [float(i) for i in header['APPROX POSITION XYZ'].split()]
    header['# / TYPES OF OBSERV'] = header['# / TYPES OF OBSERV'].split()
    header['# / TYPES OF OBSERV'][0] = int(header['# / TYPES OF OBSERV'][0])
    header['INTERVAL'] = float(header['INTERVAL'])
        
    headlines=[]
    obstimes=[]
    sats=[]
    svset=set()
    i=eoh+1
    while True:
        if not lines[i]: break
        if not int(lines[i][28]):
            #no flag or flag=0
            headlines.append(i)
            obstimes.append(_obstime([lines[i][1:3],lines[i][4:6],
                                   lines[i][7:9],lines[i][10:12],
                                   lines[i][13:15],lines[i][16:26]]))
            numsvs = int(lines[i][30:32])
            if(numsvs>12):
                sats.append([int(lines[i][33+s*3:35+s*3]) for s in range(12)])
                i += 1
                sats.append([int(lines[i][33+s*3:35+s*3]) for s in range(numsvs-12)])
            else:
                sats.append([int(lines[i][33+s*3:35+s*3]) for s in range(numsvs)])
        
            i+=numsvs*int(np.ceil(header['# / TYPES OF OBSERV'][0]/5))+1
        else:
            #there was a comment or some header info
            flag=int(lines[i][28])
            if(flag!=4):
                print(flag)
            skip=int(lines[i][30:32])
            i+=skip+1
    for s in sats:
        svset = svset.union(set(s))

    return header,verRinex,headlines,obstimes,sats,svset



def processBlocks(lines,header,obstimes,svset,headlines,sats):
    obstypes = header['# / TYPES OF OBSERV'][1:]
    blocks = Panel4D(labels=obstimes,
                     items=list(svset),
                     major_axis=obstypes,
                     minor_axis=['data','lli','ssi'])
    ttime1 = 0
    ttime2 = 0
    for i in range(len(headlines)):
        linesinblock = len(sats[i])*int(np.ceil(header['# / TYPES OF OBSERV'][0]/5))
        block = ''.join(lines[headlines[i]+1:headlines[i]+linesinblock+1])
        t1 = time.time()
        bdf = _block2df(block,obstypes,sats[i],len(sats[i]))
        ttime1 += (time.time()-t1)
        t2 = time.time()
        blocks.loc[obstimes[i],sats[i]] = bdf
        ttime2 += (time.time()-t2)            
    print("{0:.2f} seconds for _block2df".format(ttime1))
    print("{0:.2f} seconds for panel assignments".format(ttime2))
    return blocks       
        

def _obstime(fol):
    year = int(fol[0])
    if 80<= year <=99:
        year+=1900
    elif year<80: #because we might pass in four-digit year
        year+=2000
    return datetime(year=year, month=int(fol[1]), day= int(fol[2]),
                    hour= int(fol[3]), minute=int(fol[4]),
                    second=int(float(fol[5])),
                    microsecond=int(float(fol[5]) % 1) *100000
                    )

def _block2df(block,obstypes,svnames,svnum):
    """
    input: block of text corresponding to one time increment INTERVAL of RINEX file
    output: 2-D array of float64 data from block. Future: consider whether best to use Numpy, Pandas, or Xray.
    """
    nobs = len(obstypes)
    stride=3

    strio = BytesIO(block.encode())
    barr = np.genfromtxt(strio, delimiter=(14,1,1)*5).reshape((svnum,-1), order='C')

    data = barr[:,0:nobs*stride:stride]
    lli  = barr[:,1:nobs*stride:stride]
    ssi  = barr[:,2:nobs*stride:stride]

    data = np.vstack(([data.T],[lli.T],[ssi.T])).T

    return data


def readRinexNav(fn,writeh5=None):
    """
    Michael Hirsch
    It may actually be faster to read the entire file via f.read() and then .split()
    and asarray().reshape() to the final result, but I did it frame by frame.
    http://gage14.upc.es/gLAB/HTML/GPS_Navigation_Rinex_v2.11.html
    """
    stem,ext = splitext(expanduser(fn))
    startcol = 3 #column where numerical data starts
    nfloat=19 #number of text elements per float data number
    nline=7 #number of lines per record

    with open(expanduser(fn),'r') as f:
        #find end of header, which has non-constant length
        while True:
            if 'END OF HEADER' in f.readline(): break
        #handle frame by frame
        sv = []; epoch=[]; raws=''
        while True:
            headln = f.readline()
            if not headln: break
            #handle the header
            sv.append(headln[:2])
            year = int(headln[2:5])
            if 80<= year <=99:
                year+=1900
            elif year<80: #good till year 2180
                year+=2000
            epoch.append(datetime(year =year,
                                  month   =int(headln[5:8]),
                                  day     =int(headln[8:11]),
                                  hour    =int(headln[11:14]),
                                  minute  =int(headln[14:17]),
                                  second  =int(headln[17:20]),
                                  microsecond=int(headln[21])*100000))
            """
            now get the data.
            Use rstrip() to chomp newlines consistently on Windows and Python 2.7/3.4
            Specifically [:-1] doesn't work consistently as .rstrip() does here.
            """
            raw = (headln[22:].rstrip() +
                   ''.join(f.readline()[startcol:].rstrip() for _ in range(nline-1))
                   +f.readline()[startcol:40].rstrip())

            raws += raw + '\n'

    raws = raws.replace('D','E')

    strio = BytesIO(raws.encode())
    darr = np.genfromtxt(strio,delimiter=nfloat)

    nav= DataFrame(np.hstack((np.asarray(sv,int)[:,None],darr)), epoch,
               ['sv','SVclockBias','SVclockDrift','SVclockDriftRate','IODE',
                'Crs','DeltaN','M0','Cuc','Eccentricity','Cus','sqrtA','TimeEph',
                'Cic','OMEGA','CIS','Io','Crc','omega','OMEGA DOT','IDOT',
                'CodesL2','GPSWeek','L2Pflag','SVacc','SVhealth','TGD','IODC',
                'TransTime','FitIntvl'])

    if writeh5:
        h5fn = stem + '.h5'
        print('saving NAV data to {}'.format(h5fn))
        nav.to_hdf(h5fn,key='NAV',mode='a',complevel=6,append=False)

    return nav

