#!/usr/bin/env python
import os, glob,getopt, sys
import matplotlib
import numpy as np
import pdb
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
import matplotlib.pyplot as plt
import Tkinter as Tk
import tkFileDialog as fd
import ConfigParser
from copy import copy
from PlottingClass import PlotClass, INIOPTIONS, str2posix, posix2str, readini
from GeoData.plotting import insertinfo
#INIOPTIONS = ['latbounds','lonbounds','timebounds','timewin','asgamma','aslim','gpslim','paramlim','reinterp','paramheight','isrlatnum','isrlonnum','wl']


class App():

    def __init__(self,root,fn=None,inputdata = {'GPS':'','ISR':'','AllSky':''}):
        
        self.root=root
        self.root.title("Mahali")
        
        self.fn = fn
        self.m = None
        self.PC = None
        self.allhands=[]
        self.cbarsax = []        
        bd = {'entries':[],'labels':[]}
        self.input = {'ISR':copy(bd),'GPS':copy(bd),'AllSky':copy(bd)}
        self.options = {}
        for op in INIOPTIONS:
            self.options[op]={}

        self.titleframe=Tk.Frame(self.root)
        self.titleframe.grid(row=0,columnspan=3,pady=20,padx=20)
        self.leb=Tk.Label(self.titleframe,text="Mahali GUI",font=("Helvetica", 24))
        self.leb.grid(row=0,columnspan=3)

        self.menubar = Tk.Menu(self.titleframe)
        self.filemenu = Tk.Menu(self.menubar, tearoff=0)
        self.filemenu.add_command(label="Load", command=self.loadfile)
        self.filemenu.add_command(label="Save", command=self.savefile)
        self.menubar.add_cascade(label="File", menu=self.filemenu)
        self.root.config(menu=self.menubar)

        
        self.plotframe=Tk.Frame(self.root,padx=10,pady=10)
        self.plotframe.grid(row=1,column=0, sticky='n')
        self.fig=plt.Figure(figsize=(8, 6), dpi=100)
        self.sp=self.fig.add_subplot(111)
        self.sp.plot(np.arange(100))
        
        self.canvas=FigureCanvasTkAgg(self.fig,master=self.plotframe)
        self.canvas.show()
        self.canvas.get_tk_widget().grid(row=1,column=0)
        self.canvas._tkcanvas.grid(row=2,column=0)

        self.toolbar = NavigationToolbar2TkAgg(self.canvas, self.plotframe)
        self.toolbar.update()
        self.canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=True)
            
        self.optionsframe=Tk.LabelFrame(self.root,text="Options",padx=10,pady=10)
        self.optionsframe.grid(row=1,column=1,sticky='n')
        for irow, field in enumerate(self.input.keys()):
            self.input[field]['entries']=Tk.Entry(self.optionsframe,width=45)
            self.input[field]['entries'].insert(0,inputdata[field])
            self.input[field]['entries'].grid(row=irow+1,column=1,columnspan=2)
            self.input[field]['labels']=Tk.Label(self.optionsframe,text=field)
            self.input[field]['labels'].grid(row=irow+1,column=0)
        # Time menu bar
        self.times={}
        self.times['var'] = Tk.StringVar(root,'None')
        self.times['options'] = ['None']
        self.times['list'] = [0]
        self.times['menu'] = Tk.OptionMenu(self.optionsframe,self.times['var'],tuple(self.times['options']))
        self.times['menu'].grid(row=len(self.input.keys())+1,column=1)
        self.times['label'] = Tk.Label(self.optionsframe,text='Choose Time/Parmameter')
        self.times['label'].grid(row=len(self.input.keys())+1,column=0)
        # Radar Parameter menu
        self.radarparam={}
        self.radarparam['var'] = Tk.StringVar(root,'None')
        self.radarparam['options'] = ['None']
        self.radarparam['list'] = [0]
        self.radarparam['menu'] = Tk.OptionMenu(self.optionsframe,self.radarparam['var'],tuple(self.radarparam['options']))
        self.radarparam['menu'].grid(row=len(self.input.keys())+1,column=2)
        
        # buttons
        self.buttons={}
        self.buttons['ReadIn'] = Tk.Button(self.optionsframe, text="Read In Data", command=self.readindata)
        self.buttons['ReadIn'].grid(row=len(self.input.keys())+2,column=1,sticky='w')
        self.buttons['ReadIn'] = Tk.Button(self.optionsframe, text="Update Plot", command=self.updateplot)
        self.buttons['ReadIn'].grid(row=len(self.input.keys())+2,column=2,sticky='w')
        # Inputs
        self.i=len(self.input.keys())+3
        for field in self.options:
            self.options[field]['entries']=[]
            self.options[field]['labels']=[]
            if field=='reinterp':
                self.options[field]['var'] = Tk.IntVar()
                self.options[field]['entries'].append(Tk.Checkbutton(self.optionsframe,variable=self.options[field]['var']))
                self.options[field]['entries'][0].grid(row=self.i,column=1)
            elif field in ['latbounds','lonbounds','timebounds','aslim','gpslim']:
                self.options[field]['entries'].append(Tk.Entry(self.optionsframe))
                self.options[field]['entries'][0].grid(row=self.i,column=1)
                self.options[field]['entries'].append(Tk.Entry(self.optionsframe))
                self.options[field]['entries'][1].grid(row=self.i,column=2)
            elif field not in ['paramheight','paramlim']:
                self.options[field]['entries'].append(Tk.Entry(self.optionsframe))
                self.options[field]['entries'][0].grid(row=self.i,column=1)
            if field not in ['paramheight','paramlim']:
                self.options[field]['labels'] = [Tk.Label(self.optionsframe,text=field)]
                self.options[field]['labels'][0].grid(row=self.i,column=0)
            self.i+=1

        self.AddParamButton = Tk.Button(master=self.optionsframe,command=self.AddParam,text='Add Parameter')
        self.AddParamButton.grid(row=self.i,column=1)
        self.i+=1
        self.options['paramheight']['labels'].append(Tk.Label(self.optionsframe,text='Parameter'))
        self.options['paramheight']['labels'][0].grid(row=self.i,column=0)
        self.options['paramheight']['labels'].append(Tk.Label(self.optionsframe,text='Height'))
        self.options['paramheight']['labels'][1].grid(row=self.i,column=1)
        self.options['paramlim']['labels'].append(Tk.Label(self.optionsframe,text='Lower Limit'))
        self.options['paramlim']['labels'][0].grid(row=self.i,column=2)
        self.options['paramlim']['labels'].append(Tk.Label(self.optionsframe,text='Upper Limit'))
        self.options['paramlim']['labels'][1].grid(row=self.i,column=3)
        self.i+=1
        self.numparams=0
        self.AddParam()
        
        # Read in stuff from command line
        if not fn is None:
            self.loadfile(fn)
            files = [len(inputdata[i]) for i in inputdata.keys()]
            if np.any(files>0):
                self.readindata()
    def AddParam(self):
        """ This will add ISR parameters and heights to the gui."""
        self.i+=1
        self.options['paramheight']['entries'].append(Tk.Entry(self.optionsframe))
        self.options['paramheight']['entries'][2*self.numparams].grid(row=self.i,column=0)
        self.options['paramheight']['entries'].append(Tk.Entry(self.optionsframe))
        self.options['paramheight']['entries'][2*self.numparams+1].grid(row=self.i,column=1)
        self.options['paramlim']['entries'].append(Tk.Entry(self.optionsframe))
        self.options['paramlim']['entries'][2*self.numparams].grid(row=self.i,column=2)
        self.options['paramlim']['entries'].append(Tk.Entry(self.optionsframe))
        self.options['paramlim']['entries'][2*self.numparams+1].grid(row=self.i,column=3)
        self.numparams+=1

    def updateplot(self,*args):
        """ This will update the plot with the new plotting params, time period or
        if a different ISR measurement is used. This is all the command for the 
        update plot button and the pull down menus."""
        curvar = self.times['var'].get()
        if curvar=='None':
            return            
        if self.PC is None:
            return
        
        timestr = self.times['var'].get()
        self.times['var'].set(timestr)
        itime = int(float(timestr.split(' ')[0]))
        caststr =self.radarparam['var'].get()
        self.radarparam['var'].set(caststr)
        icase =  int(float(caststr.split(' ')[0]))
        self.getnewparams()
        (self.allhands,self.cbarsax)=self.PC.plotsingle(self.m,self.sp,self.fig,timenum=itime,icase=icase,cbarax=self.cbarsax)
        self.canvas.draw()
        
    def update(self):
        for field in self.options:
            self.options[field]['values']=[]
            for entry in self.options[field]['entries']:
                if field == 'reinterp':
                    self.options[field]['values'].append(self.options[field]['var'].get())
                else:
                    self.options[field]['values'].append(entry.get())
                    
    def EmptyFields(self):
        for field in self.options:
            if field!='reinterp':
                for entry in self.options[field]['entries']:
                    entry.delete(0,Tk.END)
        self.numparams=0
        for box in self.options['paramheight']['entries']:
            box.destroy()
        for box in self.options['paramlim']['entries']:
            box.destroy()
        self.options['paramheight']['entries']=[]
        self.options['paramlim']['entries']=[]
    
    
        
        
    def readindata(self):
        """ This will create the plot class object that will be used to plot all of the data. """
        if self.fn is None:
            return
        gpsloc = self.input['GPS']['entries'].get()
        ASloc = self.input['AllSky']['entries'].get()
        ISRloc=self.input['ISR']['entries'].get()
        
        self.PC = PlotClass(self.fn,GPSloc=gpsloc,ASloc=ASloc,ISRloc=ISRloc)
        self.m=self.PC.plotmap(self.fig,self.sp)
        (self.allhands,self.cbarsax)=self.PC.plotsingle(self.m,self.sp,self.fig,timenum=0,icase=0)
        self.canvas.draw()
        strlist = [insertinfo( str(j)+' $tmdy $thmsehms',posix=i[0],posixend=i[1]) for j, i in enumerate(self.PC.Regdict['Time'])]
        timearr = np.arange(len(strlist))
        paramar = np.zeros(len(strlist))
        self.times['list'] = timearr
        self.times['var'].set('')
        self.times['menu']['menu'].delete(0, 'end')
        for choice in strlist:
            self.times['menu']['menu'].add_command(label=choice, command=Tk._setit(self.times['var'], choice))
            self.times['var'].set(strlist[0])
        # deal with case with no isr data        
        if not self.PC.GDISR is None and   len(self.PC.params['paramheight'])>0:
            nparams = len(self.PC.params['paramheight'])
            paramar = np.arange(nparams)
            
            strlist2 = []
            for i, icase in enumerate(self.PC.params['paramheight']):
                stradd = str(i) +' '+ icase[0] +' at ' + str(int(icase[1])) +' km'
                strlist2.append(stradd)
            self.radarparam['list'] = paramar
            self.radarparam['var'].set('')
            self.radarparam['menu']['menu'].delete(0, 'end')
            for choice in strlist2:
                self.radarparam['menu']['menu'].add_command(label=choice, command=Tk._setit(self.radarparam['var'], choice))
            self.radarparam['var'].set(strlist2[0])
    
    def getnewparams(self):
        """ This will take all of the terms in the entries and update the param
            dictionary in the Plot Class object. If no Plot Class exists then 
            the fucntion will become a pass through."""        
        if self.PC is None:
            return
        paramtemp = {}
        for field in INIOPTIONS:
            if field == 'reinterp':
                paramtemp[field] = self.options[field]['var'].get()
            else:
                varval = [i.get() for i in self.options[field]['entries']]
               
                if field=='paramheight':
                    paramtemp[field] = [[varval[2*i],float(varval[2*i+1])] for i in np.arange(0,len(varval)/2) ]
                elif field=='paramlim':
                    paramtemp[field] = [[float(varval[2*i]),float(varval[2*i+1])] for i in np.arange(0,len(varval)/2) ]
                elif field=='timebounds':
                    splist = [i.split(' ') for i in varval ]
                    timelist = [splist[0][0],splist[0][1],splist[1][0],splist[1][1]]
                    paramtemp[field] = str2posix(timelist)
                elif len(varval)>1:
                    paramtemp[field]=[float(i) for i in varval]
                else:
                    paramtemp[field]=float(varval[0])
        self.PC.params=paramtemp  
                
    #%% Print and Save functions
    # XXX May use Plot class version
    def savefile(self):
        self.update()
        fn = fd.asksaveasfilename(title="Save File",filetypes=[('INI','.ini')])
        cfgfile = open(fn,'w')
        config = ConfigParser.ConfigParser()
        config.add_section('params')
        config.add_section('paramsnames')
        for field in self.options:
            if field=='reinterp':
                if self.options[field]['values'][0]:
                    config.set('params',field,'Yes')
                else:
                    config.set('params',field,'No')
            else:
                config.set('params',field," ".join(self.options[field]['values']))
                config.set('paramsnames',field,field)
        config.write(cfgfile)
        cfgfile.close()
    
    
                
                
    def loadfile(self,fn=None):
        self.EmptyFields()
        
        if (fn is None):
            fn = fd.askopenfilename(title="Load File",filetypes=[('INI','.ini')])
        params=readini(fn)
        
        self.fn = fn
        nparams=len(params['paramlim'])
        for n in range(nparams):
            self.AddParam()
        
        for field in params.keys():
            if not field in INIOPTIONS:
                continue
            data = params[field]

            if field=='paramheight':
                for n in range(nparams):
                    self.options['paramheight']['entries'][2*n].insert(0,str(data[n][0]))
                    self.options['paramheight']['entries'][2*n+1].insert(0,str(data[n][1]))

            elif field=='paramlim':
                for n in range(nparams):
                    self.options['paramlim']['entries'][2*n].insert(0,str(data[n][0]))
                    self.options['paramlim']['entries'][2*n+1].insert(0,str(data[n][1]))

            elif field=='timebounds':
                data = posix2str(data)
                self.options[field]['entries'][0].insert(0,data[0])
                self.options[field]['entries'][1].insert(0,data[1])

            elif not hasattr(data, "__len__"):
                if field=='reinterp':
                    if data:
                        self.options[field]['var'].set(1)
                    else:
                        self.options[field]['var'].set(0)
                else:
                    self.options[field]['entries'][0].insert(0,str(data))

            elif len(data)==2:
                self.options[field]['entries'][0].insert(0,str(data[0]))
                self.options[field]['entries'][1].insert(0,str(data[1]))
              
#%% Main function  the command line can be used to quickly populate the gui.
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
            -c Config file name.
            -i The directory that holds all of the TEC data in ionofile formats.
            -a The allsky data directory or GeoData file.
            -r The ISR data file.
            
             Example:
             python PlottingClass.py'''


    try:
        opts, args = getopt.gnu_getopt(argv,"ha:i:r:c:")
    except getopt.GetoptError:
        print(outstr)
        sys.exit(2)    
    
    gpsloc=''
    ASloc=''
    ISRloc=''
    plotdir=os.getcwd()
    inputdata = {'GPS':'','ISR':'','AllSky':''}
    inifile=None
    for opt, arg in opts:
        if opt == '-h':
            print(outstr)
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputdata['GPS'] = os.path.expanduser(arg)   
        elif opt in ("-a", "--asky"):
            inputdata['AllSky']=os.path.expanduser(arg)
        elif opt in ("-r", "--radar"):
            inputdata['ISR']=os.path.expanduser(arg)
        elif opt in ("-c", "--config"):
            inifile = os.path.expanduser(arg) 
        elif opt in ('-p','--pdir'):
            plotdir=os.path.expanduser(arg)
    root=Tk.Tk()
    app=App(root,inifile,inputdata)
    root.mainloop()
