#!/usr/bin/env python
import matplotlib
import numpy as np
import pdb
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import Tkinter as Tk
import tkFileDialog as fd
import ConfigParser
from copy import copy
from PlottingClass import PlotClass
INIOPTIONS = ['latbounds','lonbounds','timebounds','timewin','asgamma','aslim','gpslim','paramlim','reinterp','paramheight','isrlatnum','isrlonnum','wl']


class App():

    def __init__(self,root):
        
        self.root=root
        self.root.title("Mahali")
        
        self.fn = None
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
        self.fig=Figure(figsize=(6, 5), dpi=100)
        self.sp=self.fig.add_subplot(111)
        self.sp.plot(np.arange(100))
        self.sp.set_title("Plot Title")
        self.sp.set_xlabel("X")
        self.sp.set_ylabel("Y")
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
            self.input[field]['entries']=[Tk.Entry(self.optionsframe,width=45)]
            self.input[field]['entries'][0].grid(row=irow+1,column=1,columnspan=2)
            self.input[field]['labels']=[Tk.Label(self.optionsframe,text=field)]
            self.input[field]['labels'][0].grid(row=irow+1,column=0)
        self.times={}
        self.times['var'] = Tk.StringVar(root,'None')
        self.times['options'] = ['None']
        self.times['menu'] = Tk.OptionMenu(self.optionsframe,self.times['var'],tuple(self.times['options']),command=self.updateplot)
        self.times['menu'].grid(row=len(self.input.keys())+1,column=1)
        self.times['label'] = Tk.Label(self.optionsframe,text='Choose Plot')
        self.times['label'].grid(row=len(self.input.keys())+1,column=0)
        
        self.i=len(self.input.keys())+2
        for field in self.options:
            self.options[field]['entries']=[]
            self.options[field]['labels']=[]
            if field=='reinterp':
                self.options[field]['var'] = a = Tk.IntVar()
                self.options[field]['entries'].append(Tk.Checkbutton(self.optionsframe,variable=a))
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

    def AddParam(self):
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
        curvar = self.times['var'].get()
        if curvar=='None':
            return
            
        
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
    
    def plotdata(self):
        """ """
    def readindata(self):
        if self.fn is None:
            return
        
        self.PC = PlotClass(self.fn,GPSloc=gpsloc,ASloc=ASloc,ISRloc=ISRloc)
        self.m=self.PC.plotmap(self.fig,self.sp)
        (self.allhands,self.cbarsax)=self.PC.plotsingle(self.m,self.sp,self.fig,timenum=0,icase=0,cbarax=self.cbarsax)
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
    
    def loadfile(self):
        self.EmptyFields()
        fn = fd.askopenfilename(title="Load File",filetypes=[('INI','.ini')])
        config = ConfigParser.ConfigParser()
        config.read(fn)
        self.fn = fn
        data = config.get('params','paramlim').split(" ")
        nparams=len(data)/2
        for n in range(nparams):
            self.AddParam()
        
        for field in config.options('params'):
            if not field in INIOPTIONS:
                continue
            data = config.get('params',field).split(" ")

            if field=='paramheight':
                for n in range(nparams):
                    self.options['paramheight']['entries'][2*n].insert(0,data[2*n])
                    self.options['paramheight']['entries'][2*n+1].insert(0,data[2*n+1])

            elif field=='paramlim':
                for n in range(nparams):
                    self.options['paramlim']['entries'][2*n].insert(0,data[2*n])
                    self.options['paramlim']['entries'][2*n+1].insert(0,data[2*n+1])

            elif field=='timebounds':
                self.options[field]['entries'][0].insert(0,data[0]+" "+data[1])
                self.options[field]['entries'][1].insert(0,data[2]+" "+data[3])

            elif len(data)==1:
                if field=='reinterp':
                    if data[0]=='Yes':
                        self.options[field]['var'].set(1)
                    else:
                        self.options[field]['var'].set(0)
                else:
                    self.options[field]['entries'][0].insert(0,data[0])

            if len(data)==2:
                self.options[field]['entries'][0].insert(0,data[0])
                self.options[field]['entries'][1].insert(0,data[1])
              
    
root=Tk.Tk()
app=App(root)
root.mainloop()
