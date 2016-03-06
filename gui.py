import matplotlib
import numpy as np
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import sys
if sys.version_info[0] < 3:
    import Tkinter as Tk
else:
    import tkinter as Tk

"""
f = Figure(figsize=(5, 4), dpi=100)
a = f.add_subplot(111)

a.plot()
canvas = FigureCanvasTkAgg(f, master=root)
canvas.show()
canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

toolbar = NavigationToolbar2TkAgg(canvas, root)
toolbar.update()
canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
"""

class App():

    def __init__(self,root):
        self.root=root
        self.root.title("Mahali")

        self.titleframe=Tk.Frame(self.root)
        self.titleframe.grid(row=0,columnspan=2)
        self.leb=Tk.Label(self.titleframe,text="Mahali GUI",font=("Helvetica", 16))
        self.leb.grid(row=0,columnspan=3)
        
        self.plotframe=Tk.Frame(self.root)
        self.plotframe.grid(row=1,column=0)
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

        self.optionsframe=Tk.LabelFrame(self.root,text="options",padx=5,pady=5)
        self.optionsframe.grid(row=1,column=1)

        self.o1 = Tk.Entry(self.optionsframe)
        self.o1.grid(row=1,column=1)
        self.o1label = Tk.Label(self.optionsframe,text="option 1")
        self.o1label.grid(row=1,column=2)

        self.o2 = Tk.Entry(self.optionsframe)
        self.o2.grid(row=2,column=1)
        self.o2label = Tk.Label(self.optionsframe,text="option 2")
        self.o2label.grid(row=2,column=2)

        self.o3 = Tk.Entry(self.optionsframe)
        self.o3.grid(row=3,column=1)
        self.o3label = Tk.Label(self.optionsframe,text="option 2")
        self.o3label.grid(row=3,column=2)

root=Tk.Tk()
app=App(root)
root.mainloop()
