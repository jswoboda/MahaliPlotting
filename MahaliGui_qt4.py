from PyQt4.uic import loadUiType
from PyQt4.QtGui import QFileDialog
from PyQt4.QtCore import QSettings

from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import(
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)

from PlottingClass import PlotClass
from GeoData.plotting import insertinfo

Ui_MainWindow,QMainWindow = loadUiType('MahaliGui.ui')

class Main(QMainWindow,Ui_MainWindow):
    def __init__(self,):
        super(Main,self).__init__()
        self.setupUi(self)
        
        self.actionLoad.triggered.connect(self.loadCfg)
        self.ReadButton.clicked.connect(self.readInData)
        self.UpdateButton.clicked.connect(self.updatePlot)
        self.TimeSlider.valueChanged.connect(self.updatePlot)
        self.TimeSlider.sliderMoved.connect(self.updateSlider)
        self.figs = []
        self.axs = []
        self.updating = False

        self.addmpl(Figure())

    def updateSlider(self):
        self.TimeDisplay.setText(self.strlist[self.TimeSlider.value()])

    def updatePlot(self):
        if not self.updating:
            self.updating = True
            self.rmmpl()
            self.addmpl(self.figs[self.TimeSlider.value()])
            self.updating = False

    def loadCfg(self):
        dlg = QFileDialog()
        dlg.setFilter("Config Files (*.ini)")
        self.inifn = str(dlg.getOpenFileName(self,'Open File','.','Config Files (*.ini)'))
        if(self.inifn):
            text = open(self.inifn,'r')
            self.ConfigBox.setText(text.read())
            text.close()

        #self.settings = QSettings(self.inifn, QSettings.IniFormat) #MAYBE USE QSETTINGS FOR INI READ/WRITE?!?!?

    def addmpl(self,fig):
        self.canvas = FigureCanvas(fig)
        self.mplvl.addWidget(self.canvas)
        self.canvas.draw()

        self.toolbar = NavigationToolbar(self.canvas,
                                         self.mplwindow,
                                         coordinates=True)
        self.addToolBar(self.toolbar)

    def rmmpl(self,):
        self.mplvl.removeWidget(self.canvas)
        self.canvas.close()
        self.mplvl.removeWidget(self.toolbar)
        self.toolbar.close()

    def readInData(self):
        f = open(self.inifn,'w')
        f.write(str(self.ConfigBox.toPlainText()))
        f.close()
        
        gpsloc = str(self.GPSBox.text())
        isrloc = str(self.ISRBox.text())
        asloc = str(self.AllSkyBox.text())
        self.PC = PlotClass(self.inifn,GPSloc=gpsloc,ASloc=asloc,ISRloc=isrloc)
        
        self.strlist = [insertinfo( str(j)+' $tmdy $thmsehms',posix=i[0],posixend=i[1]) for j, i in enumerate(self.PC.Regdict['Time'])]
        self.TimeSlider.setMaximum(len(self.strlist)-1)
        self.TimeSlider.setTracking(False)
        self.TimeSlider.setTickPosition(1)

        self.figs=[]
        self.axs=[]
        for t in range(len(self.strlist)):
            print(self.strlist[t])
            self.figs.append(Figure(figsize=(16,10)))
            self.axs.append(self.figs[t].add_subplot(111))
            m=self.PC.plotmap(self.figs[t],self.axs[t])
            (allhands,cbarsax)=self.PC.plotsingle(m,self.axs[t],self.figs[t],timenum=t,icase=0)
        self.rmmpl()
        self.addmpl(self.figs[0])


if __name__=='__main__':
    import sys
    from PyQt4 import QtGui
    import numpy as np
    
    app = QtGui.QApplication(sys.argv)
    main = Main()
    main.show()
    sys.exit(app.exec_())
