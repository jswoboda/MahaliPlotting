#!/usr/bin/env python
from sys import stderr
import PyQt5
from PyQt5.uic import loadUiType
from PyQt5.QtWidgets import QFileDialog

from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import(
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)

from MahaliPlotting import PlotClass,Path
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
        dlg.setNameFilters(["Config Files (*.ini)"])
        self.inifn = Path(dlg.getOpenFileName(self,'Open File','.','Config Files (*.ini)')[0])
        if self.inifn.is_file():
            text = self.inifn.open('r')
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
        try:
            self.inifn
        except AttributeError:
            print('you must first load an .ini file',file=stderr)
            return

        f = self.inifn.open('w')
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

    app = PyQt5.QtWidgets.QApplication(sys.argv)
    main = Main()
    main.show()
    sys.exit(app.exec_())
