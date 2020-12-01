import sys

from PyQt5 import QtCore, QtWidgets
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import numpy as np
import matplotlib
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure



class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.ax = fig.add_subplot(111)
        super(MplCanvas, self).__init__(fig)
        
        self.ax.set_xlim((0,100))
        self.ax.set_ylim((0,100))
        self.ax.xaxis.set_major_locator(ticker.MultipleLocator(20))
        self.ax.xaxis.set_minor_locator(ticker.MultipleLocator(10))
        self.ax.yaxis.set_major_locator(ticker.MultipleLocator(20))
        self.ax.yaxis.set_minor_locator(ticker.MultipleLocator(10))
        self.ax.set_title("Use the left button to create the polygon,\n and the right one to close it")
        self.ax.grid(True, which='major', axis="both", linestyle ="-", c="0.8", lw=0.5)
        self.ax.grid(True, which='minor', axis="both", linestyle ="--", c="0.4", lw=0.25)

        self.xdat = []
        self.ydat = []
        self.line, = self.ax.plot([], [], ls="-", c="r", marker="s", mec="k", mfc="b")
        self.cid1 = self.ax.figure.canvas.mpl_connect("button_press_event", self.mouse_click)
    
    def mouse_click(self, event):
        if not event.inaxes: 
            return
        
        if event.button == 1:
            self.xdat.append(event.xdata)
            self.ydat.append(event.ydata)
            self.line.set_xdata(self.xdat)
            self.line.set_ydata(self.ydat)        
            event.canvas.draw()
        
        elif event.button == 3 and len(self.xdat) > 2:
            data = np.array((self.xdat, self.ydat)).T
            color = [np.random.random() for n in range(3)]
            poligon = mpatches.Polygon(data, color=color, alpha=0.6)
            event.inaxes.add_patch(poligon)
            self.xdat = []
            self.ydat = []
            self.line.set_data(self.xdat, self.ydat)
            event.canvas.draw()
                


class MainWindow(QtWidgets.QMainWindow):

    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)

        # Create the maptlotlib FigureCanvas object, 
        # which defines a single set of axes as self.axes.
        sc = MplCanvas(self, width=5, height=4, dpi=100)
        #sc.axes.plot([0,1,2,3,4], [10,1,20,3,40])
        xi = np.linspace(0, 2*np.pi)
        yi = np.sin(xi)
        sc.ax.plot(xi, yi)
        self.setCentralWidget(sc)

        self.show()


app = QtWidgets.QApplication(sys.argv)
w = MainWindow()
app.exec_()