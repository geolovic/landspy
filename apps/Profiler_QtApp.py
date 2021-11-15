#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ejemplo de utilizacion de toolbar en Qt5

Para crear una toolbar en nuestra MainWindow y añadirla, simplemente utilizaremos
    toolbar = QToolBar("My main toolbar")
    self.addToolBar(toolbar)
        
"""


# 1. Importamos módulos necesarios
from PyQt5.QtWidgets import QApplication, QMainWindow, QToolBar, QLabel, QAction, QStatusBar, QSpinBox
from PyQt5.QtWidgets import QLineEdit
from PyQt5.QtCore import Qt, QSize
from PyQt5.QtGui import QIcon
import sys
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.ticker as ticker
import numpy as np
from topopy import Channel


# 2. Creamos subclase de MainWindow

class MainWindow(QMainWindow):
       
    def __init__(self, canales=None):
        
        # Iniciamos objeto MainWindow
        super().__init__()
        self.setWindowTitle("Barra de herramientas")

        # Creamos Figure canvas como widget central
        self.fig = Figure()
        self.canvas = FigureCanvas(self.fig)
        self.ax = self.fig.add_subplot(111)
        self.ax.set_xlim((0, 100))
        self.ax.set_ylim((0, 100))

        # Hacemos que el Figure Canvas sea el widget central
        self.setCentralWidget(self.canvas)
        self._create_toolbar()
             
        self.load_channels(canales)
        
        # Modificamos tamaño inicial de MainWindoww
        self.resize(750, 550)
        self.statusBar = QStatusBar(self)
        self.setStatusBar(self.statusBar)
        
        # Creamos modos de dibujo (distintos perfiles)
        self._mode = 1 # 1. Long. profile, 2. Chi profile, 3. Area-slope profile, 4. ksn profile
        # Y picking mode
        self._pick_mode = 0 # 0. No piking points 1. Knickpoint, 2. Regression, 3. Dam remover
        
        self._draw()
    
    def load_channels(self, channels):
        self._channels = channels
        self._nchannels = len(channels)
        self._active = 0
    
    def _create_toolbar(self):
        # Creamos barra de herramientas
        toolbar = QToolBar("Action toolbar")
        toolbar.setIconSize(QSize(20,20))
        self.addToolBar(toolbar)
        
        # Creamos botones
        # Desplazamiento entre perfiles
        tb_button_prev = QAction(QIcon("icons/arrow-180.png"), "Previous", self)
        tb_button_next = QAction(QIcon("icons/arrow-000.png"), "Next", self)
        # Tipos de gráficos
        tb_button_longP = QAction(QIcon("icons/long_prof.ico"), "Longitudinal profile", self)
        tb_button_chiP = QAction(QIcon("icons/chi_prof.ico"), "Chi profile", self)
        tb_button_ASP = QAction(QIcon("icons/loglog_prof.ico"), "log(a)-log(s) profile", self)
        tb_button_ksnP = QAction(QIcon("icons/ksn_prof.ico"), "ksn profile", self)
        # Modos de captura de puntos
        self.tb_button_KP = QAction(QIcon("icons/flag.ico"), "Set Knickpoint", self)
        self.tb_button_reg = QAction(QIcon("icons/reg.ico"), "Create regression", self)
        self.tb_button_dam = QAction(QIcon("icons/dam.ico"), "Remove dam", self)        
        self.tb_button_KP.setCheckable(True)
        self.tb_button_reg.setCheckable(True)
        self.tb_button_dam.setCheckable(True)
        # Número de puntos para ksn y slope
        self.npointLabel = QLabel("N Points:")
        self.npointSpinBox = QSpinBox()
        self.tb_button_refresh = QAction(QIcon("icons/arrow-circle.png"), "Refresh", self)
        self.npointSpinBox.setFocusPolicy(Qt.NoFocus)
    
        # Los añadimos a toolbar
        toolbar.addAction(tb_button_prev)
        toolbar.addAction(tb_button_next)
        toolbar.addAction(tb_button_longP)
        toolbar.addAction(tb_button_chiP)
        toolbar.addAction(tb_button_ASP)
        toolbar.addAction(tb_button_ksnP)
        toolbar.addAction(self.tb_button_KP)
        toolbar.addAction(self.tb_button_reg)
        toolbar.addAction(self.tb_button_dam)
        toolbar.addWidget(self.npointLabel)
        toolbar.addWidget(self.npointSpinBox)
        toolbar.addAction(self.tb_button_refresh) 
        toolbar.setStyleSheet("QToolBar{spacing:2px;}")
        
        # Conectamos los 4 botones de tipos de perfiles a la función change_profile_graph
        tb_button_longP.triggered.connect(lambda x: self.change_profile_graph(1))
        tb_button_chiP.triggered.connect(lambda x: self.change_profile_graph(2))
        tb_button_ASP.triggered.connect(lambda x: self.change_profile_graph(3))
        tb_button_ksnP.triggered.connect(lambda x: self.change_profile_graph(4))
        
        # Conectamos botones de desplazamiento a función
        tb_button_prev.triggered.connect(lambda x: self.next_profile(-1))
        tb_button_next.triggered.connect(lambda x: self.next_profile(1))
        
        # Conectamos boton de refresh a función
        self.tb_button_refresh.triggered.connect(self.calculate_gradients)
        
        # Conectamos boton de KP a funcion
        self.tb_button_KP.triggered.connect(self.button_KP_click)
    
    def button_KP_click(self):
        self.statusBar.showMessage("Click en boton knickpoint")
        self.canvas.mpl_connect("pick_event", self.pick_point)
        # if self.tb_button_KP.isChecked():
            
        #     self._pick_mode = 1
        #     self.pc_id = self.canvas.mpl_connect("mouse_press_event", self.pick_point)
        #     self.statusBar.showMessage("Knickpoint mode ON", self.pc_id)
            
        # else:
        #     self.statusBar.showMessage("Knickpoint mode OFF")
        #     self.canvas.mpl_disconnect(self.pc_id)
    
    def pick_point(self, event):
        """
        Pick a point in the active profile 
        :param event: matplotlib picker event
        """
        # In the case that many points are picked (true if the profile has several points). Take the center one.
        
        # self.ax.plot(event.xdata, event.ydata, mew=0.5, mec="k", ms=10)
        # self.canvas.draw()
        # self.canvas.draw()
        if len(event.ind) > 2:
            ind = (event.ind[-1] + event.ind[0]) // 2
        else:
            ind = event.ind[0]

        self.statusBar.showMessage("Click en {}".format(ind))
        # self.statusBar.showMessage("Clicked on point {}".format(ind))
        # self._channels[self._active]._knickpoints.append(ind)
        # self._draw()
        
    
    def calculate_gradients(self):
        # Get the active profile
        canal = self._channels[self._active]
        npoints = self.npointSpinBox.value()
        # If mode == 3, calculate slope
        if self._mode == 3:
            canal.calculate_gradients(npoints, 'slp')
        # If mode == 4, calculate ksn
        elif self._mode == 4:
            canal.calculate_gradients(npoints, 'ksn')
        
        self._draw()
        
    def next_profile(self, direction):
        self._active += direction
        self._active = self._active % self._nchannels
        
        if self._mode == 3:    
            npoints = self._channels[self._active]._slp_np
            self._enable_spinBox(npoints)
        elif self._mode == 4:
            npoints = self._channels[self._active]._ksn_np
            self._enable_spinBox(npoints)
        
        self._draw()
        
    def change_profile_graph(self, graph_number):
        self._mode = graph_number
        #self.statusBar.showMessage("Has seleccionado el gráfico número " + str(graph_number))
        
        if self._mode == 1: # Longitudinal profile
            self.tb_button_dam.setEnabled(True)
            self.tb_button_reg.setEnabled(False)
            self._disable_spinBox()
            self._draw()
        
        elif self._mode == 2: # Chi profile
            self.tb_button_dam.setEnabled(True)
            self.tb_button_reg.setEnabled(True)
            self._disable_spinBox()
            self._draw()

        elif self._mode == 3: # Area-slope profile
            self.tb_button_dam.setEnabled(False)
            self.tb_button_reg.setEnabled(False)
            npoints = self._channels[self._active]._slp_np
            self._enable_spinBox(npoints) 
            self._draw()
            
        elif self._mode == 4: # Ksn profile
            self.tb_button_dam.setEnabled(False)
            self.tb_button_reg.setEnabled(False)           
            npoints = self._channels[self._active]._ksn_np
            self._enable_spinBox(npoints) 
            self._draw()
            
    def _enable_spinBox(self, value):
        self.npointSpinBox.setValue(value) 
        self.npointSpinBox.setEnabled(True)
        self.npointLabel.setEnabled(True)
        self.tb_button_refresh.setEnabled(True)
        
    def _disable_spinBox(self):
        self.npointSpinBox.setValue(0)
        self.npointSpinBox.findChild(QLineEdit).deselect()
        self.npointSpinBox.setEnabled(False)
        self.npointLabel.setEnabled(False)
        self.tb_button_refresh.setEnabled(False)
        
    def _draw(self):
        
        if self._nchannels > 0:
            canal = self._channels[self._active]
        else:
            return

        self.ax.clear()
        if self._mode == 1: # Longitudinal profile
            title = "Longitudinal profile"
            if canal.get_name():
                title += "  [{}]".format(canal.get_name())
            self.ax.set_title(title)
            self.ax.set_xlabel("Distance (m)")
            self.ax.set_ylabel("Elevation (m)")
            di = canal.get_d()
            zi = canal.get_z()
            self.ax.plot(di, zi, c="b", lw=1.25, picker=True,  pickradius=5)
            
            # # Draw knickpoints
            # if len(canal._knickpoints) > 0:
            #     for k in canal._knickpoints:
            #         self.ax.plot(di[k[0]], zi[k[0]], mew=0.5, mec="k", ms=10)
            
            
        elif self._mode == 2: # Chi profile
            title = "Chi profile ($\Theta$=0.45)"
            if canal.get_name():
                title += "  [{}]".format(canal.get_name())
            self.ax.set_title(title)
            self.ax.set_xlabel("$\chi$ (m)")
            self.ax.set_ylabel("Elevation (m)")  
            chi = canal.get_chi()
            zi = canal.get_z()
            self.ax.plot(chi, zi, c="b", lw=1.25, picker=True, pickradius=5)
        
        elif self._mode == 3: # Area-slope profile
            title = "Area-slope profile"
            if canal.get_name():
                title += "  [{}]".format(canal.get_name())
            self.ax.set_title(title)
            self.ax.set_xlabel("Area")
            self.ax.set_ylabel("Slope")
            ai= canal.get_a(cells=False)
            slp = canal.get_slope()
            self.ax.plot(ai, slp, "r.", picker=True, pickradius=5)
            self.ax.set_xscale("log")
            self.ax.set_yscale("log")
            
        elif self._mode == 4: # ksn profile
            title = "Ksn profile"
            if canal.get_name():
                title += "  [{}]".format(canal.get_name())
            self.ax.set_title(title)
            self.ax.set_title(title)
            self.ax.set_xlabel("Distance (m)")
            self.ax.set_ylabel("ksn")
            di = canal.get_d()
            ksn = canal.get_ksn()
            self.ax.plot(di, ksn, c="b", lw=1.25, picker=True, pickradius=5)
        
        self.canvas.draw()

# 3. Creamos aplicación        
app = QApplication(sys.argv)

canales = np.load("canales2.npy", allow_pickle=True)

# 4. Creamos ventana con la nueva clase y la mostramos
win = MainWindow(canales)
win.show()

# 5. Iniciamos la aplicación
app.exec_()