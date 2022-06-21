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
from PyQt5.QtWidgets import QMenu
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
        self._create_actions()
        self._create_toolbar()
        self._create_menu()
        
        # Cargamos canales  
        self.load_channels(canales)
        
        # Modificamos tamaño inicial de MainWindoww
        self.resize(750, 550)
        self.statusBar = QStatusBar(self)
        self.setStatusBar(self.statusBar)
        
        # Tipos para knickpoints
        self.kp_types = {0 : {'ls':"", 'marker':"*", 'mec':"k", 'mew':0.5, 'mfc':"r", 'ms':10}, 
                         1 : {'ls':"", 'marker':"*", 'mec':"k", 'mew':0.5, 'mfc':"b", 'ms':10},
                         2 : {'ls':"", 'marker':"o", 'mec':"k", 'mew':0.5, 'mfc':"g", 'ms':6},
                         3 : {'ls':"", 'marker':"o", 'mec':"k", 'mew':0.5, 'mfc':"y", 'ms':6}}
        
        # Modos de dibujo (distintos tipos de perfiles)
        self._mode = 1 # 1. Long. profile, 2. Chi profile, 3. Area-slope profile, 4. ksn profile
        # pick_mode >> Modo de captura de puntos >> 0. No piking points 1. Knickpoint, 2. Regression, 3. Dam remover
        self._pick_mode = 0
        self.kp_type = 0
        self.current_regression = []
        
        # Dibujamos perfiles
        self.change_profile_graph(1)
    
    def load_channels(self, channels):
        """
        Loads channels into App. 
        
        channels : np.array of topopy.Channel objects
        """
        self._channels = channels
        self._nchannels = len(channels)
        self._active = 0
        
    def _create_actions(self):
        """
        Función interna que crea QActions (botones y acciones para menu y toolbar)
        """
        # Desplazamiento entre perfiles (no ncesaria referencia, los creamos y los conectamos con funciones)
        self.tb_button_prev = QAction(QIcon("icons/arrow-180.png"), "Previous", self)
        self.tb_button_next = QAction(QIcon("icons/arrow-000.png"), "Next", self)
        
        # Tipos de gráficos
        self.tb_button_longP = QAction(QIcon("icons/long_prof.ico"), "Longitudinal profile", self)
        self.tb_button_chiP = QAction(QIcon("icons/chi_prof.ico"), "Chi profile", self)
        self.tb_button_ASP = QAction(QIcon("icons/loglog_prof.ico"), "log(a)-log(s) profile", self)
        self.tb_button_ksnP = QAction(QIcon("icons/ksn_prof.ico"), "ksn profile", self)
        
        # Modos de captura de puntos (necesaria referencia para poder cambiarles el estado)
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
        
        # Acciones del menú File
        self.menu_load = QAction("Load Channels", self) 
        self.menu_save = QAction("Save Channels", self)
        self.menu_add = QAction("Add Channel", self)
        self.menu_remove = QAction("Remove Channel", self)
        
        
        
        # ===============================================================================
        # Conectamos QActions y Widgets con funciones
        # Conectamos los 4 botones de tipos de perfiles a la función change_profile_graph
        self.tb_button_longP.triggered.connect(lambda x: self.change_profile_graph(1))
        self.tb_button_chiP.triggered.connect(lambda x: self.change_profile_graph(2))
        self.tb_button_ASP.triggered.connect(lambda x: self.change_profile_graph(3))
        self.tb_button_ksnP.triggered.connect(lambda x: self.change_profile_graph(4))
        
        # Conectamos botones de desplazamiento a función
        self.tb_button_prev.triggered.connect(lambda x: self.next_profile(-1))
        self.tb_button_next.triggered.connect(lambda x: self.next_profile(1))
        
        # Conectamos boton de refresh a función
        self.tb_button_refresh.triggered.connect(self.calculate_gradients)
        
        # Conectamos boton de KP a funcion
        self.tb_button_KP.triggered.connect(self.button_KP_click)
        
        # Conectamos boton de Regression a funcion
        self.tb_button_reg.triggered.connect(self.button_reg_click)
    
        
    def _create_menu(self):
        """
        Función interna para crear el menú de la aplicación
        """
        # Creamos barra de menú vacía
        menubar = self.menuBar()
        filemenu = menubar.addMenu("&File")
        editmenu =menubar.addMenu("&Edit")
        
        # Añadimos opciones
        filemenu.addAction(self.menu_load)
        filemenu.addAction(self.menu_save)
        filemenu.addAction(self.menu_add)
        filemenu.addAction(self.menu_remove)
        
        
        
        
        
    
    def _create_toolbar(self):
        """
        Función interna para crear la barra de herramientas principal de Aplicación
        """
        # Creamos barra de herramientas
        toolbar = QToolBar("Action toolbar")
        toolbar.setIconSize(QSize(20,20))
        self.addToolBar(toolbar)
           
        # Añadimos botones a toolbar
        toolbar.addAction(self.tb_button_prev)
        toolbar.addAction(self.tb_button_next)
        toolbar.addAction(self.tb_button_longP)
        toolbar.addAction(self.tb_button_chiP)
        toolbar.addAction(self.tb_button_ASP)
        toolbar.addAction(self.tb_button_ksnP)
        toolbar.addAction(self.tb_button_KP)
        toolbar.addAction(self.tb_button_reg)
        toolbar.addAction(self.tb_button_dam)
        toolbar.addWidget(self.npointLabel)
        toolbar.addWidget(self.npointSpinBox)
        toolbar.addAction(self.tb_button_refresh) 
        toolbar.setStyleSheet("QToolBar{spacing:2px;}")
        

    def button_reg_click(self):
        # Handler para botón de regressions
        
        if self.tb_button_reg.isChecked():
            self._pick_mode = 2 # Regression selection on
            self.statusBar.showMessage("Regression mode ON")
            self.pc_id = self.canvas.mpl_connect("pick_event", self.pick_point)
            self.tb_button_KP.setEnabled(False)
            self.tb_button_dam.setEnabled(False)
            self.current_regression = []

        else:
            self.statusBar.clearMessage()
            self.current_regression = []
            self.canvas.mpl_disconnect(self.pc_id)
            self.tb_button_KP.setEnabled(True)
            if self._mode == 2 or self._mode == 1:
                self.tb_button_dam.setEnabled(True)
    
    
    def button_KP_click(self):
        # Handler para botón de knickpoint        
        if self.tb_button_KP.isChecked():
            self._pick_mode = 1 # Knickpoint selection on
            self.statusBar.showMessage("Knickpoint selection ON - Knickpoint type: {}".format(self.kp_type))
            self.pc_id = self.canvas.mpl_connect("pick_event", self.pick_point)
            self.tb_button_reg.setEnabled(False)
            self.tb_button_dam.setEnabled(False)

        else:
            self.canvas.mpl_disconnect(self.pc_id)
            self.statusBar.clearMessage()
            if self._mode == 2:
                self.tb_button_reg.setEnabled(True)
                self.tb_button_dam.setEnabled(True)
            elif self._mode == 1:
                self.tb_button_dam.setEnabled(True)


    def pick_point(self, event):
        """
        Pick a point in the active profile 
        
        event : matplotlib picker event
        """
        # Check if App has channels (and take the active one)
        if self._nchannels > 0:
            canal = self._channels[self._active]
        else:
            return
        
        # In the case that many points are picked (true if the profile has several points). Take the middle one.
        if len(event.ind) > 2:
            ind = (event.ind[-1] + event.ind[0]) // 2
        else:
            ind = event.ind[0]
            
        # If self._pick_mode == 1 --> Selecting Knickpoints
        if self._pick_mode == 1:
            # Left button >> Add knickpoint
            if event.mouseevent.button==1:
                canal.add_kp(ind, self.kp_type)
            # Rigth button >> Remove knickpoint
            elif event.mouseevent.button==3:
                kps = canal._kp[:, 0]
                diffs = np.abs(kps - ind)
                min_kp = np.min(diffs)
                pos_min = np.argmin(diffs)
                if min_kp < 3:
                    ind_to_remove = kps[pos_min]
                    canal.remove_kp(ind_to_remove)
            # Middle button >> Change knickpoint type
            elif event.mouseevent.button==2:
                self.kp_type += 1
                self.kp_type = self.kp_type % 4
                self.statusBar.showMessage("Knickpoint selection ON - Knickpoint type: {}".format(self.kp_type))
                
        # If self._pick_mode == 2 --> Selecting regression
        elif self._pick_mode == 2:
            # Left button >> Add regression point
            if event.mouseevent.button==1:
                if len(self.current_regression) == 0:
                    # Si no hay ningun punto introducido
                    # Introducimos el punto, lo dibujamos y salimos sin llamar a self._draw()
                    self.current_regression.append(ind)
                    self.ax.plot(event.mouseevent.xdata, event.mouseevent.ydata, ls="", marker="+", ms=10)
                    self.canvas.draw()
                    return
                elif len(self.current_regression) == 1:
                    # Si hay un punto introducido, añadimos la regresión
                    self.current_regression.append(ind)
                    canal.add_regression(self.current_regression[0], self.current_regression[1])
                    self.current_regression = []
            
            # Right button >> Remove regression
            if event.mouseevent.button==3:
                canal.remove_regression(ind)
                self.current_regression = []

        self._draw()

    
    def calculate_gradients(self):
        # Handler to tb_button_refresh
        # Recalcuates gradients of the active channel with specific number of points (npointsSpinBox)
        # Get the active profile
        canal = self._channels[self._active]
        npoints = self.npointSpinBox.value()
        # If mode == 3, recalculate slope
        if self._mode == 3:
            canal.calculate_gradients(npoints, 'slp')
        # If mode == 4, recalculate ksn
        elif self._mode == 4:
            canal.calculate_gradients(npoints, 'ksn')
        
        self._draw()
        
    def next_profile(self, direction):
        # Handler to tb_button_prev and tb_button_next buttons 
        # Select the next / previous channel of the channel list (self._channels)
        self._active += direction
        self._active = self._active % self._nchannels
        self.current_regression = []
        
        if self._mode == 3:    
            npoints = self._channels[self._active]._slp_np
            self._enable_spinBox(npoints)
        elif self._mode == 4:
            npoints = self._channels[self._active]._ksn_np
            self._enable_spinBox(npoints)
        
        self._draw()
        
    def change_profile_graph(self, graph_number):
        # Changes the profile graph type and activate-deactivate specific buttons and options
        self._mode = graph_number
        self.current_regression = []
        
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
            
        if self.tb_button_KP.isChecked():
            self.tb_button_reg.setEnabled(False)
            self.tb_button_dam.setEnabled(False)
            
        if self.tb_button_reg.isChecked():
            self.tb_button_KP.setEnabled(False)
            self.tb_button_dam.setEnabled(False)
            
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
        # Function to draw the active channel in the current graphic mode (self._mode)
        
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
            self.ax.plot(di, zi, c="k", lw=1.25, picker=True,  pickradius=5)
            
            # Draw knickpoints
            if len(canal._kp) > 0:
                for k in canal._kp:
                    self.ax.plot(di[k[0]], zi[k[0]], **self.kp_types[k[1]])
            
            
        elif self._mode == 2: # Chi profile
            title = "Chi profile ($\Theta$=0.45)"
            if canal.get_name():
                title += "  [{}]".format(canal.get_name())
            self.ax.set_title(title)
            self.ax.set_xlabel("$\chi$ (m)")
            self.ax.set_ylabel("Elevation (m)")  
            chi = canal.get_chi()
            zi = canal.get_z()
            self.ax.plot(chi, zi, c="k", lw=1.25, picker=True, pickradius=5)
            
            # Draw knickpoints
            if len(canal._kp) > 0:
                for k in canal._kp:
                    self.ax.plot(chi[k[0]], zi[k[0]], **self.kp_types[k[1]])
                    
            # Draw regressions
            if len(canal._regressions) > 0:
                # Channel regressions are tuples (p1, p2, poli, R2)
                # p1, p2 >> positions of first and second point of the regression
                # poli >> Polinomial with the regression
                # R2 >> Determination coeficient of the regression
                for reg in canal._regressions:
                    chi1 = chi[reg[0]]
                    chi2 = chi[reg[1]]
                    poli = reg[2]
                    z1 = np.polyval(poli, chi1)
                    z2 = np.polyval(poli, chi2)
                    
                    self.ax.plot([chi1, chi2], [z1, z2], c="r", ls="--", lw=1.5)
    
        
        elif self._mode == 3: # Area-slope profile
            title = "Area-slope profile"
            if canal.get_name():
                title += "  [{}]".format(canal.get_name())
            self.ax.set_title(title)
            self.ax.set_xlabel("Area")
            self.ax.set_ylabel("Slope")
            ai= canal.get_a(cells=False)
            slp = canal.get_slope()
            self.ax.plot(ai, slp, marker=".", ls = "", color="k", mfc="k", ms=5, picker=True, pickradius=5)
            self.ax.set_xscale("log")
            self.ax.set_yscale("log")
            
            # Draw knickpoints
            if len(canal._kp) > 0:
                for k in canal._kp:
                    self.ax.plot(ai[k[0]], slp[k[0]], **self.kp_types[k[1]])
            
            
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
            self.ax.plot(di, ksn, c="k", lw=1.25, picker=True, pickradius=5)
            
            # Draw knickpoints
            if len(canal._kp) > 0:
                for k in canal._kp:
                    self.ax.plot(di[k[0]], ksn[k[0]], **self.kp_types[k[1]])
                    
            # Draw regressions
            if len(canal._regressions) > 0:
                # Channel regressions are tuples (p1, p2, poli, R2)
                # p1, p2 >> positions of first and second point of the regression
                # poli >> Polinomial with the regression
                # R2 >> Determination coeficient of the regression
                for reg in canal._regressions:
                    ksn = reg[2][0]
                    d1 = di[reg[0]]
                    d2 = di[reg[1]]
                    self.ax.plot([d1, d2], [ksn, ksn], c="r", ls="--", lw=1.5)
                    
                    
            
        self.canvas.draw()

# 3. Creamos aplicación        
app = QApplication(sys.argv)

canales = np.load("canales.npy", allow_pickle=True)

# 4. Creamos ventana con la nueva clase y la mostramos
win = MainWindow(canales)
win.show()

# 5. Iniciamos la aplicación
app.exec_()