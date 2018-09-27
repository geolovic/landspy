# -*- coding: utf-8 -*-

from topopy import DEM, Flow
import numpy as np
import ogr, gdal
import matplotlib.pyplot as plt


# Get DEM, Flow Direction and Flow Accumulation

dem = DEM("../data/in/jebja30.tif")
fd = Flow(dem)
thr = 1000
fac = fd.get_flow_accumulation()

row, col = np.where(np.logical_and(fac.read_array() > thr, fac.read_array() != fac.get_nodata()))
x, y = fd.cell_2_xy(row, col)
coords = np.array((x, y)).T

geot = fac.get_geotransform()
xmin = geot[0]
xmax = geot[0] + fac.get_size()[0] * geot[1]
ymax = geot[3]
ymin = geot[3] + fac.get_size()[1] * geot[5]
extent = (xmin, ymin, xmax, ymax)
fig, ax = plt.subplots()
ax.plot(x, y, "ro")

class MyApp():
    
    def __init__(self, extent, coords, ax):
        
        self.extent = extent
        self.coords = coords
        self.ax = ax
        points, = ax.plot(coords[:, 0], coords[:, 1], "ro")
        self.points, = ax.plot([], [], "gs")
        self.rpoints, = ax.plot([], [], "bs")
        self.cid = self.ax.figure.canvas.mpl_connect("key_press_event", self.mouse_click)
        
    def mouse_click(self, event):
        # Generamos 5 puntos al azar dentro del raster
        xi = np.random.randint(self.extent[0], self.extent[2], 100)
        yi = np.random.randint(self.extent[1], self.extent[3], 100)
        sxi, syi = self.snap(xi, yi)
        self.points.set_xdata(xi)
        self.points.set_ydata(yi)
        self.rpoints.set_xdata(sxi)
        self.rpoints.set_ydata(syi)
        event.canvas.draw()
      
    def snap(self, xi, yi):
        aux_xi = xi.reshape((xi.size, 1))
        aux_yi = yi.reshape((yi.size, 1))
        xci = self.coords[:, 0].reshape((1, self.coords.shape[0]))
        yci = self.coords[:, 1].reshape((1, self.coords.shape[0]))
        di = np.sqrt((aux_xi - xci)**2 + (aux_yi - yci)**2 )
        
        
        pos = np.argmin(di, axis=1)
        
        return coords[pos, 0], coords[pos, 1]
        
myapp = MyApp(extent, coords, ax)