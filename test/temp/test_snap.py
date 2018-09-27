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

fig, ax = plt.subplots()
ax.plot(x, y, "ro")

class MyApp():
    
    def __init__(self, coords, ax):
        
        #self.extent = extent
        self.coords = coords
        self.ax = ax
        points, = ax.plot(coords[:, 0], coords[:, 1], "ro")
        self.new_point, = ax.plot([], [], "gs")
        self.cid = self.ax.figure.canvas.mpl_connect("button_press_event", self.mouse_click)
        
    def mouse_click(self, event):
        if not event.inaxes: 
            return
        x = event.xdata
        y = event.ydata
        x, y = self.snap(x, y)
        self.new_point.set_xdata(x)
        self.new_point.set_ydata(y)
        event.canvas.draw()
      
    def snap(self, x, y):
        di = np.sqrt((x - self.coords[:, 0])**2 + (y - self.coords[:, 1])**2 )
        point = self.coords[np.argmin(di)]
        return point[0], point[1]  
        
myapp = MyApp(coords, ax)