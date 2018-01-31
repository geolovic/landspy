import numpy as np
import gdal
import sys
sys.path.append("../../")
from topopy import DEM, Grid, Flow
import matplotlib.pyplot as plt

BOTTOM = 1
HEAD = 2

# INITIAL PARAMETERS
DIRECTION = HEAD
threshold = 250

# Get the flow object
dem = DEM("../data/small25.tif")
flow = Flow()
flow.load_gtiff("../data/small25_fd.tif")

class GraphApp():
	def __init__(self, ax, dem, flow):
		self.ax = ax
		self.dem = dem.read_array()
		self.dem = self.dem.astype(np.float)
		self.dem[dem.get_nodata_pos()] = np.nan
		self.procc = np.empty(dem.get_dims())
		self.procc.fill(np.nan)
		fac = flow.flow_accumulation(False).read_array()
		w = fac < 250
		self.dims = flow._dims
		self.w = np.array(w, dtype=np.float)
		self.w[w] = np.nan
		self.cid = ax.figure.canvas.mpl_connect('button_press_event', self.bpress)

		w = fac > 250
		w = w.ravel()
		I  = w[flow._ix]
		self.ix  = flow._ix[I]
		self.ixc = flow._ixc[I]
		self.draw()
		self.ind = 0
		
	def draw(self):
		ax.imshow(self.dem)
		#ax.imshow(self.w)
		ax.imshow(self.procc, cmap=plt.cm.Wistia)
		ax.figure.canvas.draw()

	def process(self):
		ind = self.ix[self.ind]
		row, col = np.unravel_index(ind, self.dims)
		self.procc[row, col] = 1
		self.ind += 1
	
	def bpress(self, event):
		if event.button == 1:
			self.process()
		if event.button == 3:
			for n in range(20):
				self.process()
		self.draw()

fig, ax = plt.subplots()
mygraph = GraphApp(ax, dem, flow)

plt.show()

