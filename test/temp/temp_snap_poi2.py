# -*- coding: utf-8 -*-
"""
Editor de Spyder

Este es un archivo temporal
"""

from topopy import Flow
import numpy as np
import matplotlib.pyplot as plt



def snap_points(self, input_points, threshold, kind="channel", remove_duplicates=False):
    """
    Snap input points to channel cells or to stream POI
    
    Parameters:
    ===========
    array : *numpy.ndarray*
      Numpy 2-D ndarray, which first two columns are x and y coordinates [x, y, ...]
    threshold : *int*
      Flow accumulation threshold (in number of cells) to extract channel cells or stream POI 
    kind : *str* {'channel', 'heads', 'confluences', 'outlets'}  
        Kind of point to snap input points
    
    Returns:
    ===========
    numpy.ndarray
      Numpy ndarray with two columns [xi, yi] with the snap points
    """
    
    # Extract a numpy array with the coordinate to snap the points
    if kind in ['heads', 'confluences', 'outlets']:
        poi = self.get_stream_poi(threshold, kind, "XY")         
    else:
        fac = self.get_flow_accumulation(nodata=False, asgrid=False)
        row, col = np.where(fac >= threshold)
        x, y = self.cell_2_xy(row, col)
        poi = np.array((x, y)).T
          
    # Get array reshaped for the calculation
    xi = input_points[:, 0].reshape((input_points.shape[0], 1))
    yi = input_points[:, 1].reshape((input_points.shape[0], 1))
    xci = poi[:, 0].reshape((1, poi.shape[0]))
    yci = poi[:, 1].reshape((1, poi.shape[0]))
    
    # Calculate distances and get minimum
    di = np.sqrt((xi - xci)**2 + (yi - yci)**2 )
    pos = np.argmin(di, axis=1)
    
    # Get the rest of columns (if the case)
    if input_points.shape[1] > 2:
        aux = input_points[:, 2:]
        out_p = np.concatenate((poi[pos], aux, pos.reshape(pos.size, 1)), axis=1)
    else:
        out_p = poi[pos]

    # Remove duplicates
    if remove_duplicates:     
        idx = np.unique(out_p[:,-1], return_index=True)[1]
        out_p = out_p[idx, :-1]

    return out_p

kind="confluences"

fd = Flow("../data/in/small25_fd.tif")
poi = fd.get_stream_poi(500, kind,  "XY")

# Generate 10 random points
xmin, xmax, ymin, ymax = fd.get_extent()

np.random.seed(1234)
rnd_x = (xmax-xmin) * np.random.random(10)+xmin
rnd_y = (ymax-ymin) * np.random.random(10)+ymin

plt.scatter(rnd_x, rnd_y, c="r", s=10)
plt.scatter(poi[:,0], poi[:,1], c="b", marker="s", s=20)

input_points = np.array((rnd_x, rnd_y, np.arange(0, 1000, 100))).T

for n in range(10):
    plt.annotate(str(input_points[n,2]), (input_points[n,0], input_points[n,1]))


out_p = snap_points(fd, input_points, 500, kind, remove_duplicates=True)

plt.scatter(out_p[:,0], out_p[:,1], c="g", marker="x", s=15)
for n in range(out_p.shape[0]):
    plt.annotate(str(out_p[n,2]), (out_p[n,0], out_p[n,1]))
    
    
    