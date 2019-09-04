#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 11:10:26 2019

@author: vicen
"""

from topopy import DEM, Flow, Network
import numpy as np
import ogr, osr

infolder = "data/in"
outfolder = "data/out"
file = "tunez"

dem = DEM(infolder + "/{0}.tif".format(file))
flw = Flow(dem)
net = Network(flw)
self = net
path = "data/out/tunez_streams.shp"


#def _get_segmented_shp(self, path=""):
#"""
#Export Network channels to shapefile format. Channels will split in each confluence.
#
#path : str
#  Path to save the shapefile 
#"""
# Create shapefile
driver = ogr.GetDriverByName("ESRI Shapefile")
dataset = driver.CreateDataSource(path)
sp = osr.SpatialReference()
sp.ImportFromWkt(self._proj)
layer = dataset.CreateLayer("rivers", sp, ogr.wkbLineString)

# Get ixcix auxiliar array
ixcix = np.zeros(self._ncells, np.int)
ixcix[self._ix] = np.arange(self._ix.size)

# Get heads and sort them by elevation
heads = self.get_stream_poi("heads", "IND")
elev = self._zx[ixcix[heads]]
spos = np.argsort(-elev)
heads = heads[spos]
aux_arr = np.zeros(self._ncells, dtype = np.bool)

for head in heads:
    cell = head
    river_data = [cell]
    aux_arr[cell] == True
    processing = True
    while processing:
        next_cell = self._ixc[ixcix[cell]]
        river_data.append(next_cell)
        if ixcix[next_cell] == 0 or aux_arr[next_cell]==True:
            processing = False
        else:
            river_data.append(next_cell)
            cell = next_cell
            aux_arr[cell] = True
        
                    
    # Add feature
    feat = ogr.Feature(layer.GetLayerDefn())
    row, col = self.ind_2_cell(river_data)
    xi, yi = self.cell_2_xy(row, col)
    
    geom = ogr.Geometry(ogr.wkbLineString)
    
    for n in range(xi.size):
        geom.AddPoint(xi[n], yi[n])
        
    feat.SetGeometry(geom)
    layer.CreateFeature(feat)
   
layer = None
dataset = None