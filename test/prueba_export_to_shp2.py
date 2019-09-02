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
file = "jebja30"

dem = DEM(infolder + "/{0}.tif".format(file))
flw = Flow(dem)
net = Network(flw)
self = net
path = "data/out/jebja30_streams.shp"


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
layer.CreateField(ogr.FieldDefn("order", ogr.OFTInteger))

# Get ixcix auxiliar array
ixcix = np.zeros(self._ncells, np.int)
ixcix[self._ix] = np.arange(self._ix.size)

# Get heads, confluences and orders
heads = self.get_stream_poi("heads", "IND")
confs = self.get_stream_poi("confluences", "IND")
ch_ord = self.get_stream_orders(asgrid=False).ravel()

# Get confluences where strahler index increases
strahler_confs = []
for conf in confs:
    conf_order = ch_ord[conf]
    givers = self._ix[np.where(self._ixc==conf)]
    giv_orders = ch_ord[givers]
    if giv_orders.max() < conf_order:
        strahler_confs.append(conf)
        
# Append strahler confluences to heads
heads = np.append(heads, np.array(strahler_confs))    

# Iterate heads
for head in heads:
    cell = head
    river_data = [cell]
    processing = True
    while processing:
        next_cell = self._ixc[ixcix[cell]]
        river_data.append(next_cell)
        if next_cell in confs:
            if ch_ord[next_cell] > ch_ord[cell]:
                processing = False
            else:
                river_data.append(next_cell)
                cell = next_cell
        elif ixcix[next_cell] == 0:
            processing = False
        else:
            river_data.append(next_cell)
            cell = next_cell    
                    
    # Add feature
    feat = ogr.Feature(layer.GetLayerDefn())
    row, col = self.ind_2_cell(river_data)
    xi, yi = self.cell_2_xy(row, col)
    geom = ogr.Geometry(ogr.wkbLineString)
    for n in range(xi.size):
        geom.AddPoint(xi[n], yi[n])
        
    feat.SetGeometry(geom)
    chanorder = ch_ord[cell]
    feat.SetField("order", int(chanorder))
    layer.CreateFeature(feat)
   
layer = None
dataset = None