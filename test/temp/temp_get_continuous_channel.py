# -*- coding: utf-8 -*-
"""
Editor de Spyder

Este es un archivo temporal
"""

from topopy import Network, Channel, Flow, BNetwork
import ogr, osr
import numpy as np

net = Network("../data/in/tunez_net.dat")

self = net
path = "../data/out/pruebacon_channels.shp"

# Create shapefile
driver = ogr.GetDriverByName("ESRI Shapefile")
dataset = driver.CreateDataSource(path)
sp = osr.SpatialReference()
sp.ImportFromWkt(self._proj)
layer = dataset.CreateLayer("rivers", sp, ogr.wkbLineString25D)
layer.CreateField(ogr.FieldDefn("order", ogr.OFTInteger))

# Get ixcix auxiliar array
ixcix = np.zeros(self.get_ncells(), np.int)
ixcix[self._ix] = np.arange(self._ix.size)

# Get heads, confluences, outlets and orders
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
                   
        if ixcix[next_cell] == 0: # next_cell is an outlet
            processing = False
        elif next_cell in strahler_confs: # next_cell is in strahler_confs
            processing = False
        else:
            cell = next_cell    
                    
    # Add feature
    feat = ogr.Feature(layer.GetLayerDefn())
    row, col = self.ind_2_cell(river_data)
    xi, yi = self.cell_2_xy(row, col)
    pos = ixcix[river_data]
    zi = self._zx[pos]            
    
    geom = ogr.Geometry(ogr.wkbLineString25D)
    for n in range(xi.size):
        geom.AddPoint(xi[n], yi[n], zi[n])
        
    feat.SetGeometry(geom)
    chanorder = ch_ord[cell]
    feat.SetField("order", int(chanorder))
    layer.CreateFeature(feat)
   
layer = None
dataset = None