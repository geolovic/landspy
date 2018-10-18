# -*- coding: utf-8 -*-
"""
Editor de Spyder

Este es un archivo temporal
"""
import numpy as np
from topopy import DEM, Flow, Network
import matplotlib.pyplot as plt
import gdal, ogr, osr

dem = DEM("../data/in/jebja30.tif")
fd = Flow("../data/in/jebja30_fd.tif")
net = Network(dem, fd, 2000)
path = "../data/out/jebja_streams.shp"

# Prepare auxiliar arrays
seg_array = net.get_stream_segments(False).ravel()
seg_array = seg_array[net._ix]
ord_array = net.get_stream_order('strahler', False).ravel()
ord_array = ord_array[net._ix]
ixcix = np.zeros(net._ncells, np.int)
ixcix[net._ix] = np.arange(net._ix.size)

# Create output shapefile
driver = ogr.GetDriverByName("ESRI Shapefile")
dataset = driver.CreateDataSource(path)
sp = osr.SpatialReference()
sp.ImportFromWkt(net._proj)
layer = dataset.CreateLayer("rios", sp, ogr.wkbLineString)
layer.CreateField(ogr.FieldDefn("sid", ogr.OFTInteger))
layer.CreateField(ogr.FieldDefn("to", ogr.OFTInteger))
layer.CreateField(ogr.FieldDefn("order", ogr.OFTInteger))

# Iterate by segments id and take channel cells
for sid in np.unique(seg_array):
    
    # Get segment network cells
    seg_pos = np.where(seg_array == sid)
    seg_ind = net._ix[seg_pos]
    mouth = net._ixc[ixcix[seg_ind[-1]]]
    seg_ind = np.append(seg_ind, mouth)
    head = seg_ind[0]
    flow_to = seg_array[ixcix[mouth]]
    order = ord_array[ixcix[head]]
    row, col = net.ind_2_cell(seg_ind)
    x, y = net.cell_2_xy(row, col)
    
    # Create and add new feature to layer
    feat = ogr.Feature(layer.GetLayerDefn())
    geom = ogr.Geometry(ogr.wkbLineString) 
    for n in range(x.size):
        geom.AddPoint(x[n], y[n])
    feat.SetGeometry(geom)
    feat.SetField("sid", int(sid))
    feat.SetField("to", int(flow_to))
    feat.SetField("order", int(order))
    layer.CreateFeature(feat)
   
layer = None
dataset = None    