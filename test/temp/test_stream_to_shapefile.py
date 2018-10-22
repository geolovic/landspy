#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 17:53:42 2018

@author: vicen
"""

from topopy import Network
import numpy as np
import matplotlib.pyplot as plt
import ogr, osr

net = Network("../data/in/jebja30_network.net")

self = net
#fig, ax = plt.subplots()

def export_to_shp(self, path=""):
    
    # Create shapefile
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataset = driver.CreateDataSource(path)
    sp = osr.SpatialReference()
    sp.ImportFromWkt(self._proj)
    layer = dataset.CreateLayer("rivers", sp, ogr.wkbLineString)
    layer.CreateField(ogr.FieldDefn("segid", ogr.OFTInteger))
    layer.CreateField(ogr.FieldDefn("order", ogr.OFTInteger))
    layer.CreateField(ogr.FieldDefn("flowto", ogr.OFTInteger))
    
    # Get channel segments and orders
    ch_seg = self.get_stream_segments(False).ravel()
    ch_ord = self.get_stream_order(asgrid=False).ravel()
    ch_seg = ch_seg[self._ix]
    ch_ord = ch_ord[self._ix]
    
    # Get ixcix auxiliar array
    ixcix = np.zeros(self._ncells, np.int)
    ixcix[self._ix] = np.arange(self._ix.size)
    
    seg_ids = np.unique(ch_seg)
    for idx in seg_ids:
        # skip zero channel (no channel)
        if idx == 0:
            continue   
        # Get givers for the segment id
        pos = np.where(ch_seg == idx)[0]
        ch_ix = self._ix[pos]
        
        # Add last point
        first = ch_ix[0]
        last = self._ixc[ixcix[ch_ix[-1]]]
        ch_ix = np.append(ch_ix, last)
        first = ch_ix[0]
        
        # Get segment order and receiver segment
        order = ch_ord[ixcix[first]]
        if ixcix[last] == 0:
            flowto = 0
        else:
            flowto = ch_seg[ixcix[last]]
        
        # Add feature
        feat = ogr.Feature(layer.GetLayerDefn())
        feat.SetField("segid", int(idx))
        feat.SetField("order", int(order))
        feat.SetField("flowto", int(flowto))
        row, col = self.ind_2_cell(ch_ix)
        xi, yi = self.cell_2_xy(row, col)
        
        geom = ogr.Geometry(ogr.wkbLineString)
        
        for n in range(xi.size):
            geom.AddPoint(xi[n], yi[n])
            
        feat.SetGeometry(geom)
        layer.CreateFeature(feat)
       
    layer = None
    dataset = None


def export_to_shp2(self, path=""):
    
    # Create shapefile
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataset = driver.CreateDataSource(path)
    sp = osr.SpatialReference()
    sp.ImportFromWkt(self._proj)
    layer = dataset.CreateLayer("rivers", sp, ogr.wkbLineString)
    layer.CreateField(ogr.FieldDefn("segid", ogr.OFTInteger))
    layer.CreateField(ogr.FieldDefn("order", ogr.OFTInteger))

    # Get heads and confluences
    heads = self.get_stream_poi(kind="heads", coords="IND")
    confs = self.get_stream_poi(kind="confluences", coords="IND")
    # Get channel orders
    ch_ord = self.get_stream_order(asgrid=False).ravel()
    ch_ord = ch_ord[self._ix]
    
    # Remove confluences of different orders
    confs_to_remove = []
    for conf in confs:
        givs = ch_ord[np.where(self._ixc == conf)]
        if np.unique(givs).size > 1:
            confs_to_remove.append(conf)    
    confs = confs.tolist()
    for conf in confs_to_remove:
        confs.remove(conf)     
    confs = np.array(confs)

    # Merge heads and confluences
    confs = np.append(heads, confs)

    # Get ixcix auxiliar array
    ixcix = np.zeros(self._ncells, np.int)
    ixcix[self._ix] = np.arange(self._ix.size)
    segid = 0
    
    # Iterate confluences
    for conf in confs:
        ch_idx = [conf]
        order = ch_ord[ixcix[conf]]
        next_cell = conf
        while ixcix[next_cell] != 0:
            next_cell = self._ixc[ixcix[next_cell]]
            ch_idx.append(next_cell)
            if ch_ord[ixcix[next_cell]] > order:
                break
        # Add feature
        feat = ogr.Feature(layer.GetLayerDefn())
        feat.SetField("segid", int(segid))
        feat.SetField("order", int(order))
        row, col = self.ind_2_cell(ch_idx)
        xi, yi = self.cell_2_xy(row, col)
        geom = ogr.Geometry(ogr.wkbLineString)
        for n in range(xi.size):
            geom.AddPoint(xi[n], yi[n])
            
        feat.SetGeometry(geom)
        layer.CreateFeature(feat)
        segid += 1
    