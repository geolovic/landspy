# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 17:36:31 2018

@author: Usuario
"""
import numpy as np
import ogr, osr

basedir = "/Users/vicen/Desktop/tunez/gisdata/BNetworks"
med_profiles = np.load(basedir + "/med_profiles.npy")
out_profiles = np.load(basedir + "/out_profiles.npy")

perfiles = np.append(med_profiles, out_profiles)

driver = ogr.GetDriverByName("ESRI Shapefile")
dataset = driver.CreateDataSource("/Users/vicen/Desktop/tunez/gisdata/main_channels.shp")
sp = osr.SpatialReference()
sp.ImportFromWkt(perfiles[0]._srs)
layer = dataset.CreateLayer("rios", sp, ogr.wkbLineString)
layer.CreateField(ogr.FieldDefn("id", ogr.OFTInteger))

for perfil in perfiles:
    feat = ogr.Feature(layer.GetLayerDefn())
    feat.SetField("id", int(perfil.rid))
    geom = ogr.Geometry(ogr.wkbLineString)
    xi = perfil.get_x()
    yi = perfil.get_y()
    for n in range(xi.size):
        geom.AddPoint(xi[n], yi[n])
        
    feat.SetGeometry(geom)
    layer.CreateFeature(feat)
    
layer = None
dataset = None