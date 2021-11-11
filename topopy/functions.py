# -*- coding: utf-8 -*-

# network.py
# Jose Vicente Perez Pena
# Dpto. Geodinamica-Universidad de Granada
# 18071 Granada, Spain
# vperez@ugr.es // geolovic@gmail.com
#
# MIT License (see LICENSE file)
# Version: 1.0
# December 26, 2017
#
# Last modified February 24, 2020

import numpy as np
from osgeo import ogr

def extract_points(path, idfield=""):
    """
    Extract coordinates from a point shapefile
    
    Parameters:
    ================
    path : str
      Path to the shapefile
    idfield : str
      Shapefile field with point ids
      
    Return:
    ================
    coords : np.dnarray
      Numpy array with 2 or 3 columns with points X and Y (thrid column will contain point ids
      if idfield was specified)
    """
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataset = driver.Open(path)
    layer = dataset.GetLayer()
    geom_type = layer.GetGeomType()
    lydef = layer.GetLayerDefn()
    id_fld = lydef.GetFieldIndex(idfield)
    geom_type = layer.GetGeomType()
    if geom_type != 1:
        return(np.array([0, 0, 0]).reshape(1, 3))
    points = []
    for feat in layer:
        geom = feat.GetGeometryRef()
        if geom is None:
            continue
        
        if id_fld >= 0:
            idx = feat.GetField(idfield)
            points.append((geom.GetX(), geom.GetY(), int(idx)))
        else:
            points.append((geom.GetX(), geom.GetY()))

    return np.array(points)


def rivers_to_channels(path, net, idfield=""):
    
    # Open que river shapefile
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataset = driver.Open(path)
    layer = dataset.GetLayer()
    geom_type = layer.GetGeomType()
    lydef = layer.GetLayerDefn()
    id_fld = lydef.GetFieldIndex(idfield)
    geom_type = layer.GetGeomType()
    