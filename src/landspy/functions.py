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


def shp_to_channels(path, net, id_field="", name_field=""):
    """
    Extracts channels from river shapefile. Shapefile geometry must be polyline. Only the first point (head)
    and the last point (mouth) of each line are used, the channel is extracted from the Network object. Only no-multipart
    polylines will be used. 
    
    Parameters:
    ================
    path : str
      Path to the polyline shapefile
    net : landspy.Network
      Network instance
    id_field : str 
      Field with the OId for the channel. If not provided, the oid will be an integer number (order of the line)
    name_field : str 
      Field with the label for the channel. If not provided, the name will be an integer number (order of the line)
      
    Return:
    ================
    channel : list of landspy.Channel objects
      List of landspy.Channel objects
    """
    
    # Get driver and open shapefile
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataset = driver.Open(path)
    layer = dataset.GetLayer()
    lydef = layer.GetLayerDefn()
    namefield = lydef.GetFieldIndex(name_field)
    idfield = lydef.GetFieldIndex(id_field)
    
    
    channels = []
    
    for n, feat in enumerate(layer):
        geom = feat.GetGeometryRef()
        
        # If the feature is multigeometry (more than one geometry), discard that feature
        if geom.GetGeometryCount() > 1:
            continue
        
        # Get first and last points
        head = geom.GetPoint(0)
        mouth = geom.GetPoint(geom.GetPointCount() - 1)
        
        # If namefield is in the shapefile, get it
        if namefield >= 0:
            name = feat.GetField(namefield)
        else:
            name = str(n)
            
        if idfield >= 0:
            oid = feat.GetField(idfield)
        else:
            oid = n
        
        channel = net.getChannel(head, mouth, name, oid)
        channels.append(channel)

    return channels