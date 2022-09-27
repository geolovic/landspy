#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 20 January, 2020
Testing suite for Network.get_chi_shapefile() function
@author: J. Vicente Perez
@email: geolovic@hotmail.com
@date: 20 January, 2020
@last_modified: 19 september, 2022
"""

import unittest, sys, os
# Add to the path code folder and data folder
sys.path.append("../src/")
from landspy import Network
import numpy as np
from osgeo import ogr, osr
import sys

infolder = "data/in"
outfolder = "data/out"

def canales_to_shapefile(path, canales):
    """
    Funcion de ayuda para guardar una serie de canales a un shapefile
    """
    # Creamos shapefile
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataset = driver.CreateDataSource(path)
    sp = osr.SpatialReference()
    proj = canales[0].getCRS()
    sp.ImportFromWkt(proj)
    
    # Creamos layer
    layer = dataset.CreateLayer("canales", sp, ogr.wkbLineString)
    
    # Add fields
    campos = ["oid", "name", "flowto", "thetaref", "chi0", "slp_np", "ksn_np"]
    tipos = [0, 4, 0, 2, 2, 0, 0]
    
    for n in range(len(tipos)):
        layer.CreateField(ogr.FieldDefn(campos[n], tipos[n]))
    
    # Add channels to shapefile
    for canal in canales:
        
        feat = ogr.Feature(layer.GetLayerDefn())
        feat.SetField("oid", int(canal.getOid()))
        feat.SetField("name", str(canal.getName()))
        feat.SetField("flowto", int(canal.getFlow()))
        feat.SetField("thetaref", float(canal._thetaref))
        feat.SetField("chi0", float(canal._chi0))
        feat.SetField("slp_np", int(canal._slp_np))
        feat.SetField("ksn_np", int(canal._ksn_np))
        geom = ogr.Geometry(ogr.wkbLineString)
        xy = canal.getXY()
        for row in xy:
            geom.AddPoint(row[0], row[1])
            
        feat.SetGeometry(geom)
        layer.CreateFeature(feat)

class NetworkGetChannelTest(unittest.TestCase):
    
    def test_get_channel_01(self):
        """
        Este test crea 5 puntos aleatorios dentro de la extension del objeto Network 
        (ids: 0, 2, 3, 5, 6) y dos fuera de la extension (ids: 1 y 4), y a partir de 
        los mismos crea objetos Channel hasta el outlet más cercano. Guarda 
        los canales generados como un shapefile con el nombre rdn_channels_{file}.shp
        """
        files = ["tunez", "jebja30"]
        for file in files:
            # Cargamos objeto network guardado previamente
            net_path = outfolder +  "/{0}_net.dat".format(file)
            net = Network(net_path)
            
            # Generamos 5 puntos aleatorios dentro de extension
            xmin, xmax, ymin, ymax = net.getExtent()
            puntos = []
            for n in range(5):
                x = np.random.randint(int(xmin), int(xmax))
                y = np.random.randint(int(ymin), int(ymax))
                puntos.append((x, y))
            
            # Creamos dos puntos fuera de extension
            pto = [xmin - (xmax - xmin) * 0.1, ymin - (ymax - ymin)*0.1]
            puntos.insert(1, pto)
            
            pto = [xmax + (xmax - xmin) * 0.1, ymax + (ymax - ymin)*0.1]
            puntos.insert(4, pto)
            
            canales = []
            
            for idx, pto in enumerate(puntos):
                canal = net.getChannel(pto, None, str(idx), idx)
                if canal:
                    canales.append(canal)

            out_shp = outfolder + "/rdn_channels_{}.shp".format(file)
            canales_to_shapefile(out_shp, canales)
            
            computed = os.path.exists(out_shp)
            self.assertEqual(computed, True)


    def test_get_channel_02(self):
        """
        Este test lee un archivo shapefile con la red de drenaje (creado anteriormente), 
        elige 5 lineas aleatorias y crea 5 objetos canal coincidiendo con esas líneas.
        Guarda los canales generados como un shapefile con el nombre shp_channels_{file}.shp
        """
        files = ["tunez", "jebja30"]
        for file in files:
            # Cargamos objeto network guardado previamente
            net_path = outfolder +  "/{0}_net.dat".format(file)
            net = Network(net_path)
            
            # Generamos un shapefile con todos los canales
            out_shp = outfolder +  "/{0}_str.shp".format(file)
            net.exportShp(out_shp)
            
            # Abrimos el shapefile creado
            driver = ogr.GetDriverByName("ESRI Shapefile")
            dataset = driver.Open(out_shp)
            layer = dataset.GetLayer(0)
            
            # Obtenemos el numero de entidades y seleccionamos 5 aleatorios
            nfeat = layer.GetFeatureCount()
            oids = np.random.choice(np.arange(nfeat), 5)
            
            canales = []
            # Obtenemos las entidades con esos Ids
            for oid in oids:
                feat = layer.GetFeature(oid)
                geom = feat.GetGeometryRef()
                if geom.GetGeometryCount() > 1:
                    continue
                
                head = geom.GetPoint(0)
                mouth = geom.GetPoint(geom.GetPointCount()- 1)
                canal = net.getChannel(head, mouth, str(oid), oid )
                if canal:
                    canales.append(canal)
            out_shp = outfolder + "/shp_channels_{}.shp".format(file)
            canales_to_shapefile(out_shp, canales)
            computed = os.path.exists(out_shp)
            self.assertEqual(computed, True)
          
            
if __name__ == "__main__":
    unittest.main()