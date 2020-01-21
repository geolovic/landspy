#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 18:13:34 2020

@author: vicen
"""

from topopy import Network
import numpy as np
import matplotlib.pyplot as plt
import ogr, osr

def cells_to_points(cells, net, outfile):
    row, col = net.ind_2_cell(cells)
    x, y = net.cell_2_xy(row, col)
    np.savetxt(outfile, np.array((x, y)).T, delimiter=";", header="X;Y", comments="" )

aux_cells = []


files = ["small25", "morocco", "tunez", "jebja30"]
file = files[2]

net_path = "../data/in/{0}_network.net".format(file)
net = Network(net_path)
self = net

out_shp = "../data/out/{0}_chishp.shp".format(file)
distance = 250

#def get_chi_shapefile(self, out_shp, distance):
"""
This method export network data to a shapelife. It calculates segments of a given
distance and calculate chi, ksn, slope, etc. for the segment. The shapefile will 
have the following fields:
    id_profile : Profile identifier. Profiles are calculated from heads until outlets or
    L : Lenght from the segment middle point to the head
    area_e6 : Drainage area in the segment mouth (divided by E6, to avoide large numbers)
    z : Elevation of the middle point of the segment
    chi : Mean chi of the segment
    ksn : Ksn of the segment (calculated by linear regression)
    slope : slope of the segment (calculated by linear regression)
    rksn : R2 of the ksn linear regression
    rslope : R2 of the slope linear regression
    
Parameters:
===========
out_shp : srt
  Output shapefile
distance : float
  Segment distance. 
    
"""
# Create shapefile
driver = ogr.GetDriverByName("ESRI Shapefile")
dataset = driver.CreateDataSource(out_shp)
sp = osr.SpatialReference()
sp.ImportFromWkt(self._proj)
layer = dataset.CreateLayer("rivers", sp, ogr.wkbLineString)

# Add fields
campos = ["id_profile", "L", "area_e6", "z", "chi", "ksn", "rksn", "slope", "rslope"]
tipos = [0, 2, 2, 2, 2, 2, 2, 2, 2]
for n in range(len(campos)):
    layer.CreateField(ogr.FieldDefn(campos[n], tipos[n]))

# Get ixcix auxiliar array
ixcix = np.zeros(self._ncells, np.int)
ixcix[self._ix] = np.arange(self._ix.size)

# Get heads and sort them by elevation and iterate them
heads = self.get_stream_poi("heads", "IND")
zpos = np.argsort(self._zx[ixcix[heads]])
heads = heads[zpos][::-1]

aux_arr = np.zeros(self._ncells, np.bool)
id_profile = 0
for head in heads:
    id_profile += 1
    processing = True
    cell = head
    segment_cells = [cell]
    segment_distance = 0
    profile_length = self._dx[ixcix[head]]
    
    while processing:
        add_segment = False
        # Take the next cell downstream and add it to the array
        next_cell = self._ixc[ixcix[cell]]
        segment_cells.append(next_cell)
        segment_distance += self._dd[ixcix[cell]]

        if segment_distance >= distance:
            # If segment distance is reached, add the current segment
            add_segment = True            
        
        if aux_arr[next_cell] == True or ixcix[next_cell] == 0:
            # If the cell is processed or an outlet
            # End processing and add the segment (even if distance not reached)
            processing = False
            add_segment = True
                    
        aux_arr[next_cell] = True
        cell = next_cell

        # Add segment to the shapefile (if has more than 3 cells)
        if add_segment and len(segment_cells) >= 3:
            # Get point coordinates and values
            row, col = self.ind_2_cell(segment_cells)
            xi, yi = self.cell_2_xy(row, col)
            pos = ixcix[segment_cells][::-1]
            mouth_cell = ixcix[segment_cells[-1]]
            mid_cell = ixcix[segment_cells[int(len(segment_cells)/2)]]
            dx = self._dx[pos]
            zx = self._zx[pos]
            chi = self._chi[pos]
            area = self._ax[mouth_cell]
            lenght = profile_length - self._dx[mid_cell]
            mid_z = self._zx[mid_cell]
            mid_chi = self._chi[mid_cell]      
            
            # Mean slope for the segment (by min. sqares)
            poli, SCR = np.polyfit(dx, zx, deg = 1, full = True)[:2]
            slope = poli[0]
            if slope == 0: 
                slope = 0.000001
            if zx.size * zx.var() == 0:
                rslope = 1
            else:
                rslope = float(1 - SCR/(zx.size * zx.var()))      
            
            # Mean ksn for the segment (by min. sqares)
            poli, SCR = np.polyfit(chi, zx, deg = 1, full = True)[:2]
            ksn = poli[0]
            if ksn == 0: 
                ksn = 0.000001
            if zx.size * zx.var() == 0:
                rksn = 1
            else:
                rksn = float(1 - SCR/(zx.size * zx.var()))
            
            # Create feature and add attributes
            feat = ogr.Feature(layer.GetLayerDefn())
            feat.SetField("id_profile", int(id_profile))
            feat.SetField("L", float(lenght))
            feat.SetField("area_e6", float(area/1000000))
            feat.SetField("z", float(mid_z))
            feat.SetField("chi", float(mid_chi))
            feat.SetField("ksn", float(ksn))
            feat.SetField("rksn", float(rksn))
            feat.SetField("slope", float(slope))
            feat.SetField("rslope", float(rslope))
            
            # Create geometry
            geom = ogr.Geometry(ogr.wkbLineString)
            for n in range(xi.size):
                geom.AddPoint(xi[n], yi[n])
            feat.SetGeometry(geom)
            # Add segment feature to the shapefile
            layer.CreateFeature(feat)
            # Reset variables
            segment_cells = [next_cell]
            segment_distance = 0

layer = None
dataset = None

