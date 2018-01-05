#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 17:15:30 2017

@author: vicen
"""
import gdal
import numpy as np

MY_GRID = "small25.tif"

def grid_to_xyz(grid, out_file):
    
    # Read input raster and get numpy array
    raster = gdal.Open(grid)
    banda = raster.GetRasterBand(1)
    nodata = banda.GetNoDataValue()
    geot = raster.GetGeoTransform()
    datarr = banda.ReadAsArray()
    xsize = banda.XSize
    ysize = banda.YSize
    
    # Open output text file and writes head
    outf = open(out_file, "w")
    outf.write("id;row;col;X;Y;Z\n")
    
    # Read raster band cell by cell and output xyz (in raster coordinates)
    n = 0
    for value in datarr.ravel():
        z = value
        row, col = np.unravel_index(n, (ysize, xsize))
        x = geot[0] + geot[1] * col + geot[1] / 2
        y = geot[3] - geot[1] * row - geot[1] / 2
        
        linea = "{0};{1};{2};{3:.2f};{4:.2f};{5:.2f}\n".format(n, row, col, x, y, z)
        outf.write(linea)
            
        n += 1
    
    banda = None
    datarr = None
    raster = None
    outf.close()
    

if __name__ == "__main__":
    grid_to_xyz(MY_GRID, "MY_GRID.txt")
    