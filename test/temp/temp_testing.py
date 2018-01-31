#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 17:14:39 2018
Temp testing script for Network
@author: vicen
"""
import sys
sys.path.append("../../")
from topopy import DEM, Network
from scipy.sparse import csc_matrix
import numpy as np

dem = DEM("../data/small25.tif")
fd = Network()
fd.load_gtiff("sample_fd.tif")
fac = fd.flow_accumulation()
fac.save("small25_fac.tif")
nodata = fac.get_nodata_pos()
fac = fac.read_array()
fac[nodata] = 0

# Get pixels larger than threshold
w = fac > 1000
siz = fd._dims
nrc = fd._ncells

w = w.ravel()
I   = w[fd._ix]
ix  = fd._ix[I]
ixc = fd._ixc[I]

aux_vals = np.ones(ix.shape, dtype=np.int8)

sparr = csc_matrix((aux_vals, (ix, ixc)), shape=(nrc, nrc))

heads =  (np.sum(sparr, 0) == 0) & w
heads = heads.reshape(fd._dims)

row, col = np.where(heads)

xi, yi = dem.cell_2_xy(row, col)

outf = open("heads.txt", "w")
outf.write("X;Y\n")
for n in range(len(xi)):
    linea = "{0};{1}\n".format(xi[n], yi[n])
    outf.write(linea)
    
outf.close()
    


