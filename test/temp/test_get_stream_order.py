#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 09:42:34 2018

@author: vicen
"""
import sys
sys.path.append("../../")
from topopy import Flow, Grid
import numpy as np

fd = Flow()
fd.load_gtiff("../data/tunez_fd.tif")

def get_stream_order(fd, threshold, kind="strahler", asgrid=True):
    """
    This function extract stream orders by using a determined area threshold

    Parameters:
    ===========
    threshold : *int* 
      Area threshold to initiate a channels (in cells)
    kind : *str* {'strahler', 'shreeve'}
    asgrid : *bool*
      Indicates if the network is returned as topopy.Grid (True) or as a numpy.array
    """
    if kind not in ['strahler', 'shreeve']:
        return
    
    fac = fd.get_flow_accumulation(nodata=False, asgrid=False)
    w = fac > threshold
    w = w.ravel()
    I   = w[fd._ix]
    ix  = fd._ix[I]
    ixc = fd._ixc[I]

    str_ord = np.copy(w).astype(np.int8)
    visited = np.zeros(fd._ncells, dtype=np.int8)

    if kind == 'strahler':
        for n in range(len(ix)):
            if (str_ord[ixc[n]] == str_ord[ix[n]]) & visited[ixc[n]]:
                str_ord[ixc[n]] = str_ord[ixc[n]] + 1                
            else:
                str_ord[ixc[n]] = max(str_ord[ix[n]], str_ord[ixc[n]])
                visited[ixc[n]] = True
    elif kind == 'shreeve':
        for n in range(len(ix)):
            if visited[ixc[n]]:
                str_ord[ixc[n]] = str_ord[ixc[n]] + str_ord[ix[n]]
            else:
                str_ord[ixc[n]] = max(str_ord[ix[n]], str_ord[ixc[n]])
                visited[ixc[n]] = True
    str_ord = str_ord.reshape(fd._dims)
    
    if asgrid:
        return fd._create_output_grid(str_ord, nodata_value=0)
    else:
        return str_ord

str_strahler = get_stream_order(fd, 1000, kind="strahler")
str_strahler.save("../data/tunez_strahler.tif")

str_shreeve = get_stream_order(fd, 1000, kind="shreeve")
str_shreeve.save("../data/tunez_shreeve.tif")