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
from scipy.sparse import csc_matrix

fd = Flow()
fd.load_gtiff("../data/tunez_fd.tif")
threshold = 1000

def get_stream_segments(self, threshold, asgrid=True):
    """
    TODO - Write documentation
    """
    # Create the flow accumulation file
    fac = fd.get_flow_accumulation(nodata = False, asgrid=False)
    w = fac > threshold    
    # Build a sparse array with giver-receivers cells
    w = w.ravel()
    I   = w[fd._ix]
    ix  = fd._ix[I]
    ixc = fd._ixc[I]
    aux_vals = np.ones(ix.shape, dtype=np.int8)
    
    # We don't use get_stream_poi because it implies creating the flow accumulation 3x
    sp_arr = csc_matrix((aux_vals, (ix, ixc)), shape=(fd._ncells, fd._ncells))
    
    del I, aux_vals, fac # Clean up
    
    # Get heads
    sum_arr = np.asarray(np.sum(sp_arr, 0)).ravel()
    head_pos = (sum_arr == 0) & w
    head_ind = np.where(head_pos)
    
    # Get confluences
    sum_arr = np.asarray(np.sum(sp_arr, 0)).ravel()
    conf_pos = sum_arr > 1
    conf_ind = np.where(conf_pos)
    
    del sum_arr, sp_arr, w, conf_pos, head_pos # Clean up
    
    # Merge all heads and confluences
    all_ind = np.append(head_ind, conf_ind)
    
    del conf_ind, head_ind # Clean up
    
    # We created a zeros arrays and put in confuences and heads their id
    # Those id will be consecutive numbers starting in one
    seg_arr = np.zeros(fd._ncells, dtype=np.int32)
    for n, inds in enumerate(all_ind):
        seg_arr[inds] = n+1
    
    # Move throught giver list. If receiver is 0, give the same id that giver.
    # If a receiver is not 0, that means that we are in a confluence. 
    for n in range(len(ix)):
        if seg_arr[ixc[n]] == 0:
            seg_arr[ixc[n]] = seg_arr[ix[n]]
    
    # Reshape and output
    seg_arr = seg_arr.reshape(fd._dims)
    if asgrid:
        return fd._create_output_grid(seg_arr, 0)
    else:
        return seg_arr

links = get_stream_segments(threshold, True)
links.save("tunez_segments.tif")