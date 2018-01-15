# -*- coding: utf-8 -*-

# c_distance.py
# Jose Vicente Perez Pena
# Dpto. Geodinamica-Universidad de Granada
# 18071 Granada, Spain
# vperez@ugr.es // geolovic@gmail.com
# 
# This module contains functions to calculate distances and cost-distances
# It uses cython to improve efficiency
#
# Version: 1.0
# January 15, 2018
#
# Last modified January 15, 2018

import numpy as np
cimport numpy as np

cdef int ROWADD[8] 
cdef int COLADD[8]
cdef double SQRT_2 = np.sqrt(2)

ROWADD[:] = [-1, 0, 1, 0, -1, -1, 1, 1]
COLADD[:] = [0, 1, 0, -1, -1, 1, 1, -1]

def distance(in_arr):
    valid_arr = in_arr.astype(np.int)
    res = _distance(valid_arr)
    return res
    
def _distance(np.ndarray[np.int_t, ndim=2] in_arr):
    """
    This function receives a numpy int array and calculate distances for those
    values different from zero
    """
    cdef int nrow, ncol
    cdef int row, col, nb_row, nb_col
    cdef size_t n
    cdef np.ndarray[np.float64_t, ndim=2] dist_arr
    cdef double mydist, dist_add
    
    nrow = in_arr.shape[0]
    ncol = in_arr.shape[1]
    
    # Distance array and queue
    dist_arr = np.zeros((nrow, ncol), dtype=np.float64)
    queue = []
    
    # Add seed values to que queue
    for row in range(nrow):
        for col in range(ncol):
            if in_arr[row, col] > 0:
                queue.append(row)
                queue.append(col)
                dist_arr[row, col] =  0.00001
                
    # Start loop that calculate distances
    while len(queue) > 0:
        
        # Get cell from queue
        row = queue.pop(0)
        col = queue.pop(0)
        mydist = dist_arr[row, col]
        
        # Get neighbors
        for n in range(8):
            nb_row = row + ROWADD[n]
            nb_col = col + COLADD[n]
            if n < 4:
                dist_add = 1.
            else:
                dist_add = SQRT_2
                
            # If cell is valid (border control)
            if nb_row >= 0 and nb_col >= 0 and nb_row < nrow and nb_col < ncol:                 
                if dist_arr[nb_row, nb_col] == 0.0:
                    # If neighbor is empty, enqueue and set distance
                    dist_arr[nb_row, nb_col] = mydist + dist_add
                    queue.append(nb_row)
                    queue.append(nb_col)
                
                elif dist_arr[nb_row, nb_col] > mydist + dist_add:
                    # Check if distance is larger that would be (whether cell processed before)
                    dist_arr[nb_row, nb_col] = mydist + dist_add 
            

    # Put again zero values in seed locations
    for row in range(nrow):
        for col in range(ncol):
            if in_arr[row, col] > 0:
                in_arr[row, col] = 0
                
    return np.asarray(dist_arr)

def cost(in_arr, cost_arr):
    valid_in_arr = in_arr.astype(np.int)
    valid_cost_arr = cost_arr.astype(np.float)
    res = _cost(valid_in_arr, valid_cost_arr)
    return res

def _cost(np.ndarray[np.int_t, ndim=2] in_arr, np.ndarray[np.float_t, ndim=2] cost_arr):
    """
    This function receives a numpy int array and calculate distances for those
    values different from zero
    """
    cdef int nrow, ncol
    cdef int row, col, nb_row, nb_col
    cdef size_t n
    cdef np.ndarray[np.float_t, ndim=2] dist_arr
    cdef double mydist, dist_add, coste
    
    nrow = in_arr.shape[0]
    ncol = in_arr.shape[1]
    
    # Distance array and queue
    dist_arr = np.zeros((nrow, ncol), dtype=np.float)
    queue = []
    
    # Add seed values to que queue
    for row in range(nrow):
        for col in range(ncol):
            if in_arr[row, col] > 0:
                queue.append(row)
                queue.append(col)
                dist_arr[row, col] =  0.00001
                
    # Start loop that calculate distances
    while len(queue) > 0:
        
        # Get cell from queue
        row = queue.pop(0)
        col = queue.pop(0)
        mydist = dist_arr[row, col]
        
        # Get neighbors
        for n in range(8):
            nb_row = row + ROWADD[n]
            nb_col = col + COLADD[n]
            
            # If cell is valid (border control)
            if nb_row >= 0 and nb_col >= 0 and nb_row < nrow and nb_col < ncol:  
                coste = (cost_arr[nb_row, nb_col] + cost_arr[row, col]) / 2
                if n < 4:
                    dist_add = coste
                else:
                    dist_add = SQRT_2 * coste
                                  
                if dist_arr[nb_row, nb_col] == 0.0:
                    # If neighbor is empty, enqueue and set distance
                    dist_arr[nb_row, nb_col] = mydist + dist_add
                    queue.append(nb_row)
                    queue.append(nb_col)
                
                elif dist_arr[nb_row, nb_col] > mydist + dist_add:
                    # Check if distance is larger that would be (whether cell processed before)
                    dist_arr[nb_row, nb_col] = mydist + dist_add 

    # Put again zero values in seed locations
    for row in range(nrow):
        for col in range(ncol):
            if in_arr[row, col] > 0:
                in_arr[row, col] = 0
                
    return np.asarray(dist_arr)