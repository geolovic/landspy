# -*- coding: utf-8 -*-

# p_distance.py
# Jose Vicente Perez Pena
# Dpto. Geodinamica-Universidad de Granada
# 18071 Granada, Spain
# vperez@ugr.es // geolovic@gmail.com
# 
# This module contains functions to calculate distances and cost-distances
#
# Version: 1.0
# January 15, 2018
#
# Last modified January 15, 2018

import numpy as np

ROWADD = [-1, 0, 1, 0, -1, -1, 1, 1]
COLADD = [0, 1, 0, -1, -1, 1, 1, -1]

    
def distance(in_arr):
    """
    This function calculates the distance for an array from seed locations.
    
    Parameters:
    ===========
    in_arr : numpy.array (INT data type). 
      Array with seed locations (values different to zero are considered seeds)
    """
    nrow, ncol = in_arr.shape
    # Distance array and queue
    dist_arr = np.zeros(in_arr.shape, dtype=np.float)
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
            if n <4:
                dist_add = 1
            else:
                dist_add = np.sqrt(2)
                
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
                dist_arr[row, col] = 0.0
                
    return dist_arr


#def cost(in_arr, cost_arr):
#    """
#    This function calculates the cost-distance for an array.
#    
#    Parameters:
#    ===========
#    in_arr : numpy.array (INT data type). 
#      Array with seed locations (values different to zero are considered seed)
#    cost_arr : numpy.array (FLOAT data type)
#      Array with costs (cells values will be used as cost for distance calculation)
#        
#    """
#    nrow, ncol = in_arr.shape
#    # Distance array, cost array and queue
#    dist_arr = np.zeros(in_arr.shape, dtype=np.float)
#    queue = []
#    
#    # Add seed values to que queue
#    for row in range(nrow):
#        for col in range(ncol):
#            if in_arr[row, col] > 0:
#                queue.append(row)
#                queue.append(col)
#                dist_arr[row, col] =  0.00001
#                
#    # Start loop that calculate distances
#    while len(queue) > 0:
#        
#        # Get cell from queue
#        row = queue.pop(0)
#        col = queue.pop(0)
#        mydist = dist_arr[row, col]
#        
#        # Get neighbors
#        for n in range(8):
#            nb_row = row + ROWADD[n]
#            nb_col = col + COLADD[n]
#            
#            # If cell is valid (border control)
#            if nb_row >= 0 and nb_col >= 0 and nb_row < nrow and nb_col < ncol:  
#                coste = (cost_arr[nb_row, nb_col] + cost_arr[row, col]) / 2
#                if n < 4:
#                    dist_add = coste
#                else:
#                    dist_add = coste * np.sqrt(2)
#                    
#                if dist_arr[nb_row, nb_col] == 0.0:
#                    # If neighbor is empty, enqueue and set distance
#                    dist_arr[nb_row, nb_col] = mydist + dist_add
#                    queue.append(nb_row)
#                    queue.append(nb_col)
#                
#                elif dist_arr[nb_row, nb_col] > mydist + dist_add:
#                    # Check if distance is larger that would be (whether cell processed before)
#                    dist_arr[nb_row, nb_col] = mydist + dist_add 
#                
#
#    # Put again zero values in seed locations
#    for row in range(nrow):
#        for col in range(ncol):
#            if in_arr[row, col] > 0:
#                dist_arr[row, col] = 0.0
#                
#    return dist_arr
    
def cost(in_arr, cost_arr):
    """
    This function calculates the cost-distance for an array.
    
    Parameters:
    ===========
    in_arr : numpy.array (INT data type). 
      Array with seed locations (values different to zero are considered seed)
    cost_arr : numpy.array (FLOAT data type)
      Array with costs (cells values will be used as cost for distance calculation)
        
    """
    # Get array dimensions as two integer variables
    nrow, ncol = in_arr.shape
    
    # Define distance array, aux-array and queue
    dist_arr = np.zeros(in_arr.shape, dtype=np.float)
    aux_arr = np.zeros(in_arr.shape, dtype=np.int8)
    queue = []
    
    # Add seed values to the queue and mark them as processed
    for row in range(nrow):
        for col in range(ncol):
            if in_arr[row, col] > 0:
                dist_arr[row, col] = 0.0001
                aux_arr[row, col] = 1
                queue.append(row)
                queue.append(col)
                               
    # Start loop that calculate distances
    while len(queue) > 0:
        # Get cell from queue and mark position as processed
        row = queue.pop(0)
        col = queue.pop(0)
        mydist = dist_arr[row, col]
        
        # Get cell neighbors
        for n in range(8):
            nb_row = row + ROWADD[n]
            nb_col = col + COLADD[n]
            
            # Check whether the neighbor is valid (border control)
            if not (nb_row >= 0 and nb_col >= 0 and nb_row < nrow and nb_col < ncol):
                continue

            # If the neighbor has not been processed, add it to queue
            if aux_arr[nb_row, nb_col] == 0:
                queue.append(nb_row)
                queue.append(nb_col)    
                aux_arr[nb_row, nb_col] = 1
            
            # Calculate the cost and the distance between cell and neighbor
            coste = (cost_arr[nb_row, nb_col] + cost_arr[row, col]) / 2
            if n < 4:
                dist_add = coste
            else:
                dist_add = coste * np.sqrt(2)

            # Update distance raster
            if dist_arr[nb_row, nb_col] == 0.0:
                dist_arr[nb_row, nb_col] = mydist + dist_add                 
            elif dist_arr[nb_row, nb_col] > mydist + dist_add:
                dist_arr[nb_row, nb_col] = mydist + dist_add
           
    # Put again zero values in seed locations
    for row in range(nrow):
        for col in range(ncol):
            if in_arr[row, col] > 0:
                dist_arr[row, col] = 0.0
                
    return dist_arr