#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  4 18:27:18 2018

@author: vicen
"""
import numpy as np

class Array():
    def __init__(self, arr, nodata):
        self._array = arr
        self._nodata = nodata

arr = np.arange(16).reshape(4, 4)


def values_2_nodata(arr, value):
    """
    Change specific values to NoData (if Grid nodata is defined). 
    
    Parameters:
    ===========
    value : Value or values that will be changed to NoData
    """
    if arr._nodata is None:
        return
    if type(value) == int or type(value)==float:
        ind = np.where(arr._array==value)
        arr._array[ind] = arr._nodata
    else:
        for val in value:
            ind = np.where(arr._array == val)
            arr._array[ind] = arr._nodata        
 
    
my_arr = Array(np.arange(16).reshape(4, 4), -99)

values_2_nodata(my_arr, 15)
values_2_nodata(my_arr, [0, 1, 4])
values_2_nodata(my_arr, np.array([5, 10]))
print(my_arr._array)