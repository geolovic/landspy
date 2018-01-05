#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  5 18:09:39 2018

@author: vicen
"""

import numpy as np

# Leemos archivo de texto
data_arr = np.genfromtxt("MY_GRID.txt", dtype=None, delimiter=";", names=True)

# Generamos 100 puntos aleatorios
rnd_ids = np.random.randint(0,len(data_arr), 100)

arr = []
for name in data_arr.dtype.names:
    arr = data_arr[name][rnd_ids]
    np.save("MY_GRID_100rnd_" + name, arr)

