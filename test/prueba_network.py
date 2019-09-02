#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 11:10:26 2019

@author: vicen
"""

from topopy import DEM, Flow, Network
import numpy as np

infolder = "data/in"
file = "tunez"
dem_path = infolder + "/{0}.tif".format(file)
dem = DEM(dem_path)
flw = Flow(dem)
net = Network(flw)
self = net
npoints = 4
kind = 'slp' 


#def calculate_gradients(self, npoints, kind='slp'):
"""
This function calculates gradients (slope or ksn) for all channel cells. 
Gradients of each cell are calculated by linear regression using a number
of points (npoints) up and downstream.

Parameters:
===========
npoints : *int*
  Window to analyze slopes. Slopes are calculated by linear regression using a window
  of npoints * 2 + 1 pixel (using the central pixel)
  
kind : *str* {'slp', 'ksn'}
"""
if kind not in ['slp', 'ksn']:
    kind = 'slp'
winlen = npoints * 2 + 1

# Get arrays depending on type
if kind == 'slp':
    x_arr = self._dx
    y_arr = self._zx
elif kind == 'ksn':
    x_arr = self._chi
    y_arr = self._zx
    
# Get ixcix auxiliar array
ixcix = np.zeros(self._ncells, np.int)
ixcix[self._ix] = np.arange(self._ix.size)

# Get heads array and confluences dictionary
heads = self.get_stream_poi("heads", "IND")
confs = {conf:[] for conf in self.get_stream_poi("confluences", "IND")}

# Sort heads by elevation
elev = self._zx[ixcix[heads]]
spos = np.argsort(-elev)
heads = heads[spos]

# Prepare auxiliary arrays
gi = np.zeros(self._ncells)
r2 = np.zeros(self._ncells)

# Taking sequentally all the heads and compute downstream flow
for head in heads:
    processing = True
    head_cell = head
    mid_cell = self._ixc[ixcix[head_cell]]
    mouth_cell = self._ixc[ixcix[mid_cell]]

    if ixcix[mid_cell] == 0:
        # Channel type 2 (mid_cell is an outlet)
        continue
    elif ixcix[mouth_cell]== 0:
        # Channel type 3 (mouth_cell is an outlet)
        processing = False
    # Check si cabecera es una confluencia de solo dos celdas
    elif gi[mid_cell] != 0:
        processing = False
    
    # Obtenemos datos de elevacion y distancias
    win = [mouth_cell, mid_cell, head_cell]
    xi = x_arr[ixcix[win]]
    yi = y_arr[ixcix[win]]
   
    # Calculamos pendiente de celda central por regresion
    poli, SCR = np.polyfit(xi, yi, deg = 1, full = True)[:2]
    g = poli[0]
    if g == 0: 
        g = 0.000001

    # Calculamos gradient y R2
    if yi.size * yi.var() == 0:
        R2 = 1 # Puntos colineares
    else:
        R2 = float(1 - SCR/(yi.size * yi.var()))

    gi[mid_cell] = g
    r2[mid_cell] = R2
            
    while processing:
        # Verificamos si estamos al final (en un outlet)  
        if not ixcix[mouth_cell]==0:
            # Si mouth_cell no es un outlet, cogemos siguiente celda
            next_cell = self._ixc[ixcix[mouth_cell]]
            # Si la siguiente celda es el final, 
            # añadimos una celda y eliminamos la cabecera
            if ixcix[next_cell]==0:
                win.insert(0, next_cell)
                win.pop()
            # Si longitud de ventana < winlen, se añaden dos celdas 
            elif len(win) < winlen:
                win.insert(0, next_cell)
                aux_cell = self._ixc[ixcix[next_cell]]
                win.insert(0, aux_cell)
            else:
                win.insert(0, next_cell)
                win.pop() 
        else:
            # Si mouth_cell es un outlet, no se coge siguiente celda
            next_cell = mouth_cell
            win.pop()
            win.pop()
        
        head_cell = win[-1]
        mid_cell = self._ixc[ixcix[mid_cell]]
        mouth_cell = win[0]
        
        if len(win) <= 3:
            processing = False
            
        # Obtenemos datos de elevacion y distancias para calcular pendientes
        xi = x_arr[ixcix[win]]
        yi = y_arr[ixcix[win]]
        poli, SCR = np.polyfit(xi, yi, deg = 1, full = True)[:2]
        g = poli[0]
        if g == 0: 
            g = 0.000001
        
        if yi.size * yi.var() == 0:
            R2 = 1
        else:
            R2 = float(1 - SCR/(yi.size * yi.var()))
                    
        # Comprobamos si celda central ha sido procesada
        if gi[mid_cell] == 0:
            gi[mid_cell] = g
            r2[mid_cell] = R2
        else:
            # Si ha sido procesada, es una confluencia 
            processing = False
            if len(confs[mid_cell]) == 0:
                confs[mid_cell].append(gi[mid_cell])
                confs[mid_cell].append(g)
            else:
                confs[mid_cell].append(g)
                
# Calculamos valores medios en confluencias                
for cell in confs.keys():
    if len(confs[cell]) > 0:
        gi[cell] = np.mean(np.array(confs[cell]))

# Llenamos array del objeto Network
if kind == 'slp':
    self._slp = gi[self._ix]
    self._r2slp = r2[self._ix]
    self._slp_npoints = npoints
elif kind == 'ksn':
    self._ksn = gi[self._ix]
    self._r2ksn = r2[self._ix]
    self._ksn_npoints = npoints