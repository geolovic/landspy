#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 11:10:26 2019

@author: vicen
"""

from topopy import DEM, Flow, Network
import numpy as np

dem = DEM("data/in/tunez.tif")
flw = Flow(dem)
net = Network(flw, 1000)

celdas = []

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

# Get cell orders
orders = self.get_stream_orders(kind="strahler", asgrid=False).ravel()

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
    elif gi[mid_cell] != 0:
        # Celda ya procesada, es una confluencia (canal tributario con 2 celdas)
        processing = False
    
    # Obtenemos datos de elevacion y distancias
    win = [mouth_cell, mid_cell, head_cell]
    xi = x_arr[ixcix[win]]
    yi = y_arr[ixcix[win]]
   
    # Calculamos pendiente de celda central por regresion
    poli, SCR = np.polyfit(xi, yi, deg = 1, full = True)[:2]
    
    # Calculamos R2
    if yi.size * yi.var() == 0:
        R2 = 1 # Puntos colineares
    else:
        R2 = float(1 - SCR/(yi.size * yi.var()))
    # Llenamos matrices auxiliares
    g = poli[0]
    gi[mid_cell] = g
    r2[mid_cell] = R2
    
    while processing:
        # Tomamos siguiente celda del canal (next_cell)   
        if ixcix[mouth_cell]==0:
            next_cell = mouth_cell
        else:
            next_cell = self._ixc[ixcix[mouth_cell]]

        if ixcix[mouth_cell] == 0:
            # Si estamos al final del canal, win decrece (elimina celda final - head_cell)
            win.pop()
            win.pop()
        elif ixcix[next_cell] == 0:
            # Cuando next_cell llega al final del canal
            # Se inserta la nueva celda y se elimina la celda de cabecera
            win.insert(0, next_cell)
            win.pop()
        elif head_cell == head and len(win)< winlen:
            # Si estamos a principio del canal, win crece
            win.insert(0, next_cell)
            next_cell = self._ixc[ixcix[next_cell]]
            win.insert(0, next_cell)
        else:
            # Se inserta la nueva celda y se elimina la celda de cabecera
            win.insert(0, next_cell)
            win.pop()
            
        head_cell = win[-1]
        mid_cell = self._ixc[ixcix[mid_cell]]
        mouth_cell = win[0]
        
        if len(win) == 3:
            processing = False
         
        # Obtenemos datos de elevacion y distancias para calcular pendientes
        xi = x_arr[ixcix[win]]
        yi = y_arr[ixcix[win]]
        poli, SCR = np.polyfit(xi, yi, deg = 1, full = True)[:2]
        g = poli[0]
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
            value_01 = gi[mid_cell]
            value_02 = g
#            if len(confs[mid_cell]) == 0:
#                confs[mid_cell].extend((value_01, value_02))
#            else:
#                confs[mid_cell].append(value_02)
            up_cell = win[win.index(mid_cell) + 1]
            ord_02 = orders[up_cell]
            cells = list(self._ix[np.where(self._ixc == mid_cell)])
            cells.remove(up_cell)
            if len(cells) > 1:
                for cell in cells:
                    if gi[cell] != 0:
                        print(cells, cell)
                        continue
            else:
                print(cells[0])
#            for cell in cells:
#                if cell != up_cell and gi[cell]!= 0:
#                    ord_01 = orders[cell]
#                else:
#                    ord_01 = 999
#                    
#            print(value_01, ord_01, value_02, ord_02)
#                    
            

#            
#            if mid_cell in confs.keys():
#                confs[mid_cell].append((g, win[win.index(mid_cell) - 1]))
#
#if kind == 'slp':
#    self._slp = gi[self._ix]
#    self._r2slp = r2[self._ix]
#    self._slp_npoints = npoints
#elif kind == 'ksn':
#    self._ksn = gi[self._ix]
#    self._r2ksn = r2[self._ix]
#    self._ksn_npoints = npoints
    