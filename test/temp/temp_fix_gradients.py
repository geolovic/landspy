# -*- coding: utf-8 -*-
"""
Created on Mon Nov 15 12:05:10 2021

@author: Usuario
"""

# def calculate_gradients(self, npoints, kind='slp'):
#     """
#     This function calculates gradients (slope or ksn) for all channel cells. 
#     Gradients of each cell are calculated by linear regression using a number
#     of points (npoints) up and downstream.
    
#     Parameters:
#     ===========
#     npoints : *int*
#       Window to analyze slopes. Slopes are calculated by linear regression using a window
#       of npoints * 2 + 1 pixel (using the central pixel)
      
#     kind : *str* {'slp', 'ksn'}
#     """

from topopy import Network
import numpy as np

net = Network("../data/in/jebja30_net.dat")
self = net
npoints = 5
kind = "slp"



winlen = npoints * 2 + 1

# Get arrays depending on type
y_arr = self._zx
if kind == 'ksn':
    x_arr = self._chi
else:
    x_arr = self._dx
    
# Get ixcix auxiliar array
ixcix = np.zeros(self.get_ncells(), int)
ixcix[self._ix] = np.arange(self._ix.size)

# Get heads array and confluences dictionary
heads = self.get_stream_poi("heads", "IND")
#confs = {conf:[] for conf in self.get_stream_poi("confluences", "IND")}

# Sort heads by elevation
elev = self._zx[ixcix[heads]]
spos = np.argsort(-elev)
heads = heads[spos]

heads = [160285]

# Prepare auxiliary arrays
gi = np.zeros(self.get_ncells())
r2 = np.zeros(self.get_ncells())

# Taking sequentally all the heads and compute downstream flow
for head in heads:
    processing = True
    head_cell = head
    mid_cell = self._ixc[ixcix[head_cell]]
    mouth_cell = self._ixc[ixcix[mid_cell]]
    win = [mouth_cell, mid_cell, head_cell]
    
    if ixcix[mid_cell] == 0 or gi[mid_cell] != 0:
        # Channel type 2 (mid_cell is an outlet) or mid_cell is already calculated
        continue
    
    elif ixcix[mouth_cell]== 0:
        # Channel type 3 (mouth_cell is an outlet)
        processing = False
        win = [mid_cell, mid_cell, head_cell]
        
    # Obtenemos datos de elevacion y distancias
    xi = x_arr[ixcix[win]]
    yi = y_arr[ixcix[win]]
   
    g, R2 = self.polynomial_fit(xi, yi)

    gi[mid_cell] = g
    r2[mid_cell] = R2
    
    in_outlet = False
        
    while processing:
        # Verificamos si estamos al final (en un outlet)  
        if not in_outlet:
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
    
        # Verificamos si mouth_cell es un outlet
        if ixcix[mouth_cell] == 0:
            in_outlet = True
            # Si mouth_cell es un outlet, hay que eliminarla de la ventana
            # puesto que ixcix[mouth_cell] será la posición 0 del array ._ix
            win[0] = win[1]    
    
        if len(win) <= 3:
            processing = False
            
        # Obtenemos datos de elevacion y distancias para calcular pendientes
        xi = x_arr[ixcix[win]]
        yi = y_arr[ixcix[win]]
        
        g, R2 = self.polynomial_fit(xi, yi)
                    
        # Comprobamos si celda central ha sido procesada
        if gi[mid_cell] == 0:
            gi[mid_cell] = g
            r2[mid_cell] = R2
        else:
            # Si ha sido procesada, es una confluencia 
            processing = False
#             if len(confs[mid_cell]) == 0:
#                 confs[mid_cell].append(gi[mid_cell])
#                 confs[mid_cell].append(g)
#             else:
#                 confs[mid_cell].append(g)
                
# # Calculamos valores medios en confluencias                
# for cell in confs.keys():
#     if len(confs[cell]) > 0:
#         gi[cell] = np.mean(np.array(confs[cell]))
    