# -*- coding: utf-8 -*-
from topopy import DEM, Flow, Network
import numpy as np
import matplotlib.pyplot as plt

def draw_network(net, ax=None):
    if ax == None:
        fig, ax = plt.subplots()
    
    # Get ixcix auxiliar array
    ixcix = np.zeros(self._ncells, np.int)
    ixcix[self._ix] = np.arange(self._ix.size)
    
    # Get heads and confluences and sorted them by elevation
    heads = self.get_stream_poi("heads", "IND")
    confs = self.get_stream_poi("confluences", "IND")
    outls = self.get_stream_poi("outlets", "IND")
    all_p = np.append(heads, confs)
    all_p = np.append(all_p, outls)
    
    elev = self._zx[ixcix[all_p]]
    spos = np.argsort(-elev)
    all_p = all_p[spos]
    
    for point in all_p:
        if point in outls:
            continue
        cell = point
        idx = [cell]
        processing_ch = True
        while processing_ch:
            next_cell = net._ixc[ixcix[cell]]
            idx.append(next_cell)
            if next_cell in all_p:
                processing_ch = False
            cell = next_cell
        
        row, col = net.ind_2_cell(idx)
        x, y = net.cell_2_xy(row, col)
        ax.plot(x, y, color="b")
    
def draw_point(point, net, ax=None, color="r"):
    if ax == None:
        fig, ax = plt.subplots()   
    
    row, col = net.ind_2_cell(point)
    x, y = net.cell_2_xy(row, col)
    ax.plot(x, y, ls="None", marker="o", mfc=color, ms=5)


dem = DEM("../data/in/morocco.tif")
flw = Flow(dem)
net = Network(flw, 1029)

self = net
npoints = 5
kind = "slp"

#fig, ax = plt.subplots()
#draw_network(net, ax)

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

# Get arrays depending on type
if kind == 'slp':
    x_arr = self._dx
    y_arr = self._zx
elif kind == 'ksn':
    x_arr = self._chi
    y_arr = self._zx

gi = np.zeros(self._ncells)
r2 = np.zeros(self._ncells)

def walk_network(self, npoints, kind):
    if kind not in ['slp', 'ksn']:
        kind = 'slp'

    # Get arrays depending on type
    if kind == 'slp':
        x_arr = self._dx
        y_arr = self._zx
    elif kind == 'ksn':
        x_arr = self._chi
        y_arr = self._zx
    
    # Preare arrays for gradients
    gi = np.zeros(self._ncells)
    r2 = np.zeros(self._ncells)
        
    # Get ixcix auxiliar array
    ixcix = np.zeros(self._ncells, np.int)
    ixcix[self._ix] = np.arange(self._ix.size)
    
    # Get heads and sorted them by elevation
    heads = self.get_stream_poi("heads", "IND")
    elev = self._zx[ixcix[heads]]
    spos = np.argsort(-elev)
    heads = heads[spos]
    
    # Get outlets
    outlets = net.get_stream_poi("outlets", "IND")
    
    # Get window size
    winlen = npoints * 2 + 1
    
    # Loop trhough all the heads
    for head in heads:
        # Hacemos ventana con 3 primeras celdas
        lcell = head
        mcell = self._ixc[ixcix[head]]
        fcell = self._ixc[ixcix[mcell]]
        win = [fcell, mcell, lcell]
    
        # Si canal tiene solo 2 o menos celdas, no es válido, se toma siguiente head
        if (mcell in outlets) or (lcell in outlets):
            continue
        # Si canal tiene 3 o 4 celdas solamente, se calculan los valores
        elif (fcell in outlets) or (self._ixc[ixcix[fcell]] in outlets):
            # Get slope of central cell by linear regression
            xi = x_arr[ixcix[win]]
            yi = y_arr[ixcix[win]]
            poli, SCR = np.polyfit(xi, yi, deg = 1, full = True)[:2]
            if yi.size * yi.var() == 0:
                R2 = 1
            else:
                R2 = float(1 - SCR/(yi.size * yi.var()))
            
            continue
            
        # Comenzamos el procesado de los canales válidos
        processing = True
        while processing:
            # Movemos celda central y tomamos celda siguiente
            mcell = self._ixc[ixcix[mcell]]
            next_cell = self._ixc[ixcix[fcell]]            
            if fcell in outlets:
                win.pop()
                win.pop()
                if len(win) <= 3:
                    processing=False                                
            elif len(win) < winlen:
                win.insert(0, next_cell)
                win.insert(0, self._ixc[ixcix[next_cell]])                 
            else:
                win.insert(0, next_cell)
                win.pop()
            fcell = win[0]        
            yield win

for win in walk_network(net, npoints):
    mcell = win[int(len(win)/2)]
    print(win, mcell)

# Taking sequentally all the heads and compute downstream flow
#for head in heads:
#    processing = True
#    lcell = head
#    mcell = self._ixc[ixcix[head]]
#    fcell = self._ixc[ixcix[mcell]]
#
#    if ixcix[fcell] == 0 or ixcix[self._ixc[ixcix[fcell]]] == 0:
#        continue
#
#    # Get slope of central cell by linear regression
#    win = [fcell, mcell, lcell]
#    xi = x_arr[ixcix[win]]
#    yi = y_arr[ixcix[win]]
#    poli, SCR = np.polyfit(xi, yi, deg = 1, full = True)[:2]
#    
#    # Fix to avoid issues with horizontal colinear points
#    if yi.size * yi.var() == 0:
#        R2 = 1
#    else:
#        R2 = float(1 - SCR/(yi.size * yi.var()))
#    
#    # Fill gradient and R2 arrays
#    g = poli[0]
#    gi[mcell] = g
#    r2[mcell] = R2
#    
#    while processing:
#        # Take the next cell of the network
#        fcell = win[0]
#        next_cell = self._ixc[ixcix[fcell]]
#        
#        # If next cell is not 0
#        if ixcix[next_cell] != 0:
#            # Add cell to window
#            win.insert(0, next_cell)
#            fcell = win[0]
#            # Move central cell
#            mcell = self._ixc[ixcix[mcell]]
#            
#            # Check window lenght
#            if len(win) < winlen:
#                # if win didn't reach winlen, win grows
#                next_cell = self._ixc[ixcix[fcell]]
#                win.insert(0, next_cell)
#                fcell = win[0]
#            else:
#                # if win reached winlen, remove last cell
#                win.pop()
#        # If next cell is 0 [we reached an outlet], we did not add cells, but win decreases
#        else:
#            mcell = self._ixc[ixcix[mcell]]
#            win.pop()
#            win.pop()
#            if len(win) == 3:
#                processing = False
#                gi[fcell] = 0.00001
#                r2[fcell] = 0.00001
#                
#        # Get slope of central cell by linear regression
#        xi = x_arr[ixcix[win]]
#        yi = y_arr[ixcix[win]]
#        poli, SCR = np.polyfit(xi, yi, deg = 1, full = True)[:2]
#        if yi.size * yi.var() == 0:
#            R2 = 1
#        else:
#            R2 = float(1 - SCR/(yi.size * yi.var()))
#
#        # Fill gradient and R2 arrays
#        g = poli[0]
#            
#        if gi[mcell] == 0: # The cell wasn't processed yet
#            gi[mcell] = g
#            r2[mcell] = R2
#        else: # The cell was already processed
#            processing = False
#            draw_point(mcell, net, ax)
#            if not mcell in confs:
#                print(head, mcell)
#            
##            if len(dconfs[mcell]) > 0:
##                dconfs[mcell].append(g)
##            else:
##                dconfs[mcell].append(gi[mcell])
##                dconfs[mcell].append(g)
#
#if kind == 'slp':
#    self._slp = gi[self._ix]
#    self._r2slp = r2[self._ix]
#    self._slp_npoints = npoints
#elif kind == 'ksn':
#    self._ksn = gi[self._ix]
#    self._r2ksn = r2[self._ix]
#    self._ksn_npoints = npoints
#    
#    
#
#
#confs = net.get_stream_poi("confluences", "IND")
#draw_point(confs, net, ax) 
        
    