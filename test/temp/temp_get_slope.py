# -*- coding: utf-8 -*-
"""
Editor de Spyder

Este es un archivo temporal
"""

import numpy as np
import matplotlib.pyplot as plt
from topopy import DEM, Flow, Network

indem = "../data/in/tunez.tif"

dem = DEM(indem)
fd = Flow(dem)
net = Network(dem, fd)

# Get ixcix auxiliar array
ixcix = np.zeros(net._ncells, np.int)
ixcix[net._ix] = np.arange(net._ix.size)

# Get heads and sorted them by elevation
heads = net.get_stream_poi("heads", "IND")
elev = net._zx[ixcix[heads]]
spos = np.argsort(-elev)
heads = heads[spos]
winlen = 11
slopes = np.zeros(net._ncells)

# Taking sequentally all the heads and compute downstream flow
for head in heads:
    processing = True
    lcell = head
    mcell = net._ixc[ixcix[head]]
    fcell = net._ixc[ixcix[mcell]]

    if ixcix[fcell] == 0 or ixcix[net._ixc[ixcix[fcell]]] == 0:
        continue
    
    win = [fcell, mcell, lcell]
    zi = net._zx[ixcix[win]]
    di = net._dx[ixcix[win]]
    slp = np.polyfit(di, zi, 1)[0]
    slopes[mcell] = slp
    
    while processing:
        # Cogemos la siguiente celda del canal (next_cell)
        fcell = win[0]
        next_cell = net._ixc[ixcix[fcell]]
        
        # Comprobamos la siguiente celda del canal
        # Si la siguiente celda del canal no es cero
        if ixcix[next_cell] != 0:
            # Añadimos siguiente celda
            win.insert(0, next_cell)
            fcell = win[0]
            # Movemos celda central
            mcell = net._ixc[ixcix[mcell]]

            if len(win) < winlen:
                # Si estamos al principio del canal, win crece
                next_cell = net._ixc[ixcix[fcell]]
                win.insert(0, next_cell)
                fcell = win[0]
            else:
                # Si no estamos al principio, eliminamos celda final
                win.pop()
        # Si la siguiente celda es cero, no añadimos celdas, sino que win decrece
        else:
            mcell = net._ixc[ixcix[mcell]]
            win.pop()
            win.pop()
            if len(win) == 3:
                processing = False
                
        # Calculamos la pendiente de la celda central
        zi = net._zx[ixcix[win]]
        di = net._dx[ixcix[win]]
        slp = np.polyfit(di, zi, 1, full=True)[0][0]
        
        if slopes[mcell] == 0:
            slopes[mcell] = slp
        else:
            processing = False
        

#
#elev = net._zx[ixcix[heads]]

#
#demarr = dem.read_array().ravel()
#fig, ax = plt.subplots()
## Get heads
#heads = net.get_stream_poi("heads", "IND")
#
#visited = np.zeros(net._ncells, np.bool)
#
##for head in heads:
#
#head = 20860
#
#chcell = head
#ixcix = np.zeros(net._ncells, np.int)
#ixcix[net._ix] = np.arange(net._ix.size)
#
#elev = net._zx[ixcix[heads]]
#
## Sort heads by elevation
#sorted_pos = np.argsort(-elev)
#heads = heads[sorted_pos]
#
#
#for head in heads:
#    chcell = head
#    win_points = [chcell]
#    while ixcix[chcell] != 0:
#        chcell = net._ixc[ixcix[chcell]]
#        elev = net._zx[ixcix[chcell]]
#        print(chcell, elev, demarr[chcell])
#        if not visited[chcell]:
#            visited[chcell] = True
#        else:
#            win_points.append(chcell)
#            break
#        win_points.append(chcell)
#    
#    row, col = net.ind_2_cell(win_points)
#    x, y = net.cell_2_xy(row, col)
#    
#        
#    ax.plot(x, y)
