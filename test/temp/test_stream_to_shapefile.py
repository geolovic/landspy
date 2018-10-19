#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 17:53:42 2018

@author: vicen
"""

from topopy import Network
import numpy as np
import matplotlib.pyplot as plt

net = Network("../data/in/jebja30_network.net")

self = net
fig, ax = plt.subplots()

#def export_to_shp(self, path="", continous=False):
    
# Get channel segments and orders
ch_seg = self.get_stream_segments(False).ravel()
ch_ord = self.get_stream_order(asgrid=False).ravel()
ch_seg = ch_seg[self._ix]
ch_ord = ch_ord[self._ix]

# Get ixcix auxiliar array
ixcix = np.zeros(self._ncells, np.int)
ixcix[self._ix] = np.arange(self._ix.size)

channels = []
conexions = []
seg_ids = np.unique(ch_seg)
for idx in seg_ids:
    # skip zero channel (no channel)
    if idx == 0:
        continue
    # Get givers for the segment id
    pos = np.where(ch_seg == idx)[0]
    ch_ix = self._ix[pos]
    
    # Add last point
    first = ch_ix[0]
    last = self._ixc[ixcix[ch_ix[-1]]]
    ch_ix = np.append(ch_ix, last)
    first = ch_ix[0]
    
    # Get segment order and receiver segment
    order = ch_ord[ixcix[first]]
    if ixcix[last] == 0:
        flowto = 0
    else:
        flowto = ch_seg[ixcix[last]]
    channels.append([idx, ch_ix])
    conexions.append([idx, flowto, order])
    
    
    colors = {1:"blue", 2:"red", 3:"green"}
    row, col = net.ind_2_cell(ch_ix)
    xi, yi = net.cell_2_xy(row, col)
    ax.plot(xi, yi, color = colors[order])

    ax.text(xi[0], yi[0], str(idx), size=12)
    
procesados = []
canales = []
conexions = np.array(conexions)

for conex in conexions:
    canal = conex[0]
    orden = conex[2]
    
    cana
    terminado = False
    while not terminado:
        next_canal = ch[1]
        next_orden = conexions np.where(conexions[0] == next_canal)
        if next_canal == 0:
            terminado=True
            
            


    