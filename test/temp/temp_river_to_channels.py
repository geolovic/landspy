# -*- coding: utf-8 -*-
"""
Editor de Spyder

Este es un archivo temporal
"""

from topopy import Channel
import numpy as np
import matplotlib.pyplot as pl


canales = np.load("../../apps/canales2.npy", allow_pickle=True)
canal = canales[0]
npoints = 5
kind = "ksn"

self = canal

# Get arrays depending on type
y_arr = self._zx
if kind == 'ksn':
    x_arr = self._chi
else:
    x_arr = self._dx

for n in range(self._ix.size):
    low = n - npoints
    high = n + npoints

    if low < 0:
        low = 0

    xi = x_arr[low:high + 1]
    yi = y_arr[low:high + 1]
    poli, SCR = np.polyfit(xi, yi, deg = 1, full = True)[:2]
    
    if yi.size * yi.var() == 0:
        R2 = 0.0
    else:
        R2 = float(1 - SCR/(yi.size * yi.var()))

    g = poli[0]
    r2 = abs(R2)

    if abs(g) < 0.001:
       g = 0.001
       
       
    if kind == 'ksn':
        self._ksn[n] = g
        self._r2ksn[n] = r2
        self._ksn_np = npoints
    else:
        self._slp[n] = g
        self._r2slp[n] = r2
        self._slp_np = npoints



