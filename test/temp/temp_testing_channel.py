# -*- coding: utf-8 -*-
"""
Editor de Spyder

Este es un archivo temporal
"""

from topopy import Network, Channel, Flow, BNetwork
import ogr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cmap


canales = np.load("../data/out/canales-b3_jebja30.npy", allow_pickle=True)

canal = canales[34]

fig = plt.figure()
ax = fig.add_subplot()
d, z = canal.get_d(), canal.get_z()
ax.plot(d, z)

fig = plt.figure()
ax = fig.add_subplot()
chi, z = canal.get_chi(head=False), canal.get_z(head=False)
ax.plot(chi, z)

fig = plt.figure()
ax = fig.add_subplot()
a, slp = canal.get_a(cells=False), canal.get_slope()
ax.plot(a, slp, "b.")
ax.set_xscale("log")
ax.set_yscale("log")

fig = plt.figure()
ax = fig.add_subplot()
ksn, d = canal.get_ksn(), canal.get_d(tohead=True)
ax.plot(d, ksn)
