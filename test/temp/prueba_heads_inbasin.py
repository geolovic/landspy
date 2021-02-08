#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  7 18:53:07 2021

@author: vicen
"""

from topopy import Basin
import numpy as np

basin = Basin("../data/in/jebja30_sbasin.tif")

xmin, xmax, ymin, ymax = basin.get_extent()

xmin -= 4000
ymin -= 4000
xmax += 4000
ymax += 4000

xi = np.random.randint(xmin, xmax, 100)
yi = np.random.randint(ymin, ymax, 100)

np.savetxt("../data/out/puntos_total.txt", np.array((xi, yi)).T, delimiter=";", header="X;Y\n",
           comments="")


# puntos dentro
inside = basin.is_inside(xi, yi, True)

x1 = xi[inside]
y1 = yi[inside]

np.savetxt("../data/out/puntos_dentro.txt", np.array((x1, y1)).T, delimiter=";", header="X;Y\n",
           comments="")

# puntos rectangulo
inside = basin.is_inside(xi, yi, False)

x2 = xi[inside]
y2 = yi[inside]
np.savetxt("../data/out/puntos_rect.txt", np.array((x2, y2)).T, delimiter=";", header="X;Y\n",
           comments="")

