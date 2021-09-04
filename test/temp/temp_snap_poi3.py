# -*- coding: utf-8 -*-
"""
Editor de Spyder

Este es un archivo temporal
"""

from topopy import Flow
import numpy as np
import matplotlib.pyplot as pl

fd = Flow("../data/in/small25_fd.tif")

puntos = np.array([[470511.3, 4116528.4, 344],
                   [470740.9, 4116920.4, 36]])

snap = fd.snap_points(puntos, 500)
cuencas =   fd.get_drainage_basins(snap)
cuencas.save("../data/out/prueba_cuencas.tif")
                
    
    