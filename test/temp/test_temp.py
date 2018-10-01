# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 12:14:38 2018

@author: VicenPC
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
# Add to the path code folder and data folder
sys.path.append("../")
from topopy import Grid, Flow, Network

# Get Flow and Network
fd = Flow("data/fd_tunez.tif")    
st = Network(fd, 1000)

# Get 5 random points inside the Grid
rrows = []
rcols = []
for n in range(5):
    rrows.append(np.random.randint(0, fd.get_dims()[0]-1))
    rcols.append(np.random.randint(0, fd.get_dims()[1]-1))
    
fig, ax = plt.subplots()

srows, scols = st.snap_points(rrows, rcols, "outlets")
streams = st.get_streams()
streams.plot(ax)
ax.plot(rcols, rrows, "go")
ax.plot(scols, srows, "ro")

    
