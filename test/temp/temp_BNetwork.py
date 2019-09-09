# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 15:06:28 2019

@author: Usuario
"""
from topopy import DEM, Flow, Network, BNetwork
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker

dem = DEM("../data/in/morocco.tif")
fld = Flow(dem)
basin = fld.get_drainage_basins([580188.393, 508786.190])
net = Network(fld)

bnet = BNetwork(net, basin)
canales = bnet.get_streams()

#path = "../data/out/bnet.net"
##bnet.save(path)
#
#bnet = BNetwork(path)
#fig = plt.figure()
#for n, value in enumerate(np.arange(0.2, 0.6, 0.05)):
#    bnet.calculate_chi(value)
#    ax = fig.add_subplot(2, 4, n+1)
#    bnet.chi_plot(ax)

def chi_analysis(self, start=0, stop=1, step=0.01, graphic=True):
    """
    This function perform a Chi analysis by using different thetaref values. It returns
    the thetaref that yield the minimum dispersion (lower R2) for all basin channels.
    
    Parameters:
    ============
    start : First thetaref value to try
    stop : Last thetaref value to try
    step : Step to try different thetaref values
    graphic : *bool* If true, this function will draw a graphic with the results. 
    
    Return:
    ============
    (best_theta, R2) : Tuple with the best theta and its R2 coefficient. 
    """

    desv = []
    
    for thetaref in np.arange(start, stop, step):
        bnet.calculate_chi(thetaref)
        chi = bnet._chi
        zi = bnet._zx
        poli, SCR = np.polyfit(chi, zi, deg = 1, full = True)[:2]
        r2 = float(1 - SCR/(zi.size * zi.var()))
        desv.append((thetaref, r2))
        
    desv = np.array(desv)
    
    if graphic:
        fig, ax = plt.subplots()
        
        
        ax.set_xlabel("$m/n$")
        ax.set_ylabel("$R^2$")
        
        ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
        ax.yaxis.set_major_locator(ticker.MultipleLocator(0.1))
        
        # Get maximum R2
        pos = np.argmax(desv[:,1])
        
        mindesv = desv[pos, 0]
        mintheta = desv[pos, 1]
        miny = np.floor(np.min(desv[:,1]) * 10) / 10
        
        ax.set_xlim(xmin=0)
        ax.set_ylim(ymin=miny)
        
        ax.plot([0, mindesv, mindesv], [mintheta, mintheta, miny], linestyle="-", c="orange", lw=0.8)
        ax.grid(True, which='major', axis="both", linestyle ="--", c="0.7", lw=0.5)
        ax.plot(desv[:,0], desv[:, 1])
        
    return (mindesv, desv[pos,1])
        
        
chi_analysis(bnet)