# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 17:36:31 2018

@author: Usuario
"""
from topopy import TProfile, BNetwork, canal_2_profile
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import ogr, osr
import os

basedir = "/Users/vicen/Desktop/tunez/gisdata/BNetworks/out_files"
perfiles = np.load(basedir + "/perfiles.npy")

class NProfiler:
    
    def __init__(self, di, zi, kps=[], name = ""):
        
        ddi = np.linspace(di.min(), di.max(), 1024)
        zzi = np.interp(ddi, di, zi)
        self.di = (ddi - ddi.min()) / (ddi.max() - ddi.min())
        self.zi = (zzi - zzi.min()) / (zzi.max() - zzi.min())
        self.name = name
        
        if len(kps) > 0:
            kp_x = (kps[:, 0] - ddi.min()) / (ddi.max() - ddi.min())
            kp_y = (kps[:, 1] - zzi.min()) / (zzi.max() - zzi.min())
            self.kps = np.array((kp_x, kp_y)).T
        else:
            self.kps = np.array([])
        
        self.calculate_parameters()
        
    def calculate_parameters(self):
        
        area = 0
        for n in range(self.di.size - 1):
            base = self.di[n + 1] - self.di[n]
            z1 = min([self.zi[n + 1], self.zi[n]])
            z2 = max([self.zi[n + 1], self.zi[n]])
            area += base * z1
            area += (base * (z2-z1)) / 2
        
        area = 0.5 - area
        
        if area > 0:
            self.ct = (area/0.5) * 100
        else:
            self.ct = (area/0.5) * 100
        
        linzi = np.interp(self.di, np.array([0, 0.5, 1]), np.array([1, 0.5, 0]))
        pos = np.argmax(linzi - self.zi)
        self.pos = pos
        self.cmax = linzi[pos] - self.zi[pos]
        self.lmax = self.di[pos]
    
    def plot(self, ax=None):
        
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            
        ax.plot(self.di, self.zi)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        
rids = [1, 24, 28, 0]
pperfiles = []    
for perfil in perfiles:
    if perfil.rid in rids:
        pperfiles.append(perfil)
        
fig = plt.figure(figsize = (8.3, 11.7))
axs = fig.subplots(1, 4)

for n, ax in enumerate(axs):
    perfil = pperfiles[n]
    kp = perfil._knickpoints
    if len(kp) > 0:
        xi = perfil.get_l()
        yi = perfil.get_z()
        kpxi = []
        kpyi = []
        for k in kp:
            kpxi.append(xi[k[0]])
            kpyi.append(yi[k[0]])
            
        kps = np.array((kpxi, kpyi)).T
        nPerfil = NProfiler(perfil.get_l(), perfil.get_z(), kps = kps, name=str(perfil.rid))
    else:
        nPerfil = NProfiler(perfil.get_l(), perfil.get_z(), name=str(perfil.rid))
    
    ax.plot(nPerfil.di, nPerfil.zi)
    ax.plot([0, 1], [1, 0], c="red", lw=0.7)
    if len(nPerfil.kps) > 0:
        ax.plot(nPerfil.kps[:, 0], nPerfil.kps[:,1], "ro", ms=4)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_title(nPerfil.name)
    ax.set_aspect("equal")
    s = "C: {:.2f}%\n".format(nPerfil.ct)
    s += "MaxC: {0:.2f}\n".format(nPerfil.cmax)
    s += "MaxL: {0:.2f}\n".format(nPerfil.lmax)
    ax.text(0.96, 0.96, s, va = "top", ha = "right", fontsize=8)
    ax.set_xticks([0, 1])
    ax.set_yticks([0, 1])
    ax.set_xticklabels(["0", "1"])

plt.tight_layout(pad=0.4)





