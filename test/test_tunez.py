# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 17:36:31 2018

@author: Usuario
"""

from topopy import DEM, Flow, Network, BNetwork
import numpy as np
import matplotlib.pyplot as plt


bnet = BNetwork("data/in/ttunez_net")

#bnet.calculate_distances()
#
#canal = bnet.get_main_channel()
