#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 11:42:41 2018

@author: vicen
"""

from distutils.core import setup
from Cython.Build import cythonize

setup(
    ext_modules = cythonize("c_distance.pyx"))