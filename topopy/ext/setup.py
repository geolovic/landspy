#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 11:42:41 2018

@author: vicen
"""
import numpy
from distutils.core import setup
from Cython.Build import cythonize

setup(
    ext_modules = cythonize("distance.pyx"),
    include_dirs = numpy.get_include())