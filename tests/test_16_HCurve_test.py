#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 11 February, 2020
Testing suite for HCurve class
@author: J. Vicente Perez
@email: geolovic@gmail.com
@date: 01 October 2022
@last_modified: 15 october, 2025
"""

import unittest
import numpy as np
from landspy import DEM, Basin, Grid, Flow, HCurve

import sys, os
# Forzar el directorio actual al del archivo
os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.getcwd())
infolder = "data/in"
outfolder = "data/out"

class HCurveTest(unittest.TestCase):

    def setUp(self):
        # Create the datasets (if do not exists)
        self.files = ["jebja30", "tunez"]
        for file in self.files:
            flw_path = "{}/{}_fd.tif".format(outfolder, file)
            basin_path = "{}/{}_basins.tif".format(outfolder, file)
            dem_path = "{}/{}.tif".format(infolder, file)
            dem = DEM(dem_path)
            if not os.path.exists(flw_path):
                fd = Flow(dem)
                fd.save(flw_path)
            else:
                fd = Flow(flw_path)

            if not os.path.exists(basin_path):
                basins = fd.drainageBasins(min_area=0.0025)
                basins.save(basin_path)
            else:
                basins = Grid(basin_path)

    def test_HCurve_01(self):
        # Create hypsometric curves from a DEM and a basin Grid
        for file in self.files:
            dem_path = "{}/{}.tif".format(infolder, file)
            basin_path = "{}/{}_basins.tif".format(outfolder, file)
            dem = DEM(dem_path)
            basins = Grid(basin_path)

            # Obtenemos la cuenca mayor y una aleatoria
            bids, counts = np.unique(basins.readArray(), return_counts=True)
            if 0 in bids:
                bids = bids[1:]
                counts = counts[1:]

            idmax = bids[np.argmax(counts)]
            idrdn = np.random.choice(bids)
            
            bids = [idmax, idrdn]
            lbls = ["max", "random"]

            for n, bid in enumerate(bids):
                hcurve = HCurve(dem, basins, bid)
                # Check is HCurve instance
                self.assertIsInstance(hcurve, HCurve)
                # Check all moments
                self.assertNotEqual(hcurve.getHI(), 0)
                self.assertNotEqual(hcurve.getHI2(), 0)
                self.assertNotEqual(hcurve.getKurtosis(), 0)
                self.assertNotEqual(hcurve.getSkewness(), 0)
                self.assertNotEqual(hcurve.getDensityKurtosis(), 0)
                self.assertNotEqual(hcurve.getDensitySkewness(), 0)

    def test_HCurve_01b(self):
        # Create hipsometric curves from a DEM and a path to a basins tiff file
        for file in self.files:
            dem_path = "{}/{}.tif".format(infolder, file)
            basin_path = "{}/{}_basins.tif".format(outfolder, file)
            dem = DEM(dem_path)
            basins = Grid(basin_path)

            # Obtenemos la cuenca mayor y una aleatoria
            bids, counts = np.unique(basins.readArray(), return_counts=True)
            if 0 in bids:
                bids = bids[1:]
                counts = counts[1:]

            idmax = bids[np.argmax(counts)]
            idrdn = np.random.choice(bids)

            bids = [idmax, idrdn]
            lbls = ["max", "random"]

            for n, bid in enumerate(bids):
                hcurve = HCurve(dem, basin_path, bid)
                # Check is HCurve instance
                self.assertIsInstance(hcurve, HCurve)
                # Check all moments
                self.assertNotEqual(hcurve.getHI(), 0)
                self.assertNotEqual(hcurve.getHI2(), 0)
                self.assertNotEqual(hcurve.getKurtosis(), 0)
                self.assertNotEqual(hcurve.getSkewness(), 0)
                self.assertNotEqual(hcurve.getDensityKurtosis(), 0)
                self.assertNotEqual(hcurve.getDensitySkewness(), 0)

    def test_HCurve_02(self):
        # Create hipsometric curves from a Basin
        for file in self.files:
            dem_path = "{}/{}.tif".format(infolder, file)
            basin_path = "{}/{}_basins.tif".format(outfolder, file)
            dem = DEM(dem_path)
            basins = Grid(basin_path)

            # Obtenemos la cuenca mayor
            bids, counts = np.unique(basins.readArray(), return_counts=True)
            if 0 in bids:
                bids = bids[1:]
                counts = counts[1:]

            idmax = bids[np.argmax(counts)]

            # Creamos objeto Cuenca (Basin)
            basin = Basin(dem, basins, idmax)

            # Creamos HCurve
            hcurve = HCurve(basin)
            # Check is HCurve instance
            self.assertIsInstance(hcurve, HCurve)
            # Check all moments
            self.assertNotEqual(hcurve.getHI(), 0)
            self.assertNotEqual(hcurve.getHI2(), 0)
            self.assertNotEqual(hcurve.getKurtosis(), 0)
            self.assertNotEqual(hcurve.getSkewness(), 0)
            self.assertNotEqual(hcurve.getDensityKurtosis(), 0)
            self.assertNotEqual(hcurve.getDensitySkewness(), 0)

    def test_HCurve_03(self):
        # Create hipsometric curve, save and load it
        for file in self.files:
            dem_path = "{}/{}.tif".format(infolder, file)
            basin_path = "{}/{}_basins.tif".format(outfolder, file)
            dem = DEM(dem_path)
            basins = Grid(basin_path)

            # Obtain biggest basin (id)
            bids, counts = np.unique(basins.readArray(), return_counts=True)
            if 0 in bids:
                bids = bids[1:]
                counts = counts[1:]

            idmax = bids[np.argmax(counts)]

            # Create hypsometric curve
            hcurve = HCurve(dem, basins, idmax)

            # Save it and load it
            hcurve.save(outfolder + "/temp_curve.txt")
            curva2 = HCurve(outfolder + "/temp_curve.txt")

            # Compare both curves
            self.assertEqual(hcurve.moments, curva2.moments)
            self.assertEqual(np.array_equal(hcurve._data, curva2._data), True)

if __name__ == "__main__":
    unittest.main()