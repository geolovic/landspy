#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 11 February, 2020
Testing suite for HCurve class
@author: J. Vicente Perez
@email: geolovic@gmail.com
@date: 06 October 2022
@last_modified: 15 October 2025
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
        flw_path = "{}/{}_fd.tif".format(outfolder, "small25")
        basin_path = "{}/{}_basins.tif".format(outfolder, "small25")
        dem_path = "{}/{}.tif".format(infolder, "small25")
        dem = DEM(dem_path)
        if not os.path.exists(flw_path):
            fd = Flow(dem)
            fd.save(flw_path)
        else:
            fd = Flow(flw_path)

        if not os.path.exists(basin_path):
            basins = fd.drainageBasins(min_area=0)
            basins.save(basin_path)
        else:
            basins = Grid(basin_path)

    def test_HCurve_few_points(self):
        # Create hypsometric curves with few points
        dem_path = "{}/{}.tif".format(infolder, "small25")
        basin_path = "{}/{}_basins.tif".format(outfolder, "small25")
        dem = DEM(dem_path)
        basins = Grid(basin_path)

        # Get all basins (some of them are really small, less than 50 points)
        bids, counts = np.unique(basins.readArray(), return_counts=True)
        if 0 in bids:
            bids = bids[1:]
            counts = counts[1:]

        for n, bid in enumerate(bids):
            hcurve = HCurve(dem, basins, bid)
            # Check is HCurve instance
            self.assertIsInstance(hcurve, HCurve)
            # Check if its an empty curve
            ncells = counts[n]
            if ncells < 50:
                self.assertEqual(hcurve._data.size, 4)
                self.assertEqual(hcurve.getHI(), 0.5)
                self.assertEqual(hcurve.getHI2(), 0.5)
                self.assertEqual(hcurve.getKurtosis(), 0)
                self.assertEqual(hcurve.getSkewness(), 0)
                self.assertEqual(hcurve.getDensityKurtosis(), 0)
                self.assertEqual(hcurve.getDensitySkewness(), 0)
            else:
                self.assertNotEqual(hcurve.getHI(), 0)
                self.assertNotEqual(hcurve.getHI2(), 0)
                self.assertNotEqual(hcurve.getKurtosis(), 0)
                self.assertNotEqual(hcurve.getSkewness(), 0)
                self.assertNotEqual(hcurve.getDensityKurtosis(), 0)
                self.assertNotEqual(hcurve.getDensitySkewness(), 0)



if __name__ == "__main__":
    unittest.main()


