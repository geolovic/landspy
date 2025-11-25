# -*- coding: utf-8 -*-

# swath.py
# Jose Vicente Perez Pena
# Dpto. Geodinamica-Universidad de Granada
# 18071 Granada, Spain
# vperez@ugr.es // geolovic@gmail.com
#
# MIT License (see LICENSE file)
# Version: 1.0
# 15 October, 2025
#
# Last modified 15 October, 2025

import numpy as np
import os
import shapely
from . import DEM
from shapely import LineString
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

class SwathProfile:
    def __init__(self, center_line=None, dem=None, width=0, n_lines=0, step_size=0, name=""):
        """
        Class to create a swath profile object and related parameters

        :param center_line: shapely.geometry.LineString - LineString the swath profile center line
        :param dem: landspy.DEM - Digital Elevation Model
        :param width: float - Half width of the swath profile (in data units)
        :param n_lines: int - Number of elevation profiles of the SWATH at each side of center line
        :param step_size: float - Step-size to get elevation points along the profile
        :param name: str - Name of the profile
        """
        # Creates an empty SwathProfile Object
        if center_line is None:
            self.name = ""
            self.width = 0
            self.step_size = 0
            self.n_lines = 0
            self.center_line = LineString()
            self.data = np.array([])
            return

        elif type(center_line) == str:
            # Loads the Swath profile
            # try:
            self.load_swath(center_line)
            return
            # except:
            #     raise LoadError("Error opening the Geotiff")

        self.name = str(name)
        self.width = float(width)
        self.center_line = center_line

        # Get step size (By default dem.getCellsize if was not specified)
        if step_size == 0 or step_size < dem.getCellSize()[0]:
            self.step_size = dem.getCellSize()[0]
        else:
            self.step_size = step_size

        # Get number of lines (By default width/dem.getCellsize)
        if n_lines == 0 or n_lines > int(width / dem.getCellSize()[0]):
            self.n_lines = int(width / dem.getCellSize()[0])
        else:
            self.n_lines = n_lines

        # Get distance between lines
        line_distance = self.width / self.n_lines

        # Get the number of points for each swath line
        n_points = int(self.center_line.length // self.step_size) + 1

        # Create the elevation data array with the first line (baseline)
        self.data = self._get_zi(self.center_line, dem, n_points)

        # Simplify baseline
        sline = self.center_line.simplify(tolerance=dem.getCellSize()[0] * 5)

        # Create the elevation data for the Swath
        for n in range(self.n_lines):
            dist = line_distance * (n + 1)
            left_line = sline.parallel_offset(dist, side="left")
            right_line = sline.parallel_offset(dist, side="right")
            # Sometimes parallel_offset produces MultiLineStrings
            if left_line.geom_type == "MultiLineString":
                left_line = self._combine_multilines(left_line)
            if right_line.geom_type == "MultiLineString":
                right_line = self._combine_multilines(right_line)


            l_elev_data = self._get_zi(left_line, dem, n_points)
            r_elev_data = self._get_zi(right_line, dem, n_points)
            self.data = np.append(self.data, r_elev_data, axis=1)
            self.data = np.append(self.data, l_elev_data, axis=1)

        self._get_parameters()

    def set_name(self, name=""):
        self.name = str(name)

    def get_name(self):
        return self.name

    def _get_zi(self, line, dem, n_points):
        """
        Get elevations along a line in n_points equally spaced. If any point of the line falls
        outside the DEM or in a NoData cell, a np.nan value will be assigned.
        :param line : Shapely.LineString object. Input LineString
        :param dem : pRaster object. DEM with elevations.
        :param n_points : int. Number of points along the line to get elevations
        :return zi : Numpy.ndarray. Array with size (n_points, 1) with elevations
        """
        step_size = 1.0 / n_points
        elev_data = []

        for idx in range(n_points):
            pt = line.interpolate(step_size * idx, normalized=True)
            x, y = list(pt.coords)[0]
            if not dem.isInside(x, y):
                z = np.nan
                elev_data.append((x, y, z))
                continue
            row, col = dem.xyToCell(x, y)
            z = dem.getValue(row, col)
            if z == dem.getNodata() or not z:
                z = np.nan
            elev_data.append((x, y, z))

        return np.array(elev_data, dtype="float")

    def _combine_multilines(self, line):
        """
        Combines all the parts of a MultiLineString in a single LineString
        :param line : Shapely.LineString object. Input MultiLineString
        :return line : Shapely.LineString object. Ouput LineString
        """
        xyarr = np.array([], dtype="float32").reshape((0, 2))
        for n in range(len(line.geoms)):
            xyarr = np.append(xyarr, np.array(line.geoms[n].coords), axis=0)
        return LineString(xyarr)

    def draw_swath(self, ax, q1=False, q3=False, max=False, min=False, mean=False, central=True, data='RAW', legend=False, styles=None):
        """
        Draw the swat profile in an matplotlib Axe object
        :param ax : Axe object where the profile will be painted. Its cleared before drawing
        :param q1 : boolean. Draw Q1 profile
        :param q3 : boolean. Draw Q3 profile
        :param max : boolean. Draw maximum elevation and Q3 profiles
        :param min : boolean. Draw minimum elevation profile
        :param mean : boolean. Draw mean elevation profile
        :param central : boolean. Draw central line profile (input line).
        :param data: str. String to select raw data draw mode. 'RAW' draw all profiles, 'POLYGON' draw only boundary polygon, 'NONE' does not draw raw data
        :param legend: boolean. Show the legend.
        :kwargs : Dicctionary with line styles (linewidth - linestyle - color)
        """
        ax.clear()
        base_styles = {"q1": {'lw': 1.5, 'ls': '-', 'color': (0., 0.75, 1.)},
                       "q3": {'lw': 1.5, 'ls': '-', 'color': (0., 0.75, 1.)},
                       "max": {'lw': 1.5, 'ls': '-', 'color': (1., 0., 0.)},
                       "min": {'lw': 1.5, 'ls': '-', 'color': (0., 0., 1.)},
                       "mean": {'lw': 1.5, 'ls': '-', 'color': (0.93, 0.64, 0.)},
                       "central": {'lw': 1.5, 'ls': '-', 'color': 'k'},
                       "data": {'lw': 0.75, 'ls': '-', 'color': '0.6'}}

        if styles:
            base_styles.update(styles)
        styles = base_styles

        if self.center_line.length == 0:
            return

        # Draw raw data
        if data == 'RAW':
            for n in range(self.data.shape[1]):
                ax.plot(self.li, self.data[:, n], lw=styles['data']['lw'], ls=styles['data']['ls'], color=styles['data']['color'])
        elif data == 'POLYGON':
            poly = mpatches.Polygon(self.bg_dat, facecolor=styles['data']['color'])
            ax.add_patch(poly)
            drawdata = False

        # Draw q1, q3, max, min, mean and central line
        if q1:
            ax.plot(self.li, self.q1, lw=styles['q1']['lw'], ls=styles['q1']['ls'], color=styles['q1']['color'], label="Q1")
        if q3:
            ax.plot(self.li, self.q3, lw=styles['q3']['lw'], ls=styles['q3']['ls'], color=styles['q3']['color'], label="Q3")
        if max:
            ax.plot(self.li, self.maxz, lw=styles['max']['lw'], ls=styles['max']['ls'], color=styles['max']['color'], label="max")
        if min:
            ax.plot(self.li, self.minz, lw=styles['min']['lw'], ls=styles['min']['ls'], color=styles['min']['color'], label="min")
        if mean:
            ax.plot(self.li, self.meanz, lw=styles['mean']['lw'], ls=styles['mean']['ls'], color=styles['mean']['color'], label="mean")
        if central:
            ax.plot(self.li, self.data[:, 2::3][:, 0], lw=styles['central']['lw'], ls=styles['central']['ls'], color=styles['central']['color'], label="central line")

        ax.set_xlabel("Distance [m]")
        ax.set_ylabel("Elevation [m]")

        ax.set_title(self.name)

        # QGIS Adjustment (to make the graphic nicer)
        dz = (self.maxz.max() - self.minz.min()) * 0.05
        ax.set_xlim(0, self.length)
        ax.set_ylim(self.minz.min() - dz, self.maxz.max() + dz)

        if legend:
            legend = ax.legend()
            for tx in legend.texts:
                tx.set_fontsize(10)

    def draw_thi(self, ax, enhanced=False):
        """
        Draws the THI profile in an input Axe

        :param ax : matplotlib.Axe object to draw the THI profile
        :param enhanced : boolean. Specify if the enhanced THI (THI*) is calculated
        """
        if self.center_line.length == 0:
            return

        if enhanced:
            hi = (self.HI - 0.2) / 0.6
        else:
            hi = self.HI

        max_relief = float(np.nanmax(self.relief))
        wi = 0.2 * np.log(self.relief / max_relief) + 1
        thi = (hi - 0.5) * wi + 0.5

        ax.plot(self.li, thi, c="k", linewidth=1.2)

        ax.plot([0, self.length], [0.5, 0.5], linestyle="--",
                linewidth=0.75, color=(0.4, 0.4, 0.4))

        ax.set_ylim((0.0, 1.0))
        ax.set_xlim((0.0, self.length))
        ax.set_xlabel("Distance [m]")

        if enhanced:
            label = "THI*"
        else:
            label = "THI"

        ax.set_ylabel(label)
        ax.set_yticks((0.0, 0.5, 1.0))

    def save_swath(self, path):
        """
        Save the SWATH profile in two files (*.dat and *.npy)
        Args:
            path: Path to save the SWATH profile

        """
        # Get the base name
        base_name = os.path.splitext(path)[0]

        # Save swath data as numpy array
        np.save(base_name + ".npy", self.data)

        # Save properties in text format
        fw = open(base_name + ".swt", "w")
        fw.write(self.name + "\n")
        fw.write(str(self.width) + "\n")
        fw.write(str(self.step_size) + "\n")
        fw.write(str(self.n_lines) + "\n")
        fw.write(shapely.to_wkt(self.center_line))
        fw.close()

    def load_swath(self, path):
        """
        Load the SWATH profile from disk
        """
        # Get the base name
        base_name = os.path.splitext(path)[0]

        if not os.path.exists(base_name + ".npy") or not os.path.exists(base_name + ".swt"):
            return

        # Loads swath data as numpy array
        self.data = np.load(base_name + ".npy")

        # Loads properties from the *.swt file
        fr = open(path, "r")
        lines = fr.readlines()
        self.name = lines[0].strip()
        self.width = float(lines[1].strip())
        self.step_size = float(lines[2].strip())
        self.n_lines = int(lines[3].strip())
        self.center_line = shapely.from_wkt(lines[4].strip())
        fr.close()

        self._get_parameters()

    def _get_parameters(self):

        if self.center_line.length == 0:
            return

        # Get li
        self.li = np.arange(0., self.center_line.length, self.step_size)

        # View of the elevation data to get max, min, mean and q1-q3
        zdat = self.data[:, 2::3]
        self.idx_max = np.nanargmax(zdat, axis=1)
        self.idx_min = np.nanargmin(zdat, axis=1)
        q1 = np.nanquantile(zdat, 0.25, axis=1)
        q3 = np.nanquantile(zdat, 0.75, axis=1)
        self.idx_q1 = np.nanargmin(np.abs(zdat - q1[:, None]), axis=1)
        self.idx_q3 = np.nanargmin(np.abs(zdat - q3[:, None]), axis=1)

        # Get parameters (max, min, mean, q1, q3, HI, relief)
        self.maxz = zdat[np.arange(zdat.shape[0]), self.idx_max]
        self.minz = zdat[np.arange(zdat.shape[0]), self.idx_min]
        self.meanz = np.nanmean(zdat, axis=1)
        self.q1 = zdat[np.arange(zdat.shape[0]), self.idx_q1]
        self.q3 = zdat[np.arange(zdat.shape[0]), self.idx_q3]
        self.HI = (self.meanz - self.minz) / (self.maxz - self.minz)
        self.relief = self.maxz - self.minz

        # Get a background polygon for the data
        xi = np.append(self.li, self.li[::-1])
        yi = np.append(self.maxz, self.minz[::-1])
        xi = xi.reshape((xi.size, 1))
        yi = yi.reshape((yi.size, 1))
        self.bg_dat = np.append(xi, yi, axis=1)

        # Length of the swath
        self.length = self.li[-1]
