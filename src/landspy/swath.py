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
            return

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
        self.line_distance = self.width / self.n_lines

        # Get profile distances for the center line (these will be x coordinates of the swath profile)
        self.li = np.arange(0., self.center_line.length, self.step_size)

        # Get the number of points for each swath line
        self.npoints = self.li.shape[0]

        # Create the elevation data array with the first line (baseline)
        self.data = self._get_zi(self.center_line, dem, self.npoints)

        # Simplify baseline
        sline = self.center_line.simplify(tolerance=dem.getCellSize()[0] * 5)
        self.lines = []

        # Create the elevation data for the Swath
        for n in range(self.n_lines):
            dist = self.line_distance * (n + 1)
            left_line = sline.parallel_offset(dist, side="left")
            right_line = sline.parallel_offset(dist, side="right")
            # Sometimes parallel_offset produces MultiLineStrings
            if left_line.geom_type == "MultiLineString":
                left_line = self._combine_multilines(left_line)
            if right_line.geom_type == "MultiLineString":
                right_line = self._combine_multilines(right_line)
            self.lines.append(left_line)
            self.lines.append(right_line)


            l_elev = self._get_zi(left_line, dem, self.npoints)
            r_elev = self._get_zi(right_line, dem, self.npoints)
            self.data = np.append(self.data, r_elev, axis=1)
            self.data = np.append(self.data, l_elev, axis=1)

        # Get parameters (max, min, mean, q1, q3, HI, relief)
        self.maxz = np.nanmax(self.data, axis=1)
        self.minz = np.nanmin(self.data, axis=1)
        self.meanz = np.nanmean(self.data, axis=1)
        self.q1 = np.nanpercentile(self.data, q=25, axis=1)
        self.q3 = np.nanpercentile(self.data, q=75, axis=1)
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

    def set_name(self, name=""):
        self.name = str(name)

    def get_name(self):
        return self.name

    def _get_zi(self, line, dem, npoints):
        """
        Get elevations along a line in npoints equally spaced. If any point of the line falls
        outside the DEM or in a NoData cell, a np.nan value will be assigned.
        :param line : Shapely.LineString object. Input LineString
        :param dem : pRaster object. DEM with elevations.
        :param npoints : int. Number of points along the line to get elevations
        :return zi : Numpy.ndarray. Array with size (npoints, 1) with elevations
        """
        step_size = 1.0 / npoints
        zi = []
        for idx in range(npoints):
            pt = line.interpolate(step_size * idx, normalized=True)
            x, y = list(pt.coords)[0]
            if not dem.isInside(x, y):
                z = np.nan
                zi.append(z)
                continue
            row, col = dem.xyToCell(x, y)
            z = dem.getValue(row, col)
            if z == dem.getNodata() or not z:
                z = np.nan
            zi.append(z)

        return np.array(zi, dtype="float").reshape((len(zi), 1))

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

        # Draw raw data
        if data == 'RAW':
            for n in range(self.data.shape[1]):
                ax.plot(self.li, self.data[:, n], lw=styles['data']['lw'], ls=styles['data']['ls'], color=styles['data']['color'])
        elif data == 'POLYGON':
            poly = mpatches.Polygon(self.bg_dat, facecolor="0.85")
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
            ax.plot(self.li, self.data[:, 0], lw=styles['central']['lw'], ls=styles['central']['ls'], color=styles['central']['color'], label="central line")

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
                tx.set_fontsize(12)

    def draw_thi(self, ax, enhanced=False):
        """
        Draws the THI profile in an input Axe

        :param ax : matplotlib.Axe object to draw the THI profile
        :param enhanced : boolean. Specify if the enhanced THI (THI*) is calculated
        """
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

