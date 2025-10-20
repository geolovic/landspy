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

class SwathProfile:
    def __init__(self, center_line=None, dem=None, width=0, n_lines=None, step_size=0, name=""):
        """
        Class to create a swath profile object and related parameters

        :param center_line: shapely.geometry.LineString - LineString the swath profile center line
        :param dem: landspy.DEM - Digital Elevation Model
        :param width: float - With of the swath profile (in data units)
        :param n_lines: int - Number of elevation profiles of the SWATH. The swath profile will have n_lines/2 elevation profiles at each side of the center line
        :param step_size: float - Step-size to get elevation points along the profile
        :param name: str - Name of the profile
        """

        self.name = str(name)
        self.width = float(width)
        self.center_line = center_line

        # Get step size (By default dem.cellsize if was not specified)
        if step_size == 0 or step_size < dem.getCellSize[0]:
            self.step_size = dem.getCellsize[0]
        else:
            self.step_size = step_size

        # Get number of lines (By default width/dem.cellsize)
        if n_lines == 0 or n_lines > int(width / dem.getCellSize[0]):
            self.n_lines = int(width / dem.getCellSize[0])
        else:
            self.n_lines = n_lines

        # Get distance between lines
        self.line_distance = self.width / self.n_lines

        # Get profile distances
        self.li = np.arange(0., self.center_line.length, self.step_size)

        # Creates an empty SwathProfile Object
        if center_line is None:
            return



        # Get distance between lines
        line_distance = float(width) / n_lines




        # Get the number of points for each swath line
        npoints = self.li.shape[0]

        # Create the elevation data array with the first line (baseline)
        self.data = self._get_zi(line, dem, npoints)

        # Simplify baseline
        sline = self.center_line.simplify(tolerance=dem.cellsize * 5)

        # Create the elevation data for the Swath
        for n in range(n_lines):
            dist = line_distance * (n + 1)
            left_line = sline.parallel_offset(dist, side="left")
            right_line = sline.parallel_offset(dist, side="right")
            # Sometimes parallel_offset produces MultiLineStrings
            if left_line.type == "MultiLineString":
                left_line = self._combine_multilines(left_line)
            if right_line.type == "MultiLineString":
                right_line = self._combine_multilines(right_line)

            right_line = self._flip(right_line)
            l_elev = self._get_zi(left_line, dem, npoints)
            r_elev = self._get_zi(right_line, dem, npoints)
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

    def _get_zi(self, line, dem, npoints):
        """
        Get elevations along a line in npoints equally spaced. If any point of the line falls
        outside the DEM or in a NoData cell, a np.nan value will be asigned.
        :param line : Shapely.LineString object. Input LineString
        :param dem : pRaster object. DEM with elevatations.
        :param npoints : int. Number of points along the line to get elevations
        :return zi : Numpy.ndarray. Array with size (npoints, 1) with elevations
        """
        step_size = 1.0 / npoints
        zi = []
        for idx in range(npoints):
            pt = line.interpolate(step_size * idx, normalized=True)
            xy = list(pt.coords)[0]
            z = dem.get_xy_value(xy)
            if z == dem.nodata or not z:
                z = np.nan
            zi.append(z)

        return np.array(zi, dtype="float").reshape((len(zi), 1))

    def _flip(self, line):
        """
        Flips a LineString object. Returns the new line flipped
        :param line : Shapely.LineString object. Input LineString
        :return line : Shapely.LineString object. Fliped LineString
        """
        coords = list(line.coords)
        coords = np.array(coords)[::-1]
        return LineString(coords)

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

    def draw_swath(self, ax, legend=False, drawdata=False, drawbg=False, q=False, **kwargs):
        """
        Draw the swat profile in an Axe
        :param ax : Axe where the profile will be painted. Its cleared before drawing
        :param legend : boolean. Draw the legend
        :param drawdata: boolean. Draw the data (all the profiles)
        :param drawbg: boolean. Draw a background instead of the data
        :param p : boolean. Draw the q1, q3 quartiles
        :kwargs : Colors and line widths for the profiles. See colors and linew dictionaries
        """
        ax.clear()
        colors = {"data": (0.75, 0.75, 0.75),
                  "min": (202. / 255, 111. / 255, 30. / 255),
                  "max": (169. / 255, 50. / 255, 38. / 255),
                  "mean": (22. / 255, 160. / 255, 133. / 255),
                  "q1": (0. / 255, 191. / 255, 255. / 255),
                  "q3": (0. / 255, 191. / 255, 255. / 255)}

        colors.update(kwargs)
        linew = {"dataw": 0.65, "linesw": 1.}
        linew.update(kwargs)

        if drawbg:
            poly = mpatches.Polygon(self.bg_dat, facecolor="0.85")
            ax.add_patch(poly)

        if drawdata:
            for n in range(self.data.shape[1]):
                ax.plot(self.li, self.data[:, n], lw=linew["dataw"], color=colors["data"])

        ax.plot(self.li, self.maxz, lw=linew["linesw"], color=colors["max"], label="max")
        ax.plot(self.li, self.minz, lw=linew["linesw"], color=colors["min"], label="min")
        ax.plot(self.li, self.meanz, lw=linew["linesw"], color=colors["mean"], label="mean")

        if q:
            ax.plot(self.li, self.q1, lw=linew["linesw"], color=colors["q1"], label="q1")
            ax.plot(self.li, self.q3, lw=linew["linesw"], color=colors["q3"], label="q3")

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

