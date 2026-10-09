# -*- coding: utf-8 -*-

# network.py
# Jose Vicente Perez Pena
# Dpto. Geodinamica-Universidad de Granada
# 18071 Granada, Spain
# vperez@ugr.es // geolovic@gmail.com
#
# MIT License (see LICENSE file)
# Version: 1.0
# October 1st, 2022
#
# Last modified October 1st, 2022

from . import DEM, Basin, Grid
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np


class HCurve():
    """
    Class to represent and manipulate hypsometric curves

    Parameters:
    -----------
    dem : *landspy.DEM* | *landspy.Basin* | *str*
      DEM, Basin instance or path to a previous saved hypsometric Curve. For a saved path, all remaining arguments are ignored. For a Basin,
      basin and bid are ignored, while name is used.
      If it's a DEM, a Grid with the basin and the bid are needed.

    basin : None, str, Grid
        Drainage basin. If dem is a Basin or a str, this argument is ignored. Needs to have the same dimensions and cellsize than the input DEM.

    bid : int
        ID value that identifies the basin cells.

    name : str, default ""
        Name stored on a newly computed curve; ignored when loading a file.
     """

    def __init__(self, dem, basin=None, bid=1, name=""):

        # If the first parameter is a str (previously saved HCurve), load it
        """Build from a Basin or a DEM and basin mask, or load a saved curve.

        See the class docstring for parameters. Fewer than 50 selected cells
        produce a two-point placeholder with HI=HI2=0.5 and zero moments."""
        if isinstance(dem, str):
            self._load(dem)
            return

        # If is a Basin instance, just take the elevation values
        elif isinstance(dem, Basin):
            values = dem.readArray()
            pos = np.where(dem.getNodata() != values)
            elev_values = values[pos]

        # If not, it is a DEM
        elif isinstance(dem, DEM):
            if isinstance(basin, str):
                basin = Grid(basin)
            elif isinstance(basin, Grid):
                basin = basin
            else:
                raise HypsometryError("Wrong basin grid!!")

            pos = np.where(basin.readArray() == bid)
            elev_values = dem.readArray()[pos]

        if elev_values.size < 50:
            # Not enough values to get a hypsometric curve
            # Just create an empty one
            self._data = np.array([[0, 1], [1, 0]])
            self._name = name
            self._HI = 0.5
            self._HI2 = 0.5
            self.moments = [0, 0, 0, 0, 0]
            return

        elev_values = np.sort(elev_values)[::-1]
        hh = (elev_values - elev_values[-1]) / (elev_values[0] - elev_values[-1])
        aa = (np.arange(elev_values.size) + 1) / (elev_values.size)
        max_h = np.max(elev_values)
        min_h = np.min(elev_values)
        mean_h = np.mean(elev_values)

        # simplify elevation and area data to 10 bits (1024 values)
        if hh.size < 1024:
            self._data = np.array((aa, hh)).T
        else:
            x = np.linspace(0, 1, 1024)
            y = np.interp(x, aa, hh)
            self._data = np.array((x, y)).T

        # get hi values and name
        self._HI = (mean_h - min_h) / (max_h - min_h)
        self._HI2 = self._calculate_hi2()
        self._name = name

        # get statistical moments
        self.moments = self._get_moments()

    def _get_moments(self):
        """
        Fit a cubic to the curve and return five derived moments.
        Order: fitted integral, curve skewness, curve kurtosis, density
        skewness, density kurtosis.
        This code has been translated from a FORTRAN code developed by Harlin
        Harlin, J. M. (1978). Statistical Moments of the Hypsometric Curve and Its Density
        Function. Mathematical Geology, Vol. 10 (1), 59-72.
        """
        coeficients = np.polyfit(self._data[:, 0], self._data[:, 1], 3)[::-1]
        moments = []
        ia = 0
        m10 = 0
        tmp1 = 0
        tmp2 = 0
        tmp3 = 0
        for idx, coef in enumerate(coeficients):
            ia += coef / (idx + 1)
            m10 += coef / (idx + 2)
            tmp1 += coef / (idx + 3)
            tmp2 += coef / (idx + 4)
            tmp3 += coef / (idx + 5)

        tmp1 = tmp1 / ia
        tmp2 = tmp2 / ia
        tmp3 = tmp3 / ia
        moments.append(ia)

        m10 = m10 / ia
        m20 = tmp1 - (m10 ** 2)
        m30 = tmp2 - (3 * m10 * tmp1) + (2 * m10 ** 3)
        m40 = tmp3 - (4 * m10 * tmp2) + (6 * (m10 ** 2) * tmp1) - (3 * m10 ** 4)
        moments.append(m30 / (m20 ** 1.5))
        moments.append(m40 / (m20 ** 2))

        demon = 0
        for coef in coeficients[1:]:
            demon += coef

        tmp1 = 0
        tmp2 = 0
        ex = 0
        ex4 = 0
        for n in range(3):
            ex += (n + 1) * coeficients[n + 1] / (n + 2)
            ex4 += (n + 1) * coeficients[n + 1] / (n + 5)
            tmp1 = tmp1 + (n + 1) * coeficients[n + 1] / (n + 3)
            tmp2 = tmp2 + (n + 1) * coeficients[n + 1] / (n + 4)

        ex = ex / demon
        tmp1 = tmp1 / demon
        tmp2 = tmp2 / demon
        ex2 = tmp1 - ex ** 2
        ex3 = tmp2 - 3.0 * ex * tmp1 + 2.0 * ex ** 3
        tk = ex2 ** 0.5
        ts = tk ** 3
        sdk = ex3 / ts
        ex4 = ex4 / demon - 4.0 * ex * tmp2 + 6.0 * ex ** 2 * tmp1 - 3.0 * ex ** 4
        krd = ex4 / ex2 ** 2

        moments.append(sdk)
        moments.append(krd)

        return moments

    def _calculate_hi2(self):
        """
        Integrate consecutive curve samples using trapezoidal areas.

        Current limitation: the loop omits the last sample interval. The
        historical result is retained pending a numerical-behaviour decision.
        """
        area_accum = 0
        for n in range(len(self._data[:, 0]) - 2):
            a1 = self._data[n, 0]
            a2 = self._data[n + 1, 0]
            h1 = self._data[n, 1]
            h2 = self._data[n + 1, 1]

            if h1 < h2:
                aux = h1
                h1 = h2
                h2 = aux

            area_accum += (a2 - a1) * h2 + ((a2 - a1) * (h1 - h2)) / 2
        return area_accum

    def getA(self):
        """
        Return relative area values of the hypsometric curve
        """
        return self._data[:, 0]

    def getH(self):
        """
        Return relative elevation values of the hypsometric curve
        """
        return self._data[:, 1]

    def getName(self):
        """
        Return the name of the hypsometric curve
        """
        return self._name

    def setName(self, name):
        """
        Sets the name of the hypsometric curve
        """
        self._name = name

    def getHI(self):
        """
        Return the hypsometric integral calculated by the elevation-relief ratio (hmean-hmin) / (hmax-hmin)
        """
        return self._HI

    def getHI2(self):
        """
        Return the hypsometric integral calculated by integrating the curve
        """
        return self._HI2

    def getKurtosis(self):
        """Return moments[1].

        Legacy naming mismatch: _get_moments stores curve skewness at this index,
        not kurtosis. The historical return value is retained."""
        return self.moments[1]

    def getSkewness(self):
        """Return moments[2].

        Legacy naming mismatch: _get_moments stores curve kurtosis at this index,
        not skewness. The historical return value is retained."""
        return self.moments[2]

    def getDensityKurtosis(self):
        """Return moments[3].

        Legacy naming mismatch: this index stores density skewness, not kurtosis.
        The historical return value is retained."""
        return self.moments[3]

    def getDensitySkewness(self):
        """Return moments[4].

        Legacy naming mismatch: this index stores density kurtosis, not skewness.
        The historical return value is retained."""
        return self.moments[4]

    def plot(self, ax=None, **kwargs):
        """
        This function plots the hypsometric curve in an Axes
        :param ax: Axes instance. If None, a new Axes will be created
        :param kwargs: Any matplotlib.Line2D plot property
        :return: None. The plot is added to the supplied or newly created Axes.
        """
        if not ax:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        kw = {"c": "olive", "label": self.getName()}
        kw.update(kwargs)
        ax.plot(self.getA(), self.getH(), **kw)

    def save(self, path):
        """Save curve data and metadata to a semicolon-delimited UTF-8 text file.

        path is used exactly as supplied. Two comment header lines store the name,
        HI, HI2 and the five moments. Returns None."""
        header = self.getName() + "\n"
        moments = [self._HI, self._HI2] + self.moments
        header += ";".join([str(moment) for moment in moments])
        np.savetxt(path, self._data, header=header, comments="#",  delimiter=";", encoding="utf8")

    def _load(self, path):
        # Open the file as normal text file to get its properties
        """Load curve data, name, integrals and moments from a file written by save()."""
        fr = open(path, "r")
        # Line 1: Name
        linea = fr.readline()[1:-1]
        self._name = linea
        # Line 2: Statistical moments
        linea = fr.readline()[1:-1]
        data = linea.split(";")
        self._HI = float(data[0])
        self._HI2 = float(data[1])
        self.moments = [float(dd) for dd in data[2:]]
        fr.close()
        # Load array data
        self._data = np.loadtxt(path, dtype=float, comments='#', delimiter=";", encoding="utf8")


class HypsometryError(Exception):
    """Error raised when a hypsometric calculation receives an invalid basin grid."""
    pass
