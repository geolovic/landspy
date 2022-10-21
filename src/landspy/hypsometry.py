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
      DEM, Basin instance or path to a previous saved hypsometric Curve. If it's a Basin or a string, the rest of the arguments will be ignored.
      If it's a DEM, a Grid with the basin and the basin_id are needed.

    basin : None, str, Grid
        Drainage basin. If dem is a Basin or a str, this argument is ignored. Needs to have the same dimensions and cellsize than the input DEM.

    bid : int
        Id value that identifies the basin cells
     """

    def __init__(self, dem, basin=None, bid=1, name=""):

        # If the first parameter is a str (previously saved HCurve), load it
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
        Get statistical moments given the 4 coefficients of the curve polynomial regression
        This code has been translated from a FORTRAN code developed by Harlin (1979)
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
        Calculates HI by integrating the hypsometric curve
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
        Return the hypsometric integral calculated by aproximation (hmed-hmin) / (hmax-hmin)
        """
        return self._HI

    def getHI2(self):
        """
        Return the hipsometric integral calculated by integrating the curve
        """
        return self._HI2

    def getKurtosis(self):
        return self.moments[1]

    def getSkewness(self):
        return self.moments[2]

    def getDensityKurtosis(self):
        return self.moments[3]

    def getDensitySkewness(self):
        return self.moments[4]

    def plot(self, ax=None, **kwargs):
        """
        This function plots the hypsometric curve in an Axe
        :param ax: Axes instance. If None, a new Axe will be created
        :param kwargs: Any matplotlib.Line2D plot property
        :return:
        """
        if not ax:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        kw = {"c":"olive", "label":self.getName()}
        kw.update(kwargs)
        ax.plot(self.getA(), self.getH(), **kw)

    def save(self, path):
        header = self.getName() + "\n"
        moments = [self._HI, self._HI2] + self.moments
        header += ";".join([str(moment) for moment in moments])
        np.savetxt(path, self._data, header=header, comments="#",  delimiter=";", encoding="utf8")

    def _load(self, path):
        # Open the file as normal text file to get its properties
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
    pass
