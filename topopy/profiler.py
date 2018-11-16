# -*- coding: iso-8859-15 -*-
#
#  profiler.py
#
#  Copyright (C) 2016  J. Vicente Perez, Universidad de Granada
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.

#  For additional information, contact to:
#  Jose Vicente Perez Pena
#  Dpto. Geodinamica-Universidad de Granada
#  18071 Granada, Spain
#  vperez@ugr.es // geolovic@gmail.com

#  Version: 3.0
#  March 02, 2017

#  Last modified 29 November, 2017

import numpy as np
import ogr
import osr
import os


PROFILE_DEFAULT = {'name': "", 'thetaref': 0.45, 'chi0': 0, 'reg_points': 4, 'srs': "", 'smooth': 0}

def canal_2_profile(canal, bnet, bid=0):
    
    di = canal[:, 3]
    mouthdist = di[-1]
    length = di[0] - di[-1]
    di -= di[-1]
    di = length - di
    canal[:, 3] = di
    mcellsize = (bnet._cellsize[0] - bnet._cellsize[1]) / 2 # Mean cellsize
    perfil = TProfile(canal, mcellsize, rid = bid, thetaref = bnet._thetaref, chi0 = canal[-1, 5], 
                  reg_points = bnet._slp_npoints, srs=bnet._proj, mouthdist=mouthdist)

    return perfil


class TProfile:
    """
    Properties:
    ============================
    self.dem_res :: *float*
      Dem resolution of the Digital elevation model used to retreive area-elevation data
    self.rid :: *int*
      Indentifier of the profile
    self.thetaref :: *float*
      Value of m/n used in area-slope calculations
    self.reg_points :: *int*
      Number of points used initially in slope and ksn regressions
    self.n_points :: *int*
      Number of vertexes of the profile
    self._data :: *numpy.array*
      11-column numpy.array with profile data

    =======   ==============================================
    Column    Description
    =======   ==============================================
    c0        X Coordinates of river profile vertex
    c1        Y Coordinates of river profile vertex
    c2        Z Elevation of the river profile vertex
    c3        L Distance from river head
    c4        A Drainage area to the vertex
    c5        Chi (Integral mode)
    c6        Slope of river profile in each vertex
    c7        ksn values for each vertex (calculated by linear regression)
    c8        Quality slope, correlation coefficient (r^2) of slope regressions
    c9        Quality ksn, correlation coefficient (r^2) of ksn regressions
    c10       Raw Z Elevation of the river profile vertex (used to reset the profile)
    =======   ==============================================
    """

    def __init__(self, chandata, dem_res=0, rid=0, thetaref=0.45, chi0=0, reg_points=4, srs="", name="", mouthdist=0, smooth=0):
        """
        Class that defines a river profile with morphometry capabilities.

        :param chandata: *numpy array* - Array with channel data (10 columns array)
        :param dem_res: *float* - Mean resolution of the DEM used to extract profiles
        :param rid: *int* - Profile Identifier
        :param thetaref: *float* - Thetaref (m/n) value used to calculate Chi and Ksn indexes
        :param chi0: *float* - Value of chi index for first point (for tributaries)
        :param name: *str* - Profile name. It will used as the profile label
        :param reg_points: *int* - Number of points (at each side) to calculate initial slope and ksn for each vertex
        :param srs: *str* - Spatial Reference system expresed as well knwon text (wkt)
        :param mouthdist: *float* - Distance from profile to the river mouth (for tributaries)
        :param smooth: *float* - Initial distance to smooth elevations (before to calculate slopes and ksn)

        pf_data param: (numpy.array) with at least 5 columns:

        =======   ==============================================
        Column    Description
        =======   ==============================================
        c0        X Coordinates of river profile vertex
        c1        Y Coordinates of river profile vertex
        c2        Z Elevation of the river profile vertex
        c3        L Distance to head (or to the first vertex)
        c4        A Drainage area of the vertex (in square meters!)
        =======   ==============================================
        """

        # Set profile properties
        self._srs = srs  # WKT with the Spatial Reference
        self._mouthdist = float(mouthdist)
        self._chi0 = chi0
        self._smooth_win = smooth
        self.dem_res = float(dem_res)
        self.rid = int(rid)
        if name == "":
            self.name = str(rid)
        else:
            self.name = str(name)
        self.thetaref = abs(thetaref)
        self.slope_reg_points = reg_points
        self.ksn_reg_points = reg_points
        self.n_points = chandata.shape[0]

       # Get channel data
        self._data = np.copy(chandata)
        # Set raw elevations as 11th column
        self._data = np.append(self._data, self._data[:, 2].reshape(self.n_points, 1), axis=1)
        
        # List for knickpoints and regressions
        self._knickpoints = [] # Tuples of (pos, type)
        self._regressions = [] # Tuples of (pos1, pos2, ksn)
        self._smoothpoints = 0


    def get_projection(self):
        """
        Returns a string with the projection as wkt
        """
        return self._srs

    def set_projection(self, projection):
        """
        Set the projection of the profile

        :param projection: str - String with the projection in wkt
        :return: None
        """
        self._srs = projection

    def length(self):
        """
        Returns the total length of the profile
        """
        return self._data[-1, 3]

    def get_x(self, head=True):
        """
        Returns x coordinates for all vertices

        :param head: boolean - Specifies if x coordinates are returned from head (True) or mouth (False)
        :return: numpy.array wiht x values for all vertices
        """
        if head:
            return np.copy(self._data[:, 0])
        else:
            return np.copy(self._data[::-1, 0])

    def get_y(self, head=True):
        """
        Returns y coordinates for all vertices

        :param head: boolean - Specifies if y coordinates are returned from head (True) or mouth (False)
        :return: numpy.array wiht y values for all vertices
        """
        if head:
            return np.copy(self._data[:, 1])
        else:
            return np.copy(self._data[::-1, 1])

    def get_z(self, head=True, relative=False):
        """
        Returns elevation values for all vertices

        :param head: boolean - Specifies if elevations are returned from head (True) or mouth (False)
        :param relative: boolean - Specifies if elevations are relative (min elevation = 0) or not
        :return: numpy.array wiht elevation values for all vertices
        """
        z_values = np.copy(self._data[:, 2])
        if relative:
            z_min = z_values[-1]
            z_values -= z_min

        if head:
            return z_values
        else:
            return z_values[::-1]

    def get_raw_z(self, head=True, relative=False):
        """
        Returns raw elevation values for all vertices

        :param head: boolean - Specifies if raw elevations are returned from head (True) or mouth (False)
        :param relative: boolean - Specifies if elevations are relative (min elevation = 0) or not
        :return: numpy.array wiht elevation values for all vertices
        """
        raw_z = np.copy(self._data[:, 10])
        if relative:
            z_min = raw_z[-1]
            raw_z -= z_min

        if head:
            return raw_z
        else:
            return raw_z[::-1]

    def get_l(self, head=True):
        """
        Returns a numpy.array with distances for all profile vertices

        :param head: boolean - Specifies if distances are returned from head (True) or mouth (False)
        If measured from mouth, a initial distance (self._mouthdist) will be added (to account for tributaries)
        :return: numpy.array with distances for all vertices (measured from head or mouth)
        """
        river_length = float(self._data[-1, 3])

        if head:
            li = np.copy(self._data[:, 3])
        else:
            li = river_length - self._data[:, 3] + self._mouthdist
            li = li[::-1]

        return li

    def get_area(self, head=True):
        """
        Returns a numpy.array with drainage area values for all vertices

        :param head: boolean - Specifies if areas are returned from head (True) or mouth (False)
        :return: numpy.array wiht area values for all vertices
        """
        areas = np.copy(self._data[:, 4])

        if head:
            return areas
        else:
            return areas[::-1]

    def get_slope(self, threshold=0, head=True, lq=False):
        """
        Returns slopes calculated by linear regression

        :param threshold: *float* R^2 threshold. (Slopes with R^2 < threshold will be in lq_slopes array)
        :param head: boolean - Specifies if slopes are returned from head (True) or mouth (False)
        :param lq: boolean - Specifies lq_slopes will be returned or not
        :return: array with slopes (lq_slopes=False) or tuple of arrays (slopes, lq_slopes) (lq_slopes=True)
         slopes --> numpy.array of slopes with R^2 >= threshold (lq_slopes will receive a np.nan value)
         lq_slopes --> numpy.array of slopes with R^2 < threshold (slopes will receive a np.nan value)
        """
        slopes = []
        lq_slopes = []
        r2_values = self.get_slope_r2()
        for n in range(len(self._data)):
            if r2_values[n] >= threshold:
                slopes.append(self._data[n, 6])
                lq_slopes.append(np.nan)
            else:
                slopes.append(np.nan)
                lq_slopes.append(self._data[n, 6])

        slopes = np.array(slopes)
        lq_slopes = np.array(lq_slopes)

        if not head:
            slopes = slopes[::-1]
            lq_slopes = lq_slopes[::-1]

        if lq:
            return slopes, lq_slopes
        else:
            return slopes

    def get_ksn(self, threshold=0, head=True, lq=False):
        """
        Returns ksn values calculated by linear regression

        :param threshold: *float* R^2 threshold. (ksn with R^2 < threshold will be in lq_ksn array)
        :param head: boolean - Specifies if ksn are returned from head (True) or mouth (False)
        :param lq: boolean - Specifies lq_ksn will be returned or not
        :return: array with ksn (lq=False) or tuple of arrays (ksn, lq_ksn) (lq=True)
         ksn --> numpy.array of ksn values with R^2 >= threshold (lq_ksn will receive a np.nan value)
         lq_ksn --> numpy.array of ksn values with R^2 < threshold (ksn will receive a np.nan value)
        """
        ksn = []
        lq_ksn = []
        ksn_r2 = self.get_ksn_r2()
        for n in range(len(self._data)):
            if ksn_r2[n] >= threshold:
                ksn.append(self._data[n, 7])
                lq_ksn.append(np.nan)
            else:
                ksn.append(np.nan)
                lq_ksn.append(self._data[n, 7])

        ksn = np.array(ksn)
        lq_ksn = np.array(lq_ksn)

        if not head:
            ksn = ksn[::-1]
            lq_ksn = lq_ksn[::-1]

        if lq:
            return ksn, lq_ksn
        else:
            return ksn

    def get_slope_r2(self, head=True):
        """
        Returns slope R2 values from linear regressions for all vertices

        :param head: boolean - Specifies if R2 values are returned from head (True) or mouth (False)
        :return: numpy.array wiht slope R2 values for all vertices
        """
        if head:
            return np.copy(self._data[:, 8])
        else:
            return np.copy(self._data[::-1, 8])

    def get_ksn_r2(self, head=True):
        """
        Returns ksn R2 values from linear regressions for all vertices

        :param head: boolean - Specifies if R2 values are returned from head (True) or mouth (False)
        :return: numpy.array wiht ksn R2 values for all vertices
        """
        if head:
            return np.copy(self._data[:, 9])
        else:
            return np.copy(self._data[::-1, 9])

    def get_chi(self, head=True, relative=False):
        """
        Returns chi values for all vertices in ascending order.

        :param head: boolean - Specifies if chi values are returned from head (True) or mouth (False)
        :param relative: boolean - Specifies if chi values are relative (min chi = 0) or not
        :return: numpy.array wiht chi values for all vertices
        """
        chi_values = np.copy(self._data[:, 5])
        if relative:
            chi0 = chi_values[-1]
            chi_values -= chi0

        if head:
            return chi_values
        else:
            return chi_values[::-1]

    def smooth(self, npoints = 0):
        """
        Smooths the elevations of the profile with a movil mean of window size. It also removes peaks and flat segments
        from the profile (to avoid problems when calculating slopes)

        :param window: Number of points (at each side) to smooth the elevations of the river profile
        :return: None
        """

        # Remove peaks and flat segments
        for n in range(len(self._data) - 1):
            if self._data[n + 1, 2] >= self._data[n, 2]:
                self._data[n + 1, 2] = float(self._data[n, 2]) - 0.001

        # Smooth elevations if window distance > 0
        if npoints > 0:
            for ind in range(len(self._data)):
                low = ind - npoints
                high = ind + npoints + 1
                if low < 0:
                    low = 0

                elevations = self._data[low:high, 10]
                self._data[ind, 2] = np.mean(elevations)

    def reset_elevations(self):
        """
        Reset smooth elevations. When reset, smooth elevations will equal to raw elevations
        """
        self._smoothpoints = 0
        for n in range(len(self._data)):
            self._data[n, 2] = np.copy(self._data[n, 10])
        
    def calculate_chi(self, a0=1, chi0=0.0):
        """
        This function creates the chi data array. Chi data will be calculated for each vertex of the river profile and
        stored in the 7th column of self._data numpy array

        :param a0: *int* - Reference area to remove dimensionality of Chi index
        :param chi0: *float* - Initial Chi value in profile mouth. Needed to calculate chi values for tributaries
        """
        # Invert area array
        ai = self._data[::-1, 4] ** self.thetaref
        a0 = a0 ** self.thetaref
        chi = [chi0]
        for n in range(len(ai)):
            if n > 0:
                dx = self._data[n, 3] - self._data[n - 1, 3]
                chi.append(chi[n - 1] + (a0 * dx / ai[n]))

        self._data[:, 5] = chi[::-1]

    def calculate_slope(self, reg_points=4, raw_z=False):
        """
        This function calculates slopes for all vertexes by linear regression of distance-elevation data.
        Slopes are stored in column c5 of self._data. Together with slopes, R^2 are calculated (column  c6)

        :param reg_points: Number of profile points before and after each vertex to calculate slope
        :param raw_z: bool Specifies if raw_z values are taken from calculate slopes (True) or not (False)
        :return: None
        """
        self.slope_reg_points = reg_points
        li = self.get_l()
        if raw_z:
            zi = self.get_raw_z()
        else:
            zi = self.get_z()

        for n in range(self.n_points):
            low = n - reg_points
            high = n + reg_points

            if low < 0:
                low = 0

            xi = li[low:high + 1]
            yi = zi[low:high + 1]
            poli, SCR = np.polyfit(xi, yi, deg = 1, full = True)[:2]
            
            if yi.size * yi.var() == 0:
                R2 = 0.0
            else:
                R2 = float(1 - SCR/(yi.size * yi.var()))

            gradient = poli[0]
            self._data[n, 8] = abs(R2)

            if abs(gradient) < 0.001:
                self._data[n, 6] = 0.001
            else:
                self._data[n, 6] = abs(gradient)

    def calculate_ksn(self, reg_points=4, raw_z=False):
        """
        This function calculates ksn for all vertexes by linear regression of chi-elevation data.

        :param reg_points: *int* - Number of profile points before and after each vertex to calculate ksn
        :param raw_z: bool Specifies if raw_z values are taken from calculate slopes (True) or not (False)
        :return: numpy.array with ksn values for all vertexes. If full is true, it returns a tuple of arrays (ksn, r^2)
        """
        self.ksn_reg_points = reg_points
        ksn_values = []
        ksn_r2_values = []
        chi = self.get_chi(False)
        if raw_z:
            zi = self.get_raw_z(False)
        else:
            zi = self.get_z(False)

        for n in range(self.n_points):
            low = n - reg_points
            high = n + reg_points

            if low < 0:
                low = 0

            xi = chi[low:high + 1]
            yi = zi[low:high + 1]

            poli, SCR = np.polyfit(xi, yi, deg = 1, full = True)[:2]

            if yi.size * yi.var() == 0:
                R2 = 0.0
            else:
                R2 = float(1 - SCR/(yi.size * yi.var()))

            gradient = poli[0]

            ksn_r2_values.append(abs(R2))

            if abs(gradient) < 0.0001:
                ksn_values.append(0.0001)
            else:
                ksn_values.append(abs(gradient))

        self._data[:, 7] = np.array(ksn_values)[::-1]
        self._data[:, 9] = np.array(ksn_r2_values)[::-1]

    def get_best_theta(self, a0=1, step=0.05):
        """
        Description
        ===========
        This function obtain the best m/n value for the profile following the approach
        proposed in Perron and Royden, 2013. This best m/n value will be the one that
        increases the linearity of the Chi-Elevation profile

        Parameters:
        ==============
        a0 :: *int (Default = 1)*
          Reference area value. By default set as 1 square meter

        step :: *float (Default = 0.05)*
          Step to test the different theta values. Recommended 0.1 or 0.05

        Returns:
        ==============
        best_theta :: *float*
          Best m/n value for the profile. The one that icreases the linearity
          for the whole Chi-Zi profile
        """
        best_r2 = 0
        best_theta = 0
        theta_values = np.arange(0, 1, step)
        zi = self._data[::-1, 2]

        for theta in theta_values:
            ai = self._data[::-1, 4] ** theta
            a0 = a0 ** theta
            chi = [0]
            for n in range(len(ai)):
                if n > 0:
                    dx = self._data[n, 3] - self._data[n - 1, 3]
                    chi.append(chi[n - 1] + (a0 * dx / ai[n]))

            # Regresion into chi-elevation space to get r^2
            a1 = np.array([chi, np.ones(len(chi))]).T
            y1 = zi
            model, resid = np.linalg.lstsq(a1, y1)[:2]
            r2 = 1 - resid / (y1.size * y1.var())
            if r2 > best_r2 and best_theta != 0:
                best_r2 = r2
                best_theta = theta

        return best_theta



def heads_from_points(point_shp, id_field=""):
    """
    Extract heads from point shapefile to a 3-column numpy ndarray [x, y, id]
    
    Parameters:
    ===========
    point_shp : *str* 
      Path to the point shapefile
    id_field : *str* 
      Id field (heads will be ordered by these field values)
    
    Returns:
    ========
    **numpy.ndarray** with 3 columns [X, Y, id]
    """

    heads = []
    dataset = ogr.Open(point_shp)
    layer = dataset.GetLayer(0)
    layerdef = layer.GetLayerDefn()
    idx = layerdef.GetFieldIndex(id_field)
    take_field = False
    if idx >= 0:
        field_defn = layerdef.GetFieldDefn(idx)
        if field_defn.GetType() in [0, 12]:
            take_field = True

    n = 0.
    for feat in layer:
        geom = feat.GetGeometryRef()
        if take_field:
            hid = feat.GetField(id_field)
        else:
            hid = n + 1
        n += 1
        heads.append((geom.GetX(), geom.GetY(), hid))
    
    return np.array(heads)


def profiles_to_shp(path, profiles, distance=0):
    """
    This function saves a list of profiles in a shapefile (point or line shapefile, depending of the distance parameter)

    :param path: str - Full path to the shapefile
    :param profiles: list - List of TProfile objects
    :param distance: float - Distance for profile segments. If distance = 0, the output shapefile will be a point shp
    :return: None
    """
    if distance == 0:
        _profiles_to_points(path, profiles)
    else:
        _profiles_to_lines(path, profiles, distance)


def _profiles_to_points(out_shp, profiles):
    """
    This function save a list of profiles in a point shapefile

    :param out_shp: str - Full path to the shapefile
    :param profiles: list - List of TProfile objects
    :return: None
    """
    # Create point shapefle
    driver = ogr.GetDriverByName("ESRI Shapefile")
    if os.path.exists(out_shp):
        dataset = driver.Open(out_shp, 1)
        for n in range(dataset.GetLayerCount()):
            dataset.DeleteLayer(n)
    else:
        dataset = driver.CreateDataSource(out_shp)
    sp = osr.SpatialReference()
    sp.ImportFromWkt(profiles[0].get_projection())
    layer = dataset.CreateLayer("perfiles", sp, ogr.wkbPoint)

    # Add fields
    campos = ["id_profile", "L", "area", "z", "chi", "ksn", "rksn", "slope", "rslope"]
    tipos = [0, 2, 12, 2, 2, 2, 2, 2, 2]
    for n in range(len(campos)):
        layer.CreateField(ogr.FieldDefn(campos[n], tipos[n]))

    # Get data from all the profiles
    for profile in profiles:
        id_perfil = profile.rid
        xi = profile.get_x()
        yi = profile.get_y()
        li = profile.get_l()
        ai = profile.get_area()
        zi = profile.get_z()
        chi = profile.get_chi()
        ksn = profile.get_ksn()
        rksn = profile.get_ksn_r2()
        slp = profile.get_slope()
        rslp = profile.get_slope_r2()
        id_arr = np.zeros(xi.size, dtype="int32")
        id_arr.fill(id_perfil)
        data = np.array((xi, yi, id_arr, li, ai, zi, chi, ksn, rksn, slp, rslp)).T
        data = data[:-1]  # Remove the last point, because it will have the area of the merging channel

        for row in data:
            row_list = row.tolist()  # To avoid numpy types
            feat = ogr.Feature(layer.GetLayerDefn())
            geom = ogr.Geometry(ogr.wkbPoint)
            geom.AddPoint(row_list[0], row_list[1])
            feat.SetGeometry(geom)

            for idx, element in enumerate(row_list[2:]):
                feat.SetField(campos[idx], element)

            layer.CreateFeature(feat)

        id_perfil += 1


def _profiles_to_lines(out_shp, profiles, distance):
    """
    This function save a list of profiles in a line shapefile

    :param out_shp: str - Full path to the shapefile
    :param profiles: list - List of TProfile objects
    :param distance: float - Distance for the profile segments
    :return: None
    """

    # Create line shapefle
    driver = ogr.GetDriverByName("ESRI Shapefile")
    if os.path.exists(out_shp):
        dataset = driver.Open(out_shp, 1)
        for n in range(dataset.GetLayerCount()):
            dataset.DeleteLayer(n)
    else:
        dataset = driver.CreateDataSource(out_shp)
    sp = osr.SpatialReference()
    sp.ImportFromWkt(profiles[0].get_projection())
    layer = dataset.CreateLayer("perfiles", sp, ogr.wkbLineString)

    # Add fields
    campos = ["id_profile", "L", "area", "z", "chi", "ksn", "rksn", "slope", "rslope"]
    tipos = [0, 2, 12, 2, 2, 2, 2, 2, 2]
    for n in range(len(campos)):
        layer.CreateField(ogr.FieldDefn(campos[n], tipos[n]))

    for perfil in profiles:
        id_perfil = perfil.rid
        # Transformamos la distancia a numero de celdas
        # dx es el numero de puntos que tendra cada segmento (numero impar)
        cellsize = perfil.dem_res
        if distance == 0:
            dx = len(perfil._data)
        else:
            dx = int((distance / cellsize) + 0.5)
            if dx % 2 == 0:
                dx += 1

        p1 = 0
        p2 = p1 + dx

        # Get profile data from mouth to head
        xi = perfil.get_x(False)
        yi = perfil.get_y(False)
        chi = perfil.get_chi(False)
        zi = perfil.get_z(False)
        li = perfil.get_l()[::-1]
        area = perfil.get_area(False)

        while p1 < len(chi) - 3:
            # Get the data for the segments
            if p2 >= len(chi):
                p2 = len(chi)
            sample_x = xi[p1:p2]
            sample_y = yi[p1:p2]
            sample_chi = chi[p1:p2]
            sample_zi = zi[p1:p2]
            sample_li = li[p1:p2]

            # Ksn by regression
            coef, resid = np.polyfit(sample_chi, sample_zi, deg=1, full=True)[:2]
            ksn = float(abs(coef[0]))
            if (sample_zi.size * sample_zi.var()) == 0:
                rksn = 0.
            else:
                rksn = float(1 - resid / (sample_zi.size * sample_zi.var()))

            # Slope by regression
            coef, resid = np.polyfit(sample_li, sample_zi, deg=1, full=True)[:2]
            slope = float(abs(coef[0]))
            if (sample_zi.size * sample_zi.var()) == 0:
                rslope = 0.
            else:
                rslope = float(1 - resid / (sample_zi.size * sample_zi.var()))

            # Get the other values
            mid_pos = p1 + int((p2 - p1) / 2)
            mid_chi = float(chi[mid_pos])
            mid_area = int(area[mid_pos])
            mid_l = float(li[mid_pos])
            mid_z = float(zi[mid_pos])

            # Get geometry for the line
            geom = ogr.Geometry(ogr.wkbLineString)
            for n in range(len(sample_x)):
                geom.AddPoint(float(sample_x[n]), float(sample_y[n]))

            # Create the feature
            feature = ogr.Feature(layer.GetLayerDefn())

            # Add values to feature
            values = [id_perfil, mid_l, mid_area, mid_z, mid_chi, ksn, rksn, slope, rslope]
            for idx, value in enumerate(values):
                feature.SetField(campos[idx], int(value))

            # Add geometry to feature and feature to layer
            feature.SetGeometry(geom)
            layer.CreateFeature(feature)

            # Last point of the line was p2-1
            p1 = p2 - 1
            p2 = p1 + dx


def version():
    return "Version: 26 October 2017"
