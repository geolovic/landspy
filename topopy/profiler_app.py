# -*- coding: iso-8859-15 -*-
#
#  ProfilerApp1.py
#
#  Copyright (C) 2017  J. Vicente Perez, Universidad de Granada
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

#  Version: 2
#  October 08, 2017

#  Last modified 29 November, 2017

import matplotlib.pyplot as plt
import numpy as np
import ogr
import osr
import os
#from . import TProfile
from profiler import TProfile


# Disable some keymap characters that interfere with graph key events
plt.rcParams["keymap.xscale"] = [""]
plt.rcParams["keymap.yscale"] = [""]
plt.rcParams["keymap.save"] = [u'ctrl+s']
key_dict = {"LEFT": -1, "RIGHT": 1, "+": 1, "-": -1}
k_type = {0: "r*", 1: "g*", 2: "r.", 3: "y."}  # Knickpoint types


class ProfilerApp:
    """
    Class to draw longitudinal, area-slope, chi and ksn profiles.
    Allows selecting knickpoints and do regressions on chi-elevation profiles
    """
    gr_types = {1: 'longitudinal', 2: "area-slope", 3: "chi", 4: "ksn"}
    
    def __init__(self, figure, profiles, basedir=""):
        """
        :param figure: Matplotlib.Figure to draw graphics
        :param profiles: numpy.array with TProfiles objects
        :param basedir: Base directory to save regressions and knickpoints
        """
        # Create the output folder (for knickpoints and regressions)
        if not basedir:
            self.basedir = "out_files"
        elif basedir[-1] == "/":
            self.basedir = basedir + "out_files"
        else:
            self.basedir = basedir + "/out_files"
        
        if not os.path.exists(self.basedir):
            os.mkdir(self.basedir)
        
        self.figure = figure
        self.ax = figure.add_subplot(111)
        self.g_type = 1
        self.profiles = profiles
        self.active = 0
        self.mode = None
        self.regression = None
        self.mode_txt = self.figure.text(0.025, 0.95, "", color="navy")
        self.n_profiles = len(profiles)
        self.dam_points = []  # Temporary list of 2 elements [ind0, ind1]
        self.activate()
        self.draw()
   
    def activate(self):
        self.kcid = self.figure.canvas.mpl_connect("key_press_event", self.key_input)      

    # DRAWING METHODS
    # ===============
    def draw(self):
        """
        Draw the selected graphic for the active profile
        """
        if self.g_type == 1:
            self._draw_long_profile()
        elif self.g_type == 2:
            self._draw_area_slope()
        elif self.g_type == 3:
            self._draw_chi_profile()
        elif self.g_type == 4:
            self._draw_ksn_profile()
        self._mode_inf()
        self.figure.canvas.draw()

    def _draw_long_profile(self):
        """
        Draw the longitudinal profile for the active profile
        """
        self.ax.clear()
        perfil = self.profiles[self.active]
        li = perfil.get_l() / 1000.
        zi = perfil.get_z()
        self.ax.set_xlabel("Distance [km]")
        self.ax.set_ylabel("Elevation [m]")
        self.ax.set_title(perfil.name)
        self.ax.plot(li, zi, c="b", lw=0.7, picker=4)
        
        # Draw knickpoints
        if len(perfil._knickpoints) > 0:
            for k in perfil._knickpoints:
                self.ax.plot(li[k[0]], zi[k[0]], k_type[k[1]], mew=0.5, mec="k", ms=10)
        
    def _draw_area_slope(self):
        """
        Draw area-slope log-log plot for the active profile
        """
        self.ax.clear()
        perfil = self.profiles[self.active]
        slopes = perfil.get_slope()
        areas = perfil.get_area()
        self.ax.set_xlabel("Area $m^2$")
        self.ax.set_ylabel("Slope (reg. points: {0})".format(perfil.slope_reg_points))
        self.ax.set_xscale("log")
        self.ax.set_yscale("log")
        self.ax.set_title(perfil.name)
        self.ax.plot(areas, slopes, "b+", mew=0.5, picker=4)

        # Draw knickpoints
        if len(perfil._knickpoints) > 0:
            for k in perfil._knickpoints:
                self.ax.plot(areas[k[0]], slopes[k[0]], k_type[k[1]], mew=0.5, mec="k", ms=10)

    def _draw_chi_profile(self):
        """
        Draw chi-elevation plot for the active profile
        """
        self.ax.clear()
        perfil = self.profiles[self.active]
        chi = perfil.get_chi()
        zi = perfil.get_z()
        self.ax.set_xlabel("Chi [m]")
        self.ax.set_ylabel("Elevation [m]")
        self.ax.set_title(perfil.name)
        self.ax.plot(chi, zi, c="b", lw=0.7, picker=4)

#        # Draw regressions
#        for r in self.regressions[self.active]:
#            self.ax.plot(r[4][:, 0], r[4][:, 1], c="r", ls="--", lw=1)
#            p1x = self.profiles[self.active].get_chi()[r[0]]
#            p1y = self.profiles[self.active].get_z()[r[0]]
#            p2x = self.profiles[self.active].get_chi()[r[1]]
#            p2y = self.profiles[self.active].get_z()[r[1]]
#            self.ax.plot([p1x, p2x], [p1y, p2y], "k+")

        # Draw knickpoints
        if len(perfil._knickpoints) > 0:
            for k in perfil._knickpoints:
                self.ax.plot(chi[k[0]], zi[k[0]], k_type[k[1]], mew=0.5, mec="k", ms=10)

    def _draw_ksn_profile(self):
        """
        Draw ksn-distance plot for the active profile
        """
        self.ax.clear()
        perfil = self.profiles[self.active]
        ksn = perfil.get_ksn()
        li = perfil.get_l() / 1000.
        self.ax.set_xlabel("Distance [km]")
        self.ax.set_ylabel("Ksn (reg. points: {0})".format(perfil.ksn_reg_points))
        self.ax.set_title(perfil.name)
        self.ax.plot(li, ksn, c="b", lw=0.7, picker=4)

        # Draw knickpoints
        if len(perfil._knickpoints) > 0:
            for k in perfil._knickpoints:
                self.ax.plot(li[k[0]], ksn[k[0]], k_type[k[1]], mew=0.5, mec="k", ms=10)

    def _mode_inf(self):
        """
        Set the text to inform about the drawing mode 
        R - Regresion mode
        K - Knickpoint mode
        D - Dam remover mode
        """
        if self.mode:
            self.mode_txt.set_text("Mode: {0}".format(self.mode))
        else:
            self.mode_txt.set_text("")
                       
    # EVENT METHODS
    # =============
    def key_input(self, event):
        """
        Controls keyboard events to add functionality to application
        :param event: keypress event
        
        Keyboard actions:
            1-4: Changes the graphic type [long, area/slope, chi, ksn]
            R: Enters in regression mode
            K: Enters in knickpoint mode
            D: Enters in dam remover mode
            C: Clears the last regression or last knickpoint (if in mode R or K)
            + / - : Increases/reduces reg points in area-slope or ksn plots
            S: Save all the knickpoints and regressions
            M: Draws a map with all the profiles and knickpoints (hightlighting the active one)
        """
        key = event.key.upper()
        if key in ["1", "2", "3", "4"]:
            self.g_type = int(key)
            self.draw()
                        
        elif key in ["LEFT", "RIGHT"]:
            self.next_profile(key_dict[key])

        elif key in ["K", "R", "D"]:
            self.change_mode(key)
        
        elif key == "UP":
            self.change_knickpoint_type()
        
        elif key == "C":
            self.remove_last()
            
        elif key in ["+", "-"] and self.g_type in [2, 4]:
            self.change_reg_points(key_dict[key])
        
        elif key == "S":
            self.save()
            
        elif key == "M":
            self.draw_map()
            
        # elif key == "P":
        #     self.print_all(4, 3)

#    def change_knickpoint_type(self):
#        if self.mode == "K" and len(self.knick_points[self.active]) > 0:
#            kp = self.knick_points[self.active][-1]
#            kp[1] = (kp[1] + 1) % len(k_type)
#            #self.knick_points[self.active][-1][1] = (self.knick_points[self.active][-1][1] + 1) % 2
#            self.draw()    
        
    def pick_point(self, event):
        """
        Pick a point in the active profile 
        :param event: matplotlib picker event
        """
        # In the case that many points are picked (true if the profile has several points). Take the center one.
        if len(event.ind) > 2:
            ind = (event.ind[-1] + event.ind[0]) // 2
        else:
            ind = event.ind[0]

        # Select the proper option according with the drawing mode
        if self.mode == "K":
            self.profiles[self.active]._knickpoints.append([ind, 0])
            self.draw()
            
#        elif self.mode == "R" and self.g_type == 3:
#            self.reg_points.append(ind)
#            chi = self.profiles[self.active].get_chi()[ind]
#            zi = self.profiles[self.active].get_z()[ind]
#            self.ax.plot(chi, zi, "k+")
#            self.figure.canvas.draw()
#            self.create_regression()   
#            
#        elif self.mode == "D" and self.g_type == 3:
#            self.dam_points.append(ind)
#            chi = self.profiles[self.active].get_chi()
#            zi = self.profiles[self.active].get_z()
#            c_point = chi[ind]
#            z_point = zi[ind]
#            self.ax.plot(c_point, z_point, "k+", markersize=7)
#            self.figure.canvas.draw()
#            self.remove_dam()
#
#        elif self.mode == "D" and self.g_type == 1:
#            self.dam_points.append(ind)
#            li = self.profiles[self.active].get_l() / 1000.
#            zi = self.profiles[self.active].get_z()
#            c_point = li[ind]
#            z_point = zi[ind]
#            self.ax.plot(c_point, z_point, "k+", markersize=7)
#            self.figure.canvas.draw()
#            self.remove_dam2()

    def remove_dam(self):
        """
        Removes points that correspond with a dam in the longitudinal profile
        It will take points from self.dam_points[] >> list with two points (positions)
        """
        if len(self.dam_points) < 2:
            return
        elif len(self.dam_points) > 2:
            self.dam_points = []
            self.draw()

        p1 = self.dam_points[0]
        p2 = self.dam_points[1]
        if p1 > p2:
            tmp = p2
            p2 = p1
            p1 = tmp
            
        perfil = self.profiles[self.active]
        
        # Create a straigt-line between both points
        chi = perfil.get_chi()
        zi = perfil.get_z()
        arr_inds = list(range(p1, p2))
        
        pto1 = (chi[p1], zi[p1])
        pto2 = (chi[p2], zi[p2])        

        m = (pto2[1] - pto1[1]) / float(pto2[0] - pto1[0])
        b = pto1[1] - (m * pto1[0])
        
        xi = chi[arr_inds]
        yi = xi * m + b      

        perfil._data[p1:p2, 2] = yi
        perfil.calculate_slope(perfil.slope_reg_points)
        perfil.calculate_ksn(perfil.ksn_reg_points)
        self.dam_points = []
        self.draw()

    def remove_dam2(self):
        """
        Removes points that correspond with a dam in the longitudinal profile
        It will take points from self.dam_points[] >> list with two points (positions)
        """
        if len(self.dam_points) < 2:
            return
        elif len(self.dam_points) > 2:
            self.dam_points = []
            self.draw()

        p1 = self.dam_points[0]
        p2 = self.dam_points[1]
        if p1 > p2:
            p1, p2 = p2, p1

        perfil = self.profiles[self.active]

        # Create a straigt-line between both points
        li = perfil.get_l() / 1000.
        zi = perfil.get_z()
        arr_inds = list(range(p1, p2))

        pto1 = (li[p1], zi[p1])
        pto2 = (li[p2], zi[p2])

        m = (pto2[1] - pto1[1]) / float(pto2[0] - pto1[0])
        b = pto1[1] - (m * pto1[0])

        xi = li[arr_inds]
        yi = xi * m + b

        perfil._data[p1:p2, 2] = yi
        perfil.calculate_slope(perfil.slope_reg_points)
        perfil.calculate_ksn(perfil.ksn_reg_points)
        self.dam_points = []
        self.draw()

    def create_regression(self):
        """
        Creates a regression in Chi-Elevation space
        It will take points from self.reg_points[] >> list with two points (positions)
        """
        if len(self.reg_points) < 2:
            return
        if len(self.reg_points) > 2:
            self.reg_points = []
            self.draw()
            return

        p1 = self.reg_points[0]
        p2 = self.reg_points[1]
        if p1 > p2:
            tmp = p2
            p2 = p1
            p1 = tmp
            
        perfil = self.profiles[self.active]
        
        # Regression in chi plot
        chi = perfil.get_chi()[p1:p2+1]
        zi = perfil.get_z()[p1:p2+1]
        poli, scr = np.polyfit(chi, zi, deg=1, full=True)[:2]
        r2 = float(1 - scr/(zi.size * zi.var()))
        xc = np.linspace(chi[0], chi[-1], num=5)
        yc = np.polyval(poli, xc)
        arr = np.array((xc, yc)).T
        self.regressions[self.active].append((p1, p2, poli[0], poli[1], arr, r2))
        
        # Once the regression is created, force redraw and clear reg_point list
        self.reg_points = []
        self.draw()

    def change_mode(self, mode):
        """
        Change the interactive mode:
            K >> Knickpoint selector
            R >> Regression mode (only works in Chi profile, self.gtype == 3)
            D >> Damm remover (only works in Longitudinal profile, self.gtype==1)
        """
        if self.mode == mode:
            self.mode = None
            self.figure.canvas.mpl_disconnect(self.pcid)
        else:
            self.mode = mode
            self.pcid = self.figure.canvas.mpl_connect("pick_event", self.pick_point)
        self.draw()

    def remove_last(self):
        """
        Remove the last knickpoint or the regression from the active profile
        """
        if self.mode == "K" and len(self.knick_points[self.active]) > 0:
            self.knick_points[self.active].pop()
            self.draw()
        elif self.mode == "R" and len(self.regressions[self.active]) > 0:
            self.regressions[self.active].pop()
            self.draw()          

    def change_reg_points(self, increment):
        """
        Changes the regressions points in slope-area and ksn profile
        """
        perfil = self.profiles[self.active]
        if self.g_type == 2:
            reg_points = perfil.slope_reg_points + increment
            if reg_points > 2:
                perfil.calculate_slope(reg_points)
        elif self.g_type == 4:
            reg_points = perfil.ksn_reg_points + increment
            if reg_points > 2:
                perfil.calculate_ksn(reg_points)
                           
        self.draw()

    def next_profile(self, idx):
        """
        Get the next profile from profile lists and set it as active profile
        """
        self.active += idx
        self.active %= self.n_profiles
        self.draw()
            
    def draw_map(self):
        """
        Draws a plot with all the analyzed profiles and highlight the active one
        It also draws knickpoints.
        """
        # Create a new Figure and Axes        
        figure = plt.figure()
        ax = figure.add_subplot(111)
        
        # Draw a map with all the profiles
        for idx, perfil in enumerate(self.profiles):
            xi = perfil.get_x()
            yi = perfil.get_y()
            ax.plot(xi, yi, c="0.6", lw=0.5)
            # Draw the nickpoints
            for k in self.knick_points[idx]:
                ax.plot(xi[k[0]], yi[k[0]], k_type[k[1]], mew=0.5, mec="k", ms=7)    
        
        # Draw the active profile with different color            
        perfil = self.profiles[self.active]
        xi = perfil.get_x()
        yi = perfil.get_y()
        ax.plot(xi, yi, c="c", lw=1)
        # Draw knickpoints of the active profile bigger
        for k in self.knick_points[self.active]:
            ax.plot(xi[k[0]], yi[k[0]], k_type[k[1]], mew=0.5, mec="k", ms=10) 
        
    # SAVING METHODS
    # =============   
    def _save_regressions(self, out_file):
        
        sp = osr.SpatialReference()
        sp.ImportFromWkt(self.profiles[self.active].get_projection())
        driver = ogr.GetDriverByName("ESRI Shapefile")
        
        # Check if dataset already exists, create dataset and get layer
        if os.path.exists(out_file):
            dataset = driver.Open(out_file, 1)
            dataset.DeleteLayer(0)
            layer = dataset.CreateLayer("regresions", sp, ogr.wkbLineString)
        else:
            dataset = driver.CreateDataSource(out_file)
            layer = dataset.CreateLayer("regresions", sp, ogr.wkbLineString)
            
        # Create fields
        layer.CreateField(ogr.FieldDefn("ksn", ogr.OFTReal))
        layer.CreateField(ogr.FieldDefn("r2", ogr.OFTReal))
         
        # Save the regressions
        for idx, perfil in enumerate(self.profiles):
            for reg in self.regressions[idx]:
                xi = perfil.get_x()[reg[0]:reg[1]]
                yi = perfil.get_y()[reg[0]:reg[1]]
                feat = ogr.Feature(layer.GetLayerDefn())
                feat.SetField("ksn", reg[2])
                feat.SetField("r2", reg[5])

                # Creamos geometria
                geom = ogr.Geometry(ogr.wkbLineString)
                for n in range(len(xi)):
                    geom.AddPoint(xi[n], yi[n])
                feat.SetGeometry(geom)
                layer.CreateFeature(feat)

    def _save_knickpoints(self, out_file):
        sp = osr.SpatialReference()
        sp.ImportFromWkt(self.profiles[self.active].get_projection())
        driver = ogr.GetDriverByName("ESRI Shapefile")
    
        # Check if dataset already exists, create dataset and get layer
        if os.path.exists(out_file):
            dataset = driver.Open(out_file, 1)
            for n in range(dataset.GetLayerCount()):
                dataset.DeleteLayer(n)
            layer = dataset.CreateLayer("knickpoints", sp, ogr.wkbPoint)
        else:
            dataset = driver.CreateDataSource(out_file)
            layer = dataset.CreateLayer("knickpoints", sp, ogr.wkbPoint)

        # Create fields
        layer.CreateField(ogr.FieldDefn("ksn", ogr.OFTReal))
        layer.CreateField(ogr.FieldDefn("type", ogr.OFTInteger))
        
        # Save knickpoints
        for idx, perfil in enumerate(self.profiles):
            for kp in self.knick_points[idx]:
                xi = perfil.get_x()[kp[0]]
                yi = perfil.get_y()[kp[0]]
                ksn = perfil.get_ksn()[kp[0]]
                tipo = int(kp[1])
                
                # Populate the field ksn
                feat = ogr.Feature(layer.GetLayerDefn())
                feat.SetField("ksn", ksn)
                feat.SetField("type", tipo)
                
                # Create geometry
                geom = ogr.Geometry(ogr.wkbPoint)
                geom.AddPoint(xi, yi)
                feat.SetGeometry(geom)
                layer.CreateFeature(feat)
    
    def save(self):
        """
        Saves selected knickpoints and regressions into shapefiles
        """
        out_knicks = self.basedir + "/knickpoints.shp"
        out_regres = self.basedir + "/regressions.shp"
        self._save_regressions(out_regres)    
        self._save_knickpoints(out_knicks)
        
        if self.mode:
            self.mode = None
            self.figure.canvas.mpl_disconnect(self.pcid)
        
        np.save(self.basedir + "/graph.npy", np.array([self]))
        np.save(self.basedir + "/graph_profiles.npy", self.profiles)
        self.mode_txt.set_text('Files saved on "{0}"'.format(self.basedir))
        self.ax.figure.canvas.draw()
    
    def print_all(self, nrow, ncol, width=8.3, height=11.7):
        """
        Print all profiles in a single graphic (or graphics)
        """
        max_ksn = 0
        for perfil in self.profiles:
            local_max = perfil.get_ksn().max()
            if local_max > max_ksn:
                max_ksn = local_max
            
        label = True
        fig = plt.figure(figsize=(width, height))   # A4 size
        first = self.profiles[0].name
        last = self.profiles[0].name
        for n, perfil in enumerate(self.profiles):
            if n > (nrow * ncol):
                fig.tight_layout(pad=0.4)
                fig.savefig(self.basedir + first + "-" + last + ".eps")
                first = perfil.name
                plt.close()
                fig = plt.figure(figsize=(width, height))
                label = True
            
            ax = fig.add_subplot(nrow, ncol, n + 1)
            ax.plot(perfil.get_l() / 1000, perfil.get_z(), linewidth=0.75, color="b")
            ax.set_title(perfil.name, fontsize=12)
            ax2 = ax.twinx()
            ax2.plot(perfil.get_l()/1000, perfil.get_ksn(), linewidth=0.75, color="olive")
            ax2.set_ylim(ymax=max_ksn)
            if label:
                ax.set_ylabel("Elevation (m)")
                ax.set_xlabel("Distance (km)")
                ax2.set_ylabel("Ksn index")
                label = False
            last = perfil.name
        fig.tight_layout(pad=0.4)
        fig.savefig(self.basedir + first + "-" + last + ".eps")
        plt.close()    
