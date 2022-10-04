# -*- coding: utf-8 -*-

# network.py
# Jose Vicente Perez Pena
# Dpto. Geodinamica-Universidad de Granada
# 18071 Granada, Spain
# vperez@ugr.es // geolovic@gmail.com
#
# MIT License (see LICENSE file)
# Version: 1.1
# October 1th, 2018
#
# Last modified 07 february 2021

import numpy as np
import os
from osgeo import ogr, osr
from scipy.sparse import csc_matrix
from . import Grid, PRaster, Basin

# This import statement avoid issues with matplotlib in Mac when using Python not as a Framework
# If matplotlib is not imported, BNetwork.chi_plot() will not work.
try: 
    import matplotlib.pyplot as plt
    PLT = True
except:
    PLT = False

class Network(PRaster):
    """
    Class to manipulate cells from a Network, which is defined by applying
    a threshold to a flow accumulation raster derived from a topological 
    sorted Flow object.
    
    Parameters:
    -----------
    flow : *morfopy.Flow* object
      Flow direccion instance
    threshold : *int*
      Number the cells to initiate a channel
    thetaref : *float*
      m/n coeficient to calculate chi values in each channel cell    
    npoints : *int*
      Number of points to calculate slope and ksn in each cell. Slope and ksn values
    gradients : *bool*
      Flag to determinate if gradients are calculated or not when creating the Network 
      (calculate gradients can be slow for big grids)
    verbose : boolean
      Boolean to show processing messages in console to known the progress. Usefull with large DEMs to se the evolution.
    verb_func : str
      Function to output verbose messages (only needed if landspy is embeded in other application)
    """
    def __init__(self, flow=None, threshold=0, thetaref=0.45, npoints=5, gradients=False, verbose=False, verb_func=print):

        # The Network object has the following properties:
        # ._ix >> Giver cells
        # ._ixc >> Receivers cells
        # ._ax >> Upstream draining area (in pixel units)
        # ._dx >> Distance to mouth (distance from pixel to nearest outlet)
        # ._zx >> Elevation
        # ._chi >> Chi index
        # ._slp >> Slope of the pixel, calculated by regression with a moving window of {npoints * 2 + 1} cells.
        # ._ksn >> Ksn index of the pixel, calculated by regression with a moving window of {npoints * 2 + 1} cells.
        # ._r2slp >> R2 Coeficient of the slope regression
        # ._r2ksn >> R2 Coeficient of the ksn regression
        # ._dd >> Giver (ix) - Receiver (ixc) distance
                      
        # If flow is empty, create an empty Network instancee
        if flow is None:
            self._create_empty()
            return
        # If flow parameter is a str, we load the Network object from path
        elif type(flow)== str:
            self._load(flow)
            return
        
        # Set PRaster properties
        self._size = flow._size
        self._geot = flow._geot
        self._proj = flow._proj
     
        # Get a threshold if not specified (Default 0.25% of the total number of cells)
        if threshold == 0:
            threshold = self.getNCells() * 0.0025
        self._threshold = int(threshold)
        
        # Get sort Nodes for channel cells and elevations
        fac = flow.flowAccumulation(nodata=False, asgrid=False)
        w = fac >= threshold
        w = w.ravel()
        I   = w[flow._ix]
        self._ix  = flow._ix[I]
        self._ixc = flow._ixc[I]
        
        # Get Area and Elevations for channel cells
        self._ax = fac.ravel()[self._ix] # Area in CELLS units!!
        self._zx = flow._zx[I]
        
        # Get distances to mouth (self._dx) and giver-receiver distances (self._dd)
        di = np.zeros(self.getNCells())
        self._dd = np.zeros(self._ix.shape) # Giver-Receiver distance
        for n in np.arange(self._ix.size)[::-1]:
            grow, gcol = self.indToCell(self._ix[n])
            rrow, rcol = self.indToCell(self._ixc[n])
            gx, gy = self.cellToXY(grow, gcol)
            rx, ry = self.cellToXY(rrow, rcol)
            d_gr = np.sqrt((gx - rx)**2 + (gy - ry)**2)
            self._dd[n] = d_gr
            di[self._ix[n]] = di[self._ixc[n]] + d_gr
        self._dx = di[self._ix]
        
        # Get chi values using the input thetaref
        self._thetaref = thetaref
        self.calculateChi(thetaref)
        
        # Calculate slopes and ksn
        if gradients:
            self.calculateGradients(npoints, 'slp')
            self.calculateGradients(npoints, 'ksn')
        else:
            self._slp = np.zeros(self._ix.size)
            self._r2slp = np.zeros(self._ix.size)
            self._slp_np = 0
            self._ksn = np.zeros(self._ix.size)
            self._r2ksn = np.zeros(self._ix.size)
            self._ksn_np = 0

    def _create_empty(self):
        """
        Creates a empty instance of the Network class
        """
        # Set the empty PRaster properties
        super().__init__()
        # Set remaining properties
        self._ix =  np.array([0])
        self._ixc = np.array([0])
        self._ax = np.array([1])
        self._dx = np.array([0])
        self._zx = np.array([0])
        self._chi = np.array([0])
        self._slp = np.array([0])
        self._ksn = np.array([0])
        self._r2slp = np.array([0])
        self._r2ksn = np.array([0])
        self._dd = np.array([1])
        self._thetaref = 0.45
        self._slp_np = 0
        self._ksn_np = 0
        self._threshold = 0
        
    def save(self, path):
        """
        Saves the Network instance to disk. It will be saved as a numpy array in text format with a header.
        The first three lines will have the information of the raster:
            Line1::   xsize; ysize; cx; cy; ULx; ULy; Tx; Ty
            Line2::   thetaref; threshold; slp_np; ksn_np
            Line3::   String with the projection (WKT format)
        xsize, ysize >> Dimensions of the raster
        cx, cy >> Cellsizes in X and Y
        Tx, Ty >> Rotation factors (for geotransformation matrix)
        ULx, ULy >> X and Y coordinates of the corner of the upper left pixel of the raster
        thetaref >>  m/n coeficient to calculate chi values in each channel cell
        threshold >> Number the cells to initiate a channel
        slp_np, ksn_np >> Number of points to calculate ksn and slope by regression. Window of {npoints * 2 + 1}
        
        Parameters:
        ===========
        path : *str*
          Path to save the network object with *.dat extension (it is not necessary to  give the extension)
        """
    
        # In case the extension is wrong or path has not extension
        path = os.path.splitext(path)[0] + ".dat"
        
        # Create header with properties
        params = [self._size[0], self._size[1], self._geot[1], self._geot[5], 
                  self._geot[0], self._geot[3], self._geot[2], self._geot[4]]
        header = ";".join([str(param) for param in params]) + "\n"  
        params = [self._thetaref, self._threshold, self._slp_np, self._ksn_np]
        header += ";".join([str(param) for param in params]) + "\n" 
        header += str(self._proj)

        # Create data array
        data_arr = np.array((self._ix, self._ixc, self._ax, self._dx, self._zx,
                             self._chi, self._slp, self._ksn, self._r2slp, 
                             self._r2ksn, self._dd)).T
        
        # Save the network instance as numpy.ndarray in text format
        np.savetxt(path, data_arr, delimiter=";", header=header, encoding="utf8", comments="#")
    
    def _load(self, path):
        """
        Loads a Network instance saved in the disk.
        
        Parameter:
        ==========
           Path to the saved network object
        """
        # Open the file as normal text file to get its properties
        fr = open(path, "r")
        # Line 1: First and last characters will be "#" and "\n"
        linea = fr.readline()[1:-1]
        data = linea.split(";")
        self._size = (int(data[0]), int(data[1]))
        self._geot = (float(data[4]), float(data[2]), float(data[6]), 
                      float(data[5]), float(data[7]), float(data[3]))       
        # Line2: First and last characters will be "#" and "\n"
        linea = fr.readline()[1:-1]
        data = linea.split(";")
        self._thetaref = float(data[0])
        self._threshold = int(data[1])
        self._slp_np = int(data[2])
        self._ksn_np = int(data[3])
        # Line3: First and last characters will be "#" and "\n"
        linea = fr.readline()[1:-1]
        self._proj = linea
        fr.close()
        
        # Load array data
        data_arr = np.loadtxt(path, dtype=float, comments='#', delimiter=";", encoding="utf8")
        # Fix to avoid errors in networks with only one cell...
        if data_arr.ndim < 2:
            data_arr = data_arr.reshape((1, data_arr.size))
        self._ix = data_arr[:, 0].astype(np.int64)
        self._ixc = data_arr[:, 1].astype(np.int64)
        self._ax = data_arr[:, 2]
        self._dx = data_arr[:, 3]
        self._zx = data_arr[:, 4]
        self._chi = data_arr[:, 5]
        self._slp = data_arr[:, 6]
        self._ksn = data_arr[:, 7]
        self._r2slp = data_arr[:, 8]
        self._r2ksn = data_arr[:, 9]
        self._dd = data_arr[:, 10]
         
    def calculateChi(self, thetaref=0.45, a0=1.0):
        """
        Function that calculates chi_values for channel cells
        
        Parameters:
        -----------
        thetaref : *float*
          m/n coeficient to calculate chi
        a0 : *float*
          Reference area to avoid dimensionality (usually don't need to be changed)
        """
        chi = np.zeros(self.getNCells())
        for n in np.arange(self._ix.size)[::-1]:
            chi[self._ix[n]] = chi[self._ixc[n]] + (a0 * self._dd[n]/self._ax[n]**thetaref)            
        self._chi = chi[self._ix]
        self._thetaref = thetaref
   
    def polynomial_fit(self, x, y):
        '''Calculate gradient and R2 for two variables''' 
       
        # Calculate slope of central cell by regression 
        poli, SCR = np.polyfit(x, y, deg = 1, full = True)[:2]
        # Calculate gradient
        g = poli[0]
        if g < 0.001: 
            g = 0.001
    
        # Calculate R2 
        if y.size * y.var() == 0:
            R2 = 1 # Puntos colineares
        else:
            R2 = float(1 - SCR/(y.size * y.var()))
    
        return (g, R2)  
    
    def calculateGradients(self, npoints, kind='slp'):
        """
        This function calculates gradients (slope or ksn) for all channel cells. 
        Gradients of each cell are calculated by linear regression using a number
        of points (npoints) up and downstream.
        
        Parameters:
        ===========
        npoints : *int*
          Window to analyze slopes. Slopes are calculated by linear regression using a window
          of (npoints * 2 + 1) pixels (using the central pixel)
          
        kind : *str* {'slp', 'ksn'}
          Kind of gradient to calculate Slope (slp) or Ksn (ksn)
        """

        winlen = npoints * 2 + 1
        
        # Get arrays depending on type
        y_arr = self._zx
        if kind == 'ksn':
            x_arr = self._chi
        else:
            x_arr = self._dx
            
        # Get ixcix auxiliar array
        ixcix = np.zeros(self.getNCells(), np.int64)
        ixcix[self._ix] = np.arange(self._ix.size)
        
        # Get heads array and confluences dictionary
        heads = self.streamPoi("heads", "IND")
        #confs = {conf:[] for conf in self.get_stream_poi("confluences", "IND")}
        
        # Sort heads by elevation
        elev = self._zx[ixcix[heads]]
        spos = np.argsort(-elev)
        heads = heads[spos]
        
        # Prepare auxiliary arrays
        gi = np.zeros(self.getNCells())
        r2 = np.zeros(self.getNCells())
        
        # Taking sequentally all the heads and compute downstream flow
        for head in heads:
            processing = True
            head_cell = head
            mid_cell = self._ixc[ixcix[head_cell]]
            mouth_cell = self._ixc[ixcix[mid_cell]]
            win = [mouth_cell, mid_cell, head_cell]
            
            if ixcix[mid_cell] == 0 or gi[mid_cell] != 0:
                # Channel type 2 (mid_cell is an outlet) or mid_cell is already calculated
                continue
            
            elif ixcix[mouth_cell]== 0:
                # Channel type 3 (mouth_cell is an outlet)
                processing = False
                win = [mid_cell, mid_cell, head_cell]
                
            # Obtenemos datos de elevacion y distancias
            xi = x_arr[ixcix[win]]
            yi = y_arr[ixcix[win]]
           
            g, R2 = self.polynomial_fit(xi, yi)
        
            gi[mid_cell] = g
            r2[mid_cell] = R2
            
            in_outlet = False
                
            while processing:
                # Verificamos si estamos al final (en un outlet)  
                if not in_outlet:
                    # Si mouth_cell no es un outlet, cogemos siguiente celda
                    next_cell = self._ixc[ixcix[mouth_cell]]
                    # Si la siguiente celda es el final, 
                    # añadimos una celda y eliminamos la cabecera
                    if ixcix[next_cell]==0:
                        win.insert(0, next_cell)
                        win.pop()
                    # Si longitud de ventana < winlen, se añaden dos celdas 
                    elif len(win) < winlen:
                        win.insert(0, next_cell)
                        aux_cell = self._ixc[ixcix[next_cell]]
                        win.insert(0, aux_cell)
                    else:
                        win.insert(0, next_cell)
                        win.pop() 
                else:
                    # Si mouth_cell es un outlet, no se coge siguiente celda
                    next_cell = mouth_cell
                    win.pop()
                    win.pop()
                
                head_cell = win[-1]
                mid_cell = self._ixc[ixcix[mid_cell]]
                mouth_cell = win[0]
            
                # Verificamos si mouth_cell es un outlet
                if ixcix[mouth_cell] == 0:
                    in_outlet = True
                    # Si mouth_cell es un outlet, hay que eliminarla de la ventana
                    # puesto que ixcix[mouth_cell] será la posición 0 del array ._ix
                    win[0] = win[1]    
            
                if len(win) <= 3:
                    processing = False
                    
                # Obtenemos datos de elevacion y distancias para calcular pendientes
                xi = x_arr[ixcix[win]]
                yi = y_arr[ixcix[win]]
                
                g, R2 = self.polynomial_fit(xi, yi)
                            
                # Comprobamos si celda central ha sido procesada
                if gi[mid_cell] == 0:
                    gi[mid_cell] = g
                    r2[mid_cell] = R2
                else:
                    # Si ha sido procesada, es una confluencia 
                    processing = False
        #             if len(confs[mid_cell]) == 0:
        #                 confs[mid_cell].append(gi[mid_cell])
        #                 confs[mid_cell].append(g)
        #             else:
        #                 confs[mid_cell].append(g)
                        
        # # Calculamos valores medios en confluencias                
        # for cell in confs.keys():
        #     if len(confs[cell]) > 0:
        #         gi[cell] = np.mean(np.array(confs[cell]))
                
        # Llenamos array del objeto Network
        if kind == 'ksn':
            self._ksn = gi[self._ix]
            self._r2ksn = r2[self._ix]
            self._ksn_np = npoints
        else:
            self._slp = gi[self._ix]
            self._r2slp = r2[self._ix]
            self._slp_np = npoints

    def hierarchy_channels(self, heads="", asgrid=True):
        """
        This function classifies channels into major and minor. By default this
        order is calculated according to the height of the heads. 
        
        Parameters
        ----------
        heads : *list*
            List with indices of the heads of the channels arranged according 
            to the hierarchy of the channels.
            
        asgrid : *bool*
          Indicates if the function is returned as landspy.Grid (True) or as a 
          numpy.array (False)
          
        """
        
        tier=0
        # Prepare auxiliary arrays
        ladder = np.zeros(self.getNCells())
        # Get ixcix auxiliar array
        ixcix = np.zeros(self.getNCells(), np.int64)
        ixcix[self._ix] = np.arange(self._ix.size)   
        
        if heads == "":
            # Sort heads by elevation
            heads = self.streamPoi("heads", "IND")
            elev = self._zx[ixcix[heads]]
            spos = np.argsort(-elev)
            heads = heads[spos]
        
        # Taking sequentally all the heads and compute downstream flow
        for head in heads:
            tier+=1
            processing = True
            ladder[head]=tier
            head_cell = head
            next_cell = self._ixc[ixcix[head_cell]]
            
            while processing:
                if ladder[next_cell] > 0 or ixcix[next_cell] == 0:
                    processing = False
                    continue
                ladder[next_cell]=tier
                next_cell = self._ixc[ixcix[next_cell]]
                
        ladder = ladder.reshape(self.getDims())
        
        # Return grid
        if asgrid:
            return self._create_output_grid(ladder, 0)
        else:
            return ladder
        
    def streamPoi(self, kind="heads", coords="CELL"):
            """
            This function finds points of interest of the drainage network. These points of interest
            can be 'heads', 'confluences' or 'outlets'.
            
            Parameters:
            -----------
            kind : *str* {'heads', 'confluences', 'outlets'}
              Kind of point of interest to return.
            coords : *str* {'CELL', 'XY', 'IND'}
              Output coordinates for the stream point of interest. 
              
            Returns:
            -----------
            numpy.ndarray
              Numpy ndarray with one (id) or two columns ([row, col] or [xi, yi] - depending on coords) 
              with the location of the points of interest 
              
            References:
            -----------
            The algoritms to extract the point of interest have been adapted to Python 
            from Topotoolbox matlab codes developed by Wolfgang Schwanghart (version of 17. 
            August, 2017). These smart algoritms use sparse arrays with giver-receiver indexes, to 
            derive point of interest in a really efficient way. Cite:
                    
            Schwanghart, W., Scherler, D., 2014. Short Communication: TopoToolbox 2 - 
            MATLAB-based software for topographic analysis and modeling in Earth 
            surface sciences. Earth Surf. Dyn. 2, 1–7. https://doi.org/10.5194/esurf-2-1-2014
            """

            # Get grid channel cells
            w = np.zeros(self.getNCells(), dtype="bool")
            w[self._ix] = True
            w[self._ixc] = True
            
            # Build a sparse array with giver-receivers cells
            aux_vals = np.ones(self._ix.shape, dtype=np.int8)
            sp_arr = csc_matrix((aux_vals, (self._ix, self._ixc)), shape=(self.getNCells(), self.getNCells()))
            
            # Get stream POI according the selected type
            if kind == 'confluences':
                # Confluences will be channel cells with two or givers
                sum_arr = np.asarray(np.sum(sp_arr, 0)).ravel()
                out_pos = sum_arr > 1
            elif kind == 'outlets':
                # Outlets will be channel cells marked only as receivers (ixc) but not as givers (ix) 
                sum_arr = np.asarray(np.sum(sp_arr, 1)).ravel()
                out_pos = np.logical_and((sum_arr == 0), w)  
            else:
                # Heads will be channel cells marked only as givers (ix) but not as receivers (ixc) 
                sum_arr = np.asarray(np.sum(sp_arr, 0)).ravel()
                out_pos = (sum_arr == 0) & w
                
            out_pos = out_pos.reshape(self.getDims())
            row, col = np.where(out_pos)
            
            if coords=="XY":
                xi, yi = self.cellToXY(row, col)
                return np.array((xi, yi)).T
            elif coords=="IND":
                return self.cellToInd(row, col)
            else:
                return np.array((row, col)).T

    def snapPoints(self, input_points, kind="channel", remove_duplicates=False):
        """
        Snap input points to channel cells or to stream POI
        
        Parameters:
        ===========
        input_points : *numpy.ndarray*
          Numpy 2-D ndarray, which first two columns are x and y coordinates [x, y, ...]
        kind : *str* {'channel', 'heads', 'confluences', 'outlets'}  
          Kind of point to snap input points
        remove_duplicates : *bool*
          Remove duplicate points. When snapping points, two points can be snapped to the same poi. If True, 
          these duplicate points will be removed. 
        
        Returns:
        ===========
        numpy.ndarray
          Numpy ndarray with two columns [xi, yi] with the snap points
        """
        
        # Extract a numpy array with the coordinate to snap the points
        if kind in ['heads', 'confluences', 'outlets']:
            poi = self.streamPoi(kind, "XY")     
        else:
            row, col = self.indToCell(self._ix)
            x, y = self.cellToXY(row, col)
            poi = np.array((x, y)).T
        
        # Get array reshaped for the calculation
        xi = input_points[:, 0].reshape((input_points.shape[0], 1))
        yi = input_points[:, 1].reshape((input_points.shape[0], 1))
        xci = poi[:, 0].reshape((1, poi.shape[0]))
        yci = poi[:, 1].reshape((1, poi.shape[0]))
        
        # Calculate distances and get minimum
        di = np.sqrt((xi - xci)**2 + (yi - yci)**2 )
        pos = np.argmin(di, axis=1)
        
        # Get the rest of columns (if the case)
        if input_points.shape[1] > 2:
            aux = input_points[:, 2:]
            out_p = np.concatenate((poi[pos], aux, pos.reshape(pos.size, 1)), axis=1)
        else:
            out_p = poi[pos]
    
        # Remove duplicates
        if remove_duplicates:     
            idx = np.unique(out_p[:,-1], return_index=True)[1]
            out_p = out_p[idx, :-1]
    
        return out_p
    
    def exportPoints(self, path):
        """
        Export channel points to a semicolon-delimited text file
        This file will contain  data of 
        id ; x ; y ; z ; distance ; area ; chi ; slope ; ksn ; r2_slope ; r2_ksn
        
        path : str
          Path for the output text file
        """
        cab = "id;x;y;z;distance;area;chi;slope;ksn;r2_slope;r2_ksn"
        row, col = self.indToCell(self._ix)
        x, y = self.cellToXY(row, col)
        
        out_arr = np.array((self._ix, x, y, self._zx, self._dx, self._ax, self._chi, 
                            self._slp, self._ksn, self._r2slp, self._r2ksn)).T
        np.savetxt(path, out_arr, delimiter=";", header=cab, comments="", encoding="utf8")
    
    def getStreams(self, asgrid=True):
        """
        This function outputs a grid representation of the Network object

        Parameters:
        ===========
        asgrid : *bool*
          Indicates if the network is returned as landspy.Grid (True) or as a numpy.array
        """
        # Get grid channel cells
        w = np.zeros(self.getNCells(), dtype=np.int8)
        w[self._ix] = 1
        w[self._ixc] = 1
        w = w.reshape(self.getDims())
        # Return grid
        if asgrid:
            return self._create_output_grid(w, 0)
        else:
            return w
    
    def getStreamSegments(self, asgrid=True):
        """
        This function extract a drainage network by using a determined area threshold
        and output the numerated stream segments.

        Parameters:
        ===========
        asgrid : *bool*
          Indicates if the network is returned as landspy.Grid (True) or as a numpy.array
        """
        # Get heads and confluences and merge them
        head_ind = self.streamPoi("heads", "IND")
        conf_ind = self.streamPoi("confluences", "IND")
        all_ind = np.append(head_ind, conf_ind)
        del conf_ind, head_ind # Clean up
        
        # We created a zeros arrays and put in confluences and heads their id
        # Those id will be consecutive numbers starting in one
        seg_arr = np.zeros(self.getNCells(), dtype=np.int32)
        for n, inds in enumerate(all_ind):
            seg_arr[inds] = n+1
        
        # Move throught channel list. If receiver is 0, give receiver the same id that giver.
        # If a receiver is not 0, that means that we are in a confluence. 
        for n in range(len(self._ix)):
            if seg_arr[self._ixc[n]] == 0:
                seg_arr[self._ixc[n]] = seg_arr[self._ix[n]]
        
        # Reshape and output
        seg_arr = seg_arr.reshape(self.getDims())
        if asgrid:
            return self._create_output_grid(seg_arr, 0)
        else:
            return seg_arr
        
    def getStreamOrders(self, kind="strahler", asgrid=True):
        """
        This function extract streams orderded by strahler or shreeve. Cell values
        will have a value acording with the order of the segment they belong
    
        Parameters:
        ===========
        kind : *str* {'strahler', 'shreeve'}
        asgrid : *bool*
          Indicates if the selfwork is returned as landspy.Grid (True) or as a numpy.array
        """
        if kind not in ['strahler', 'shreeve']:
            kind = 'strahler'
        
        # Get grid channel cells
        str_ord = np.zeros(self.getNCells(), dtype=np.int64)
        str_ord[self._ix] = 1
        str_ord[self._ixc] = 1
        visited = np.zeros(self.getNCells(), dtype=np.int64)
    
        if kind == 'strahler':
            for n in range(len(self._ix)):
                if (str_ord[self._ixc[n]] == str_ord[self._ix[n]]) & visited[self._ixc[n]]:
                    str_ord[self._ixc[n]] = str_ord[self._ixc[n]] + 1                
                else:
                    str_ord[self._ixc[n]] = max(str_ord[self._ix[n]], str_ord[self._ixc[n]])
                    visited[self._ixc[n]] = True
        elif kind == 'shreeve':
            for n in range(len(self._ix)):
                if visited[self._ixc[n]]:
                    str_ord[self._ixc[n]] = str_ord[self._ixc[n]] + str_ord[self._ix[n]]
                else:
                    str_ord[self._ixc[n]] = max(str_ord[self._ix[n]], str_ord[self._ixc[n]])
                    visited[self._ixc[n]] = True
        str_ord = str_ord.reshape(self.getDims())
        
        if asgrid:
            return self._create_output_grid(str_ord, nodata_value=0)
        else:
            return str_ord

    def exportShp(self, path, con=False):
        """
        Export Network channels to shapefile format.
        
        path : str
          Path to save the shapefile
        con : bool
          If False, channels will split in each confluence (segmented channels). If True, 
          they will split only when order changes (continuous channels).
        """
        if con:
            self._get_continuous_shp(path)
        else:
            self._get_segmented_shp(path)
    
    def _get_segmented_shp(self, path=""):
        """
        Export Network channels to shapefile format. Channels will split in each confluence.
        
        path : str
          Path to save the shapefile 
        """
        # Create shapefile
        driver = ogr.GetDriverByName("ESRI Shapefile")
        dataset = driver.CreateDataSource(path)
        sp = osr.SpatialReference()
        sp.ImportFromWkt(self._proj)
        layer = dataset.CreateLayer("rivers", sp, ogr.wkbLineString25D)
        layer.CreateField(ogr.FieldDefn("segid", ogr.OFTInteger))
        layer.CreateField(ogr.FieldDefn("strahler", ogr.OFTInteger))
        layer.CreateField(ogr.FieldDefn("shreeve", ogr.OFTInteger)) 
        layer.CreateField(ogr.FieldDefn("flowto", ogr.OFTInteger))
        
        # Get channel segments and orders
        ch_seg = self.getStreamSegments(False).ravel()
        ch_stra = self.getStreamOrders(asgrid=False).ravel()
        ch_shre = self.getStreamOrders(kind='shreeve', asgrid=False).ravel()  
        ch_seg = ch_seg[self._ix]
        ch_stra = ch_stra[self._ix]
        ch_shre = ch_shre[self._ix]  
        
        # Get ixcix auxiliar array
        ixcix = np.zeros(self.getNCells(), np.int64)
        ixcix[self._ix] = np.arange(self._ix.size)
        
        seg_ids = np.unique(ch_seg)
        for idx in seg_ids:
            # skip zero channel (no channel)
            if idx == 0:
                continue   
            # Get givers for the segment id
            pos = np.where(ch_seg == idx)[0]
            ch_ix = self._ix[pos]
            
            # Add last point
            first = ch_ix[0]
            last = self._ixc[ixcix[ch_ix[-1]]]
            ch_ix = np.append(ch_ix, last)
            first = ch_ix[0]
            
            # Get segment order and receiver segment
            stra = ch_stra[ixcix[first]]
            shre = ch_shre[ixcix[first]]
            if ixcix[last] == 0:
                flowto = idx
            else:
                flowto = ch_seg[ixcix[last]]
            
            # Add feature
            feat = ogr.Feature(layer.GetLayerDefn())
            feat.SetField("segid", int(idx))
            feat.SetField("strahler", int(stra))
            feat.SetField("shreeve", int(shre))
            feat.SetField("flowto", int(flowto))
            row, col = self.indToCell(ch_ix)
            xi, yi = self.cellToXY(row, col)
            pos = ixcix[ch_ix]
            zi = self._zx[pos]
            
            geom = ogr.Geometry(ogr.wkbLineString25D)
            
            for n in range(xi.size):
                geom.AddPoint(xi[n], yi[n], zi[n])
                
            feat.SetGeometry(geom)
            layer.CreateFeature(feat)
           
        layer = None
        dataset = None

    def _get_continuous_shp(self, path=""):
        """
        Export Network channels to shapefile format. Channels will split only when order changes.
        
        path : str
          Path to save the shapefile 
        """
        # Create shapefile
        driver = ogr.GetDriverByName("ESRI Shapefile")
        dataset = driver.CreateDataSource(path)
        sp = osr.SpatialReference()
        sp.ImportFromWkt(self._proj)
        layer = dataset.CreateLayer("rivers", sp, ogr.wkbLineString25D)
        layer.CreateField(ogr.FieldDefn("order", ogr.OFTInteger))
        
        # Get ixcix auxiliar array
        ixcix = np.zeros(self.getNCells(), np.int64)
        ixcix[self._ix] = np.arange(self._ix.size)
        
        # Get heads, confluences, outlets and orders
        heads = self.streamPoi("heads", "IND")
        confs = self.streamPoi("confluences", "IND")
        ch_ord = self.getStreamOrders(asgrid=False).ravel()
        
        # Get confluences where strahler index increases
        strahler_confs = []
        for conf in confs:
            conf_order = ch_ord[conf]
            givers = self._ix[np.where(self._ixc==conf)]
            giv_orders = ch_ord[givers]
            if giv_orders.max() < conf_order:
                strahler_confs.append(conf)
                
        # Create head array with orders (by append confluences where index increases)
        heads = np.array((heads, np.ones(heads.size, np.int64))).T
        strahler_confs = np.array([strahler_confs, ch_ord[strahler_confs]]).T
        heads = np.vstack((heads, strahler_confs))
        
        # Iterate heads
        for head in heads:
            cell = head[0]
            if ixcix[cell] == 0: # Head is also an outlet, discard and continue
                continue
            river_data = [cell]
            processing = True
            while processing:
                next_cell = self._ixc[ixcix[cell]]
                river_data.append(next_cell)
                           
                if ixcix[next_cell] == 0: # next_cell is an outlet
                    processing = False
                elif next_cell in confs: # next_cell is a confluence
                    if head[1] < ch_ord[next_cell]:
                        processing = False
                    else:
                        cell = next_cell
                else:
                    cell = next_cell    
                       
            # Add feature
            feat = ogr.Feature(layer.GetLayerDefn())
            row, col = self.indToCell(river_data)
            xi, yi = self.cellToXY(row, col)
            pos = ixcix[river_data]
            zi = self._zx[pos]            
            
            geom = ogr.Geometry(ogr.wkbLineString25D)
            for n in range(xi.size):
                geom.AddPoint(xi[n], yi[n], zi[n])
                
            feat.SetGeometry(geom)
            chanorder = ch_ord[cell]
            feat.SetField("order", int(chanorder))
            layer.CreateFeature(feat)
           
        layer = None
        dataset = None

    def _create_output_grid(self, array, nodata_value=None):
        """
        Convenience function that creates a Grid object from an input array. The array
        must have the same shape that self._dims and will maintain the Flow object 
        properties as dimensions, geotransform, reference system, etc.
        
        Parameters:
        ===========
        array : *numpy.ndarray*
          Array to convert to a Grid object
        nodata_value _ *int* / *float*
          Value for NoData values
          
        Returns:
        ========
        Grid object with the same properties that Flow
        """
        grid = Grid()
        grid.copyLayout(self)
        grid._nodata = nodata_value
        grid._array = array
        grid._tipo = str(array.dtype)
        return grid

    def getChannel(self, head, mouth=None, name="", oid=-1):
        """
        Get a channel from the head to the mouth. 
        
        Parameters
        ----------
        head : tuple
            (x, y) tuple with the coordinates of the channel head (will be snapped to a channel cell)
        mouth : tuple
            (x, y) tuple with the coordinates of the channel mouth (will be snapped to a channel cell).
            If None or a cell out of the current channel, the channel will continue until the nearest outlet.
        name : str, optional
            Name ("label") for the channel    
        oid : int, optional
            Id of the channel
            
        Returns
        -------
        Channel instance

        """
        if mouth is None:
            mouth = head
            
        # Check if head and mouth are inside the grid
        if not self.isInside(head[0], head[1]):
            return
        elif not self.isInside(mouth[0], mouth[1]):
            mouth = head
        
        # Snap head and mouth to channel cells
        snap_points = self.snapPoints(np.array((head, mouth)))
        row, col = self.xyToCell(snap_points[:,0], snap_points[:,1])
        idx = self.cellToInd(row, col)
        head = idx[0]
        mouth = idx[1]
        
        # Get ixcix auxiliar array
        ixcix = np.zeros(self.getNCells(), np.int64)
        ixcix[self._ix] = np.arange(self._ix.size)
        
        # Get channel cells
        chcells = [head]
        next_cell = self._ixc[ixcix[head]]
        while ixcix[next_cell] != 0:
            chcells.append(next_cell)
            if next_cell == mouth:
                break
            next_cell = self._ixc[ixcix[next_cell]]

        chcells = np.array(chcells)
        auxarr = np.zeros(self.getNCells()).astype("bool")
        auxarr[chcells] = True
        I = auxarr[self._ix]
        ax = self._ax[I]
        dx = self._dx[I]
        zx = self._zx[I]
        chi = self._chi[I]
        slp = self._slp[I]
        ksn = self._ksn[I]
        r2_slp = self._r2slp[I]
        r2_ksn = self._r2ksn[I]
        dd = self._dd[I]
        
        chandata = np.array([chcells, ax, dx, zx, chi, slp, ksn, r2_slp, r2_ksn, dd]).T
        return Channel(self, chandata, self._thetaref, self._chi[-1], self._slp_np, self._ksn_np, name=name, oid=oid)
        
    def chiShapefile(self, out_shp, distance):
        """
        This method export network data to a shapelife. It calculates segments of a given
        distance and calculate chi, ksn, slope, etc. for the segment. The shapefile will 
        have the following fields:
            id_profile : Profile identifier. Profiles are calculated from heads until outlets or  
            L : Lenght from the middle point of the segment to the profile head
            area_e6 : Drainage area in the segment mouth (divided by E6, to avoide large numbers) 
            z : Elevation of the middle point of the segment
            chi : Mean chi of the segment 
            ksn : Ksn of the segment (calculated by linear regression)
            slope : slope of the segment (calculated by linear regression)
            rksn : R2 of the ksn linear regression
            rslope : R2 of the slope linear regression
            
        Parameters:
        ===========
        out_shp : srt
          Output shapefile
        distance : float
          Segment distance. 
            
        """
        # Create shapefile
        driver = ogr.GetDriverByName("ESRI Shapefile")
        dataset = driver.CreateDataSource(out_shp)
        sp = osr.SpatialReference()
        sp.ImportFromWkt(self.getCRS())
        layer = dataset.CreateLayer("rivers", sp, ogr.wkbLineString25D)
        
        # Add fields
        campos = ["id_profile", "L", "area_e6", "z", "chi", "ksn", "rksn", "slope", "rslope"]
        tipos = [0, 2, 2, 2, 2, 2, 2, 2, 2]
        for n in range(len(campos)):
            layer.CreateField(ogr.FieldDefn(campos[n], tipos[n]))
        
        # Get ixcix auxiliar array
        ixcix = np.zeros(self.getNCells(), np.int64)
        ixcix[self._ix] = np.arange(self._ix.size)
        
        # Get heads and sort them by elevation and iterate them
        heads = self.streamPoi("heads", "IND")
        zpos = np.argsort(self._zx[ixcix[heads]])
        heads = heads[zpos][::-1]
        
        aux_arr = np.zeros(self.getNCells(), "bool")
        id_profile = 0
        for head in heads:
            id_profile += 1
            processing = True
            cell = head
            segment_cells = [cell]
            segment_distance = 0
            profile_length = self._dx[ixcix[head]]
            
            while processing:
                add_segment = False
                # Take the next cell downstream and add it to the array
                next_cell = self._ixc[ixcix[cell]]
                segment_cells.append(next_cell)
                segment_distance += self._dd[ixcix[cell]]

                if segment_distance >= distance:
                    # If segment distance is reached, add the current segment
                    add_segment = True            
                
                if aux_arr[next_cell] == True or ixcix[next_cell] == 0:
                    # If the cell is processed or an outlet
                    # End processing and add the segment (even if distance not reached)
                    processing = False
                    add_segment = True
                            
                aux_arr[next_cell] = True
                cell = next_cell
        
                # Add segment to the shapefile (if has more than 3 cells)
                if add_segment and len(segment_cells) >= 3:
                    # Get point coordinates and values
                    row, col = self.indToCell(segment_cells)
                    xi, yi = self.cellToXY(row, col)
                    pos = ixcix[segment_cells][::-1]
                    zi = self._zx[pos[::-1]]
                    mouth_cell = ixcix[segment_cells[-1]]
                    mid_cell = ixcix[segment_cells[int(len(segment_cells)/2)]]
                    dx = self._dx[pos]
                    zx = self._zx[pos]
                    chi = self._chi[pos]
                    area = self._ax[mouth_cell]
                    lenght = profile_length - self._dx[mid_cell]
                    mid_z = self._zx[mid_cell]
                    mid_chi = self._chi[mid_cell]      
                    
                    # Mean slope for the segment (by min. sqares)
                    poli, SCR = np.polyfit(dx, zx, deg = 1, full = True)[:2]
                    slope = poli[0]
                    if slope == 0: 
                        slope = 0.000001
                    if zx.size * zx.var() == 0:
                        rslope = 1
                    else:
                        rslope = float(1 - SCR/(zx.size * zx.var()))      
                    
                    # Mean ksn for the segment (by min. sqares)
                    poli, SCR = np.polyfit(chi, zx, deg = 1, full = True)[:2]
                    ksn = poli[0]
                    if ksn == 0: 
                        ksn = 0.000001
                    if zx.size * zx.var() == 0:
                        rksn = 1
                    else:
                        rksn = float(1 - SCR/(zx.size * zx.var()))
                    
                    # Create feature and add attributes
                    feat = ogr.Feature(layer.GetLayerDefn())
                    feat.SetField("id_profile", int(id_profile))
                    feat.SetField("L", float(lenght))
                    feat.SetField("area_e6", float(area/1000000))
                    feat.SetField("z", float(mid_z))
                    feat.SetField("chi", float(mid_chi))
                    feat.SetField("ksn", float(ksn))
                    feat.SetField("rksn", float(rksn))
                    feat.SetField("slope", float(slope))
                    feat.SetField("rslope", float(rslope))
                    
                    # Create geometry
                    geom = ogr.Geometry(ogr.wkbLineString25D)
                    for n in range(xi.size):
                        geom.AddPoint(xi[n], yi[n], zi[n])
                    feat.SetGeometry(geom)
                    # Add segment feature to the shapefile
                    layer.CreateFeature(feat)
                    # Reset variables
                    segment_cells = [next_cell]
                    segment_distance = 0
        
        layer = None
        dataset = None


class BNetwork(Network):
    """
    Class to manipulate cells from a drainage network from a single basin network. 
    This class inhereits all methods and properties of the Network class plus some 
    new methods to manipulate channels
    
    Parameters:
    -----------
    net : *landspy.Network* | *str*
      Network instance or path to a previously saved BNetwork file
    basingrid : *landspy.Basin* | *landspy.Grid*
      Numpy array or landspy Grid representing the drainage basin. If array or Grid have more than
      one basin, set the basinid properly. 
    heads : *list* or *numpy.ndarray*
      List with [x, y] coordinates for basin the main head or 2 column numpy.ndarray with head(s) 
      coordinate(s). If more than one, the first one is considered the main head (trunk channel). 
      Head(s) will be snapped to Network heads.
    bid : *int*
      Id value that identifies the basin cells in case that basin will have more than one basin.        
     """    
    def __init__(self, net, basingrid=None, heads=None, bid=1):

        # If flow is a str, load it
        if isinstance(net, str):
            self._load(net)
        
        else: 
            # Get a basin grid with net dimensions (basin cells will set to 1)
            if isinstance(basingrid, Basin):
                basin = np.zeros(net.getDims(), dtype=np.int8)
                x = basingrid.getGeot()[0] + (basingrid.getGeot()[1]/2)
                y = basingrid.getGeot()[3] + (basingrid.getGeot()[5]/2)
                row, col = net.xyToCell(x, y)
                arr = np.where(basingrid.readArray()==basingrid.getNodata(), 0, 1)
                basin[row:row+basingrid.getDims()[0], col:col+basingrid.getDims()[1]] = arr        
            
            elif isinstance(basingrid, Grid):
                basin = np.where(basingrid.readArray()==bid, 1, 0)
                basingrid = basingrid.copy()
                basingrid.setArray(basin)
                basingrid.setNodata(0)
            
            # Get limits for the input basin
            c1 = basin.max(axis=0).argmax()
            r1 = basin.max(axis=1).argmax()
            c2 = basin.shape[1] - np.fliplr(basin).max(axis=0).argmax()
            r2 = basin.shape[0] - np.flipud(basin).max(axis=1).argmax()
            
            # Cut basin by those limits
            basin_cl = basin[r1:r2, c1:c2]
            
            # Create Grid
            self._size = (basin_cl.shape[1], basin_cl.shape[0])
            self._dims = (basin_cl.shape[0], basin_cl.shape[1])
            geot = net._geot
            ULx = geot[0] + geot[1] * c1
            ULy = geot[3] + geot[5] * r1
            self._geot = (ULx, geot[1], 0.0, ULy, 0.0, geot[5])
            self._cellsize = (geot[1], geot[5])
            self._proj = basingrid._proj
            self._ncells = basin_cl.size
            self._ncells = basin_cl.size
            self._threshold = net._threshold
            self._thetaref = net._thetaref
            self._slp_np = net._slp_np
            self._ksn_np = net._ksn_np
            
            # Get only points inside the basin
            # The last receiver (ixc) will be outside of the basin
            # This can give index problems. We will throw away the last point in arrays
            basin_bool = basin.astype("bool").ravel()
            I = basin_bool[net._ix]
            self._ax = net._ax[I][:-1]
            self._zx = net._zx[I][:-1]
            self._dd = net._dd[I][:-1]
            self._dx = net._dx[I][:-1]
            self._chi = net._chi[I][:-1]
            self._slp = net._slp[I][:-1]
            self._ksn = net._ksn[I][:-1]
            self._r2slp = net._r2slp[I][:-1]
            self._r2ksn = net._r2ksn[I][:-1]
            
            # If Network object has less than 3 pixels, return an empty BNetwork
            if self._ax.size < 3:
                self._create_empty()
                return
            
            # Get new indices for the new grid
            ix = net._ix[I][:-1]
            ixc = net._ixc[I][:-1]
            # Givers (ix)
            row, col = net.indToCell(ix)
            x, y = net.cellToXY(row, col)
            newrow, newcol = self.xyToCell(x, y)
            self._ix = self.cellToInd(newrow, newcol)
            # Receivers (ixc)
            row, col = net.indToCell(ixc)
            x, y = net.cellToXY(row, col)
            newrow, newcol = self.xyToCell(x, y)
            self._ixc = self.cellToInd(newrow, newcol)
            self._heads = np.array([])
           
            if heads is not None:
                # Sort heads if "id" field is present
                if heads.shape[1] > 2:
                    pos = np.argsort(heads[:, 2])
                    heads = heads[pos]
                
                # Get heads inside the basin (taking into account nodata)
                heads = heads[basingrid.isInside(heads[:,0], heads[:,1], True)]
                
                # Snap heads to network heads
                heads = self.snapPoints(heads, "heads")                   

                # Get indexes
                row, col = self.xyToCell(heads[:,0], heads[:,1])
                idx = self.cellToInd(row, col)
                
                # Remove duplicate heads if any
                aux, pos = np.unique(idx, return_index=True)
                pos.sort()
                self._heads = idx[pos]
                
            if self._heads.size == 0:
                heads = self.streamPoi("heads", "IND")
                ixcix = np.zeros(self._ncells, np.int64)
                ixcix[self._ix] = np.arange(self._ix.size)
                pos = np.argsort(-self._zx[ixcix[heads]])
                self._heads = heads[pos][0:1] #Take only the first (highest) head as a list
            
    def _load(self, path):
        """
        Loads a BNetwork instance saved in the disk.
        
        Parameter:
        ==========
           Path to the saved BNetwork object (*.net file)
        """
        # Call to the parent Network._load() function
        super()._load(path)
        
        # Open again the *.dat file to get the heads
        fr = open(path, "r")
        for n in range(3):
            fr.readline()
        head_line = fr.readline()
        if  head_line[0] == "#":
            self._heads = np.array(head_line[1:-1].split(";")).astype(np.int64)
        # If the file hasn't got a four line in the header (with the heads) is not a BNetwork file
        else:
            raise NetworkError("The selected file is not a BNetwork objetct")
            
    def save(self, path):
        """
        Saves the Network instance to disk. It will be saved as a numpy array in text format with a header.
        The first three lines will have the information of the raster:
            Line1::   xsize; ysize; cx; cy; ULx; ULy; Tx; Tyy
            Line2::   thetaref; threshold; slp_np; ksn_np
            Line3::   String with the projection (WKT format)
            Line4::   head_1, head_2, ..., head_n
        xsize, ysize >> Dimensions of the raster
        cx, cy >> Cellsizes in X and Y
        Tx, Ty >> Rotation factors (for geotransformation matrix)
        ULx, ULy >> X and Y coordinates of the corner of the upper left pixel of the raster
        thetaref >>  m/n coeficient to calculate chi values in each channel cell
        threshold >> Number the cells to initiate a channel
        slp_np, ksn_np >> Number of points to calculate ksn and slope by regression. Window of {npoints * 2 + 1}
        
        Parameters:
        ===========
        path : *str*
          Path to save the network object with *.dat extension (it is not necessary to  give the extension)
        """
    
        # In case the extension is wrong or path has not extension
        path = os.path.splitext(path)[0] + ".dat"
        
        # Create header with properties
        params = [self._size[0], self._size[1], self._geot[1], self._geot[5], 
                  self._geot[0], self._geot[3], self._geot[2], self._geot[4]]
        header = ";".join([str(param) for param in params]) + "\n"  
        params = [self._thetaref, self._threshold, self._slp_np, self._ksn_np]
        header += ";".join([str(param) for param in params]) + "\n" 
        header += str(self._proj) + "\n"
        params = [str(head) for head in self._heads]
        header += ";".join([str(param) for param in params])
        
        # Create data array
        data_arr = np.array((self._ix, self._ixc, self._ax, self._dx, self._zx,
                             self._chi, self._slp, self._ksn, self._r2slp, 
                             self._r2ksn, self._dd)).T
        
        # Save the network instance as numpy.ndarray in text format
        np.savetxt(path, data_arr, delimiter=";", header=header, encoding="utf8", comments="#")
    
    def _create_empty(self):
        super()._create_empty()
        self._heads = np.array([0])
    
    def chiPlot(self, ax=None):
        """
        This function plot the Chi-elevation graphic for all the channels of the basin. 
        
        Parameters:
        ===========
        ax : matplotlib.Axe
          If is not defined, the function will create a new Figure and Axe
        """
        
        if not PLT:
            return

        if ax is None:
            fig, ax = plt.subplots()
            
        canales = self.getChannels(nchannels="ALL")
        main_ch = canales[0]
        
        for canal in canales:
            ax.plot(canal.getChi(), canal.getZ(), color="0.75", lw=0.75)
        
        ax.plot(main_ch.getChi(), main_ch.getZ(), color="k", lw=1)
        
        ax.set_xlim(xmin=0)
        ax.set_ylim(ymin=min(self._zx))
        
        ax.set_xlabel("Chi [m]")
        ax.set_ylabel("Elevation [m]")
        ax.set_title("Chi plot (m/n = {0:.2f})".format(self._thetaref))
   
    def chiSensitivityAnalysis(self, draw=False):
        pass

    def getChannel(self, id): 
        """
        Get the channel at specific index
        
        Parameters
        ----------
        id : int
          Channel index

        Returns
        -------
        Channel instance

        """
        channels = self.getChannels(nchannels="ALL")
        if id >= len(channels):
            raise NetworkError("Index ERROR")
        return channels[id]
        
    def getChannels(self, nchannels=None, min_length=0):
        """
        Get all channels in the basin. 

        Parameters
        ----------
        nchannels : int // None // "ALL", optional
            Number of channel that will be returned. If None (default) only channels corresponding to 
            Channel heads will be returned, otherwise, an specific number of channels will return. If nchannels 
            is greater than the Channel main heads, other channel will be returned for Network heads (sorted
            by elevation). To get all channels in the basin pass "ALL"
        min_length : float
            Minimun channel length. Channels of the BNetwork with lower lengths will be discarded
        Returns
        -------
        canales : list
            List of landspy.Channel objects

        """
        # Create auxiliary array to walk through the network cells
        aux_arr = np.zeros(self.getNCells(), "int")
        aux_arr.fill(-1)
        ixcix = np.zeros(self.getNCells(), "int")
        ixcix[self._ix] = np.arange(self._ix.size)
        # Empty list for returned Channel objects
        canales = []
        getflow = False
        
        # If nchannels is None, only Channel main heads will be returned
        if nchannels == None:
            heads = self._heads
       
        # Else, generate all heads for the basin (keeping main heads firts)
        else:
            # Get a list with all the heads and their elevations
            heads = self.streamPoi("heads", "IND")
            heads_z = self._zx[ixcix[heads]]
            
            # Sort heads by elevation (from highest to lowest)
            sorted_heads = heads[np.argsort(-heads_z)]
            
            # Merge sorted heads and main Channel heads
            all_heads = np.append(self._heads, sorted_heads)
            
            # Remove duplicates
            unique_val, unique_pos = np.unique(all_heads, return_index=True)
            
            # Sorted indexes (otherwise will be ordered according unique values)
            pos = np.sort(unique_pos)
            heads = all_heads[pos]
            
            # If nchannels == "ALL", return all the channels (also flow will be computed)
            if nchannels == "ALL":
                nchannels = heads.size
                getflow = True
        
        # Iterate heads and get Channels
        for idx, head in enumerate(heads[:nchannels]):
            flowto = -1
            chcells = [head]
            aux_arr[head] = True
            nextcell = self._ixc[ixcix[head]]
            aux_arr[nextcell] = True
            while ixcix[nextcell] != 0:
                chcells.append(nextcell)
                nextcell = self._ixc[ixcix[nextcell]]
                if aux_arr[nextcell]==-1:
                    aux_arr[nextcell] = idx
                else:
                    if getflow:
                        flowto = aux_arr[nextcell]
                    chcells.append(nextcell)
                    break
            
            row, col = self.indToCell(chcells)
            xi, yi = self.cellToXY(row, col)
            auxarr = np.zeros(self.getNCells()).astype("bool")
            auxarr[chcells] = True
            I = auxarr[self._ix]
            ai = self._ax[I]
            zi = self._zx[I]
            di = self._dx[I]
            chi = self._chi[I]
            dd = self._dd[I]
            slp = self._slp[I]
            ksn = self._ksn[I]
            r2slp = self._r2slp[I]
            r2ksn = self._r2ksn[I]
            chandata = np.array([chcells, ai, di, zi,  chi, slp, ksn, r2slp, r2ksn, dd]).T
            canal = Channel(self, chandata, self._thetaref, self._chi[-1], self._slp_np, self._ksn_np, str(idx), idx, flowto)
            if canal.getLength() >= min_length:
                canales.append(canal)
            
        return canales
     

class Channel(PRaster):

    def __init__(self, praster=None, chandata=None, thetaref=0.45, chi0=0, slp_np=5, ksn_np=5, name="", oid=-1, flowto=-1):
        """
        Class that defines a channel (cells from head to mouth)
        
        Parameters
        ----------
        praster : PRaster instance
            PRaster object to copy internal properties.
        chandata : numpy Array
            Array of 10 columns wiht channel data. These columns are:
                0 (ix)    >> Channel cells (ordered from head to mouth, in IND coordinates)
                1 (ax)    >> Upstream draining area (in pixel units)
                2 (dx)    >> Distance to mouth (distance from pixel to nearest outlet, NOT the channel mouth)
                3 (zx)    >> Elevation
                4 (chi)   >> Chi index
                5 (slp)   >> Slope of the pixel, calculated by regression with a moving window of {npoints * 2 + 1} cells.
                6 (ksn)   >> Ksn index of the pixel, calculated by regression with a moving window of {npoints * 2 + 1} cells.
                7 (r2slp) >> R2 Coeficient of the slope regression
                8 (r2ksn) >> R2 Coeficient of the ksn regression
                9 (dd)    >> Distance between channel cell and flowing cell (giver-receirver distance)
        thetaref : double, optional
            The m/n reference coeficient. The default is 0.45.
        chi0 : TYPE, optional
            Chi value of the mouth cell. The default is 0.
        slp_np : int, optional
            Number of points for slope calculation. Slope will be calculated within a moving window 
            of 2 * npoints + 1 . The default is 5.
        ksn_np : int, optional
            Number of points for ksn calculation. Ksn will be calculated within a moving window 
            of 2 * npoints + 1 . The default is 5.
            
        name : str, optional
            Name ("label") for the channel
            
        oid : int, optional
            Id of the channel
            
        flowto : int, optional
            Id of the channel where this channel flows. 

        Returns
        -------
        Channel object
        """
       
        # If praster is empty, create an empty Channel instancee
        if praster is None:
            self._create_empty()
            return
        # If praster parameter is a str, we load the Channel object from path
        elif type(praster)== str:
            self.load(praster)
            return
        
        # Initalize internal properties
        self._size = praster._size
        self._geot = praster._geot
        self._proj = praster._proj
        self._ix = chandata[:, 0].astype(np.int64)
        self._ax = chandata[:, 1]
        self._dx = chandata[:, 2]
        self._zx = chandata[:, 3]
        self._zx0 = np.copy(chandata[:, 3]) # Original elevations
        self._chi = chandata[:, 4]
        self._slp = chandata[:, 5]
        self._ksn = chandata[:, 6]
        self._r2slp = chandata[:, 7]
        self._r2ksn = chandata[:, 8]
        self._dd = chandata[:, 9]
        self._thetaref = thetaref
        self._chi0 = chi0
        self._slp_np = slp_np
        self._ksn_np = ksn_np
        self._kp = np.empty((0, 2), int)
        self._regressions = []
        self._name = name
        self._oid = oid
        self._flowto = flowto

    def _create_empty(self):
        """
        Creates an empty instance of the Channel class
        """
        # Set the empty PRaster properties
        super().__init__()
        
        # Set remaining properties
        self._ix =  np.array([0])
        self._ax = np.array([1])
        self._dx = np.array([0])
        self._zx = np.array([0])
        self._zx0 = np.array([0])
        self._chi = np.array([0])
        self._slp = np.array([0])
        self._ksn = np.array([0])
        self._r2slp = np.array([0])
        self._r2ksn = np.array([0])
        self._dd = np.array([1])
        self._thetaref = 0.45
        self._chi0 = 0
        self._slp_np = 0
        self._ksn_np = 0
        self._kp = np.empty((0, 2), int)
        self._regressions = []
        self._name = ''
        self._oid = -1
        self._flowto = -1
        
    def addKP(self, ind, tipo=0):
        """
        Adds a knickpoint at a given position

        Parameters
        ----------
        ind : int
            Knickpoint position (index of the knickpoint within the channel)
        tipo : int, optional
            Integer representing knickpoint type.

        Returns
        -------
        None.

        """
        self._kp = np.append(self._kp, [[ind, tipo]], axis=0)
        
    def removeKP(self, ind):
        """
        Remove knickoint from channel

        Parameters
        ----------
        ind : int
            Knickpoint position (index of the knickpoint within the channel)

        Returns
        -------
        None.

        """
        pos = np.where(self._kp[:,0] == ind)[0]
        self._kp = np.delete(self._kp, pos, axis=0) 
    
    def addRegression(self, p1, p2):
        """
        Add regression in chi-elevation space.
        A regression is defined by a tuple of 5 elements: [id, pos1, pos2, poly, r2ksn]
          id: index of the regression. Is the middle possition (pos2-pos1)/2
          pos1, pos2: positions of the first and last points of the regression
          poly: 1nd order polynomial with the regression >> [a, b] (y = ax + b)
          r2ksn: R2 of the regression

        Parameters
        ----------
        p1 : int
            Position of the start of the regression (index within the channel)
        p2 : int
            Position of the end of the regression (index within the channel)

        Returns
        -------
        None.

        """
        # If p2 is greater than p1, change values
        if p2 < p1:
            p2, p1 = p1, p2
        
        # Return if p1 or p2 are equal or not valid indexes
        if p1 == p2 or p1 < 0 or p2 >= self._ix.size:
            return
        
        # Get values of Chi-Elevation for regression
        chi = self._chi[p1:p2]
        zi = self._zx[p1:p2]
        
        # Calculate gradient by linear regression ()
        poli, SCR = np.polyfit(chi, zi, deg = 1, full = True)[:2]
        
        # Calculate R2 
        if zi.size * zi.var() == 0:
            R2 = 1 # Puntos colineares
        else:
            R2 = float(1 - SCR/(zi.size * zi.var()))
            
        # Get Id (mid-point of the regression)
        idx = int(p1 + (p2 - p1)/2)
        
        # Append regression and returns regression id
        self._regressions.append((idx, p1, p2, poli, R2))
        return idx
    
    def getRegression(self, ind):
        """
        Get regression from Channel based on a position. In case of multiple regressions passing through the position,
        returns the closest to the mid-point of the regression. 
        """
        if len(self._regressions) == 0:
            return
        
        regs = []
        for reg in self._regressions:
            if reg[1] <= ind <= reg[2]:
                regs.append(reg)
        min_dist = 99999999
        out_reg = None
        for reg in regs:
            diff = abs(reg[0] - ind)
            if diff < min_dist:
                out_reg = reg
        
        return out_reg
    
    def removeRegression(self, idx):
        """
        Remove regression from channel

        Parameters
        ----------
        idx : int
            Regression index

        Returns
        -------
        None.

        """
        # Get the regression to remove (equal index of idx)
        remove_reg = None
        for reg in self._regressions:
            if reg[0] == idx:
                remove_reg = reg
        
        if remove_reg:
            self._regressions.remove(remove_reg)
            return True
        else:
            return False
                  
    def save(self, path):
        """
        Saves the Channel instance to disk. It will be saved as a numpy array in text format with a header.
        The first three lines will have the information of the raster:
            Line1::   name; oid; flowto
            Line2::   xsize; ysize; cx; cy; ULx; ULy; Tx; Ty
            Line3::   thetaref; threshold; slp_np; ksn_np
            Line4::   String with the projection (WKT format)
        xsize, ysize >> Dimensions of the raster
        cx, cy >> Cellsizes in X and Y
        Tx, Ty >> Rotation factors (for geotransformation matrix)
        ULx, ULy >> X and Y coordinates of the corner of the upper left pixel of the raster
        thetaref >>  m/n coeficient to calculate chi values in each channel cell
        threshold >> Number the cells to initiate a channel
        slp_np, ksn_np >> Number of points to calculate ksn and slope by regression. Window of {npoints * 2 + 1}
        
        Parameters:
        ===========
        path : *str*
          Path to save the network object with *.dat extension (it is not necessary to  give the extension)
        """
    
        # In case the extension is wrong or path has not extension
        path = os.path.splitext(path)[0] + ".dat"
        
        # Create header with properties
        # Line 1: name; oid; flowto
        header = ";".join([self._name, str(self._oid), str(self._flowto)]) + "\n"
        
        # Line2: xsize; ysize; cx; cy; ULx; ULy; Tx; Ty
        header += ";".join([str(param) for param in [self._size[0], self._size[1], self._geot[1], self._geot[5], 
                  self._geot[0], self._geot[3], self._geot[2], self._geot[4]]]) + "\n"
        
        # Line 3: thetaref; chi0; slp_np; ksn_np
        header += ";".join([str(param) for param in [self._thetaref, self._chi0, self._slp_np, self._ksn_np]]) + "\n"  
        
        # Line 4: proj
        header += str(self._proj) + "\n" 
        
        # Line 5: knickpoints
        kp_line = ""
        for kp in self._kp:
            kp_line += str(kp[0]) + ";" + str(kp[1]) + ";"
        header += kp_line[:-1] + "\n"
        
        # Line 6: regressions
        reg_line = ""
        for reg in self._regressions:
            reg_line += str(reg[0]) + ";" + str(reg[1]) + ";"
        header += reg_line[:-1]   
            
        # Create data array
        data_arr = np.array((self._ix, self._ax, self._dx, self._zx,
                             self._chi, self._slp, self._ksn, self._r2slp, 
                             self._r2ksn, self._dd, self._zx0)).T
        
        # Save the Channel instance as numpy.ndarray in text format
        np.savetxt(path, data_arr, delimiter=";", header=header, encoding="utf8", comments="#")
        
    def load(self, path):
        """
        Loads a Network instance saved in the disk.
        
        Parameter:
        ==========
           Path to the Channel object
        """
        # Open the file as normal text file to get its properties
        fr = open(path, "r")
        
        # For all the first 6 lines, first and last characters will be "#" and "\n" respectively 
        
        # Line 1: name; oid; flowto
        linea = fr.readline()[1:-1]
        data = linea.split(";")
        self._name = data[0]
        self._oid = int(data[1])
        self._flowto = int(data[2])
        
        # Line 2: xsize; ysize; cx; cy; ULx; ULy; Tx; Ty
        linea = fr.readline()[1:-1]
        data = linea.split(";")
        self._size = (int(data[0]), int(data[1]))
        self._geot = (float(data[4]), float(data[2]), float(data[6]), 
                      float(data[5]), float(data[7]), float(data[3]))       
       
        # Line3: thetaref; chi0; slp_np; ksn_np
        linea = fr.readline()[1:-1]
        data = linea.split(";")
        self._thetaref = float(data[0])
        self._chi0 = float(data[1])
        self._slp_np = int(data[2])
        self._ksn_np = int(data[3])
        
        # Line4: Proj (wkt)
        linea = fr.readline()[1:-1]
        self._proj = linea
        
        # Line5: Knickpoints
        self._kp = np.empty((0, 2), int)
        kp_line = fr.readline()[1:-1]
            
        # Line6: Regressions
        self._regressions = []
        reg_line = fr.readline()[1:-1]
        
        fr.close()
        
        # Load array data
        data_arr = np.loadtxt(path, dtype=float, comments='#', delimiter=";", encoding="utf8")
        # Fix to avoid errors in networks with only one cell...
        if data_arr.ndim < 2:
            data_arr = data_arr.reshape((1, data_arr.size))
        self._ix = data_arr[:, 0].astype(np.int64)
        self._ax = data_arr[:, 1]
        self._dx = data_arr[:, 2]
        self._zx = data_arr[:, 3]
        self._chi = data_arr[:, 4]
        self._slp = data_arr[:, 5]
        self._ksn = data_arr[:, 6]
        self._r2slp = data_arr[:, 7]
        self._r2ksn = data_arr[:, 8]
        self._dd = data_arr[:, 9]
        self._zx0 = data_arr[:, 10]
        
        
        # Add knickpoints and regressions
        if kp_line:
            data = kp_line.split(";")
            for n in range(0, len(data), 2):
                ind = int(data[n])
                tipo = int(data[n+1])
                self.addKP(ind, tipo)
            
        if reg_line:
            data = reg_line.split(";")
            for n in range(0, len(data), 2):
                p1 = int(data[n])
                p2 = int(data[n+1])
                self.addRegression(p1, p2)
        
    def setName(self, name):
        """
        Sets the name (label) of the channel
        """
        self._name = name
        
    def getName(self):
        """
        Returns the name (label) of the channel
        """
        return self._name
    
    def getOid(self):
        """
        Gets the id of the channel. Necessary to compute flow. 
        """
        return self._oid
        
    def getFlow(self):
        """
        Gets the id of the channels where this channel flows. Necessary to compute flow. 
        """
        return self._flowto


    def getLength(self):
        """
        Returns channel lenght
        """
        return self._dx[0] - self._dx[-1]
    
    def getXY(self, head=True):
        """
        Returns channel coordiates (numpy.array with two columns, x and y)
        """
        row, col = self.indToCell(self._ix)
        x, y = self.cellToXY(row, col)
        return np.array((x, y)).T
    
    def getZ(self, head=True, relative=False):
        """
        Returns channel elevations
        
        Parameters
        ----------
        head : boolean, optional
            Returns elevations from head to mouth (True) or mouth to head (False). The default is True.
        relative : boolean, optional
            Returns relative elevations (True) or real elevations (True). The default is False.

        Returns
        -------
        1D numpy array
        """
        zi = np.copy(self._zx)
        
        if relative:
            zi = zi - zi[-1]
        if not head:
            zi = zi[::-1]
        
        return zi

    def getChi(self, head=True, relative=False):
        """
        Returns channel chi values
        
        Parameters
        ----------
        head : boolean, optional
            Returns values from head to mouth (True) or mouth to head (False). The default is True.
        relative : boolean, optional
            Returns relative chi value (True) or real chi values (True). The default is False.

        Returns
        -------
        1D numpy array
        """
        chi = np.copy(self._chi)
        
        if relative:
            chi = chi - chi[-1]
        if not head:
            chi = chi[::-1]
        
        return chi
      
    def getA(self, head=True, cells=True):
        """
        Returns channel area values (in cell units)
        
        Parameters
        ----------
        head : boolean, optional
            Returns values from head to mouth (True) or mouth to head (False). The default is True.
        cells : boolean, optional
            Returns area values in cells (True), or in length units (False)

        Returns
        -------
        1D numpy array
        """
        ai = np.copy(self._ax)
        
        if not cells:
            ai = ai * self.getCellSize()[0] * self.getCellSize()[1] * -1
        if not head:
            ai = ai[::-1]

        return ai    
    
    def calculateGradients(self, npoints, kind='slp'):
        """
        This function calculates gradients (slope or ksn) for all the cells. 
        Gradients of each cell are calculated by linear regression using a number
        of points (npoints) up and downstream.
        
        Parameters:
        ===========
        npoints : *int*
          Window to analyze slopes. Slopes are calculated by linear regression using a window
          of npoints * 2 + 1 pixel (using the central pixel)
          
        kind : *str* {'slp', 'ksn'} 
          Calculates the gradients for slope (distance-elevation) or kwn (chi-elevation)
        """
        if npoints < 2:
            return

        # Get arrays depending on type
        y_arr = self._zx
        if kind == 'ksn':
            x_arr = self._chi
        else:
            x_arr = self._dx
        
        for n in range(self._ix.size):
            low = n - npoints
            high = n + npoints
        
            if low < 0:
                low = 0
        
            xi = x_arr[low:high + 1]
            yi = y_arr[low:high + 1]
            poli, SCR = np.polyfit(xi, yi, deg = 1, full = True)[:2]
            
            if yi.size * yi.var() == 0:
                R2 = 0.0
            else:
                R2 = float(1 - SCR/(yi.size * yi.var()))
        
            g = poli[0]
            r2 = abs(R2)
        
            if abs(g) < 0.001:
               g = 0.001
               
               
            if kind == 'ksn':
                self._ksn[n] = g
                self._r2ksn[n] = r2
                self._ksn_np = npoints
            else:
                self._slp[n] = g
                self._r2slp[n] = r2
                self._slp_np = npoints
        
    def getSlope(self, head=True):
        """
        Returns channel slope values
        
        Parameters
        ----------
        head : boolean, optional
            Returns values from head to mouth (True) or mouth to head (False). The default is True.
        
        Returns
        -------
        1D numpy array
        """
        slp = np.copy(self._slp)
        if not head:
            slp = slp[::-1]
        return slp
            
    def getD(self, tohead=True, head=True):
        """
        Returns channel distante values
        
        Parameters
        ----------
        tohead : boolean, optional
            Distance are computed from cell to head (True) or from cell to mouth (False)
        head : boolean, optional
            Returns distance values from head to mouth (True) or mouth to head (False). The default is True.
        Returns
        -------
        1D numpy array
        """
        di = np.copy(self._dx)
        if tohead:
            di = self._dx[0] - di
        else:
            di = di - self._dx[-1]
        
        if not head:
            di = di[::-1]
        
        return di

    def getKsn(self, head=True):
        """
        Returns ksn values
        
        Parameters
        ----------
        head : boolean, optional
            Returns values from head to mouth (True) or mouth to head (False). The default is True.
        
        Returns
        -------
        1D numpy array
        """
        ksn = np.copy(self._ksn)
        if not head:
            ksn = ksn[::-1]
        return ksn
    
    def smoothChannel(self, winsize=0, recalculate_gradients=True):
        """
        This function smooth channel elevations win a moving window
        Parameters:
        ===========
        winsize : *double*
          Size of the moving windows (meters or profile units). The window size will be transformed to number of pixels. 
          
        recalculate_gradients : *boolean*
          Flag to indicate if recalulate gradients (slope and ksn) with the new elevation values
        """
        # Remove peaks and flat segments
        for n in range(self._zx.size - 1):
            if self._zx[n + 1] >= self._zx[n]:
                self._zx[n + 1] = self._zx[n] - 0.001

        # Smooth elevations if window distance > 0
        if winsize > 0:
            n_cells = int(int((winsize / self.getCellSize()[0]) + 0.5) / 2)
            z_copy = np.copy(self._zx)
            for ind in range(self._zx.size):
                low = ind - n_cells
                high = ind + n_cells + 1
                if low < 0:
                    low = 0
        
                elevations = z_copy[low:high]
                self._zx[ind] = np.mean(elevations)
            
            if recalculate_gradients:
                self.calculateGradients(self._ksn_np, kind='ksn')
                self.calculateGradients(self._slp_np, kind='slp')

class NetworkError(Exception):
    pass