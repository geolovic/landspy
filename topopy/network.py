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
# Last modified October 23th, 2018

import numpy as np
import os
import ogr, osr
from scipy.sparse import csc_matrix
from . import Grid, PRaster

# This import statement avoid issues with matplotlib in Mac when using Python not as a Framework
# If matplotlib is not imported, BNetwork.chi_plot() will not work.
try: 
    import matplotlib.pyplot as plt
    PLT = True
except Exception as e:
    PLT = False

class Network(PRaster):

    def __init__(self, flow=None, threshold=0, thetaref=0.45, npoints=5, gradients=True):
        """
        Class to manipulate cells from a Network, which is defined by applying
        a threshold to a flow accumulation raster derived from a topological 
        sorted Flow object.
        
        Parameters:
        -----------
        flow : *topopy.Flow* object
          Flow direccion instance
        threshold : *int*
          Number the cells to initiate a channel
        thetaref : *float*
          m/n coeficient to calculate chi values in each channel cell    
        npoints : *int*
          Number of points to calculate slope and ksn in each cell. Slope and ksn values
          are calculated with a moving window of (npoints * 2 + 1) cells.
        """
        # If flow is a str, load it
        if type(flow)== str:
            self._load(flow)
            return
        
        # Set PRaster properties
        self._size = flow.get_size()
        self._dims = flow.get_dims()
        self._geot = flow.get_geotransform()
        self._cellsize = flow.get_cellsize()
        self._proj = flow.get_projection()
        self._ncells = flow.get_ncells()
      
        # Get a threshold if not specified (Default 0.25% of the total number of cells)
        if threshold == 0:
            threshold = self._ncells * 0.0025
        self._threshold = int(threshold)
        
        # Get sort Nodes for channel cells and elevations
        fac = flow.get_flow_accumulation(nodata=False, asgrid=False)
        w = fac >= threshold
        w = w.ravel()
        I   = w[flow._ix]
        self._ix  = flow._ix[I]
        self._ixc = flow._ixc[I]
        
        # Get Area and Elevations for channel cells
        self._ax = fac.ravel()[self._ix] * self._cellsize[0] * self._cellsize[1] * -1 # Area in map units
        self._zx = flow._zx[I]
        
        # Get distances to mouth (self._dx) and giver-receiver distances (self._dd)
        di = np.zeros(self._ncells)
        self._dd = np.zeros(self._ix.shape) # Giver-Receiver distance
        for n in np.arange(self._ix.size)[::-1]:
            grow, gcol = self.ind_2_cell(self._ix[n])
            rrow, rcol = self.ind_2_cell(self._ixc[n])
            gx, gy = self.cell_2_xy(grow, gcol)
            rx, ry = self.cell_2_xy(rrow, rcol)
            d_gr = np.sqrt((gx - rx)**2 + (gy - ry)**2)
            self._dd[n] = d_gr
            di[self._ix[n]] = di[self._ixc[n]] + d_gr
        self._dx = di[self._ix]
        
        # Get chi values using the input thetaref
        self._thetaref = thetaref
        self.calculate_chi(thetaref)
        
        # Calculate slopes and ksn
        if gradients:
            self.calculate_gradients(npoints, 'slp')
            self.calculate_gradients(npoints, 'ksn')
        else:
            self._slp = np.zeros(self._ix.size)
            self._r2slp = np.zeros(self._ix.size)
            self._slp_npoints = 0
            self._ksn = np.zeros(self._ix.size)
            self._r2ksn = np.zeros(self._ix.size)
            self._ksn_npoints = 0

    def save(self, path):
        """
        Saves the Network instance to disk. It will be saved as a numpy array in text format.
        The first three rows will have the information of the raster
        
        Parameters:
        ===========
        path : *str*
          Path to save the network object, with *.net extension
        """
    
        path = os.path.splitext(path)[0]
            
        # Create *.net file with properties
        netfile = open(path + ".net", "w")
        params = [self._size, self._cellsize, self._ncells, self._ksn_npoints, 
                  self._slp_npoints, self._thetaref, self._threshold] 
        linea = ";".join([str(param) for param in params]) + "\n"
        netfile.write(linea)
        linea = ";".join([str(elem) for elem in self._geot]) + "\n"
        netfile.write(linea)
        netfile.write(str(self._proj))
        netfile.close()
        
        # Create data array
        data_arr = np.array((self._ix, self._ixc, self._ax, self._dx, self._zx,
                             self._chi, self._slp, self._ksn, self._r2slp, 
                             self._r2ksn, self._dd)).T
        
        # Save the network instance as numpy.ndarray in text format
        np.save(path + ".npy", data_arr)
    
    def _load(self, path):
        """
        Loads a Network instance saved in the disk.
        
        Parameter:
        ==========
           Path to the saved network object
        """
        # Get properties from properties file *.net
        netfile = path
        fr = open(netfile, "r")
        linea = fr.readline()[:-1]
        data = linea.split(";")
        self._size = (int(data[0].split(",")[0][1:]), int(data[0].split(",")[1][:-1]))
        self._dims = (self._size[1], self._size[0])
        self._cellsize = (float(data[1].split(",")[0][1:]), float(data[1].split(",")[1][:-1]))
        self._ncells = int(data[2])
        self._ksn_npoints = int(data[3])
        self._slp_npoints = int(data[4])
        self._thetaref = float(data[5])
        self._threshold = int(data[6])
        linea = fr.readline()[:-1]
        self._geot = tuple([float(n) for n in linea.split(";")])
        linea = fr.readline()
        self._proj = linea
        fr.close()
        
        # Load array data from the auxiliar array file *.npy
        arrfile = os.path.splitext(netfile)[0] + ".npy"
        data_arr = np.load(arrfile)
        self._ix = data_arr[:, 0].astype(np.int)
        self._ixc = data_arr[:, 1].astype(np.int)
        self._ax = data_arr[:, 2]
        self._dx = data_arr[:, 3]
        self._zx = data_arr[:, 4]
        self._chi = data_arr[:, 5]
        self._slp = data_arr[:, 6]
        self._ksn = data_arr[:, 7]
        self._r2slp = data_arr[:, 8]
        self._r2ksn = data_arr[:, 9]
        self._dd = data_arr[:, 10]
         
    def calculate_chi(self, thetaref=0.45, a0=1.0):
        """
        Function that calculate chi_values for channel cells
        
        Parameters:
        -----------
        thetaref : *float*
          m/n coeficient to calculate chi
        a0 : *float*
          Reference area to avoid dimensionality (usually don't need to be changed)
        """
        chi = np.zeros(self._ncells)
        for n in np.arange(self._ix.size)[::-1]:
            chi[self._ix[n]] = chi[self._ixc[n]] + (a0 * self._dd[n]/self._ax[n]**thetaref)            
        self._chi = chi[self._ix]
        self._thetaref = thetaref
      
    def calculate_gradients(self, npoints, kind='slp'):
        """
        This function calculates gradients (slope or ksn) for all channel cells. 
        Gradients of each cell are calculated by linear regression using a number
        of points (npoints) up and downstream.
        
        Parameters:
        ===========
        npoints : *int*
          Window to analyze slopes. Slopes are calculated by linear regression using a window
          of npoints * 2 + 1 pixel (using the central pixel)
          
        kind : *str* {'slp', 'ksn'}
        """
        if kind not in ['slp', 'ksn']:
            kind = 'slp'
        winlen = npoints * 2 + 1
        
        # Get arrays depending on type
        if kind == 'slp':
            x_arr = self._dx
            y_arr = self._zx
        elif kind == 'ksn':
            x_arr = self._chi
            y_arr = self._zx
            
        # Get ixcix auxiliar array
        ixcix = np.zeros(self._ncells, np.int)
        ixcix[self._ix] = np.arange(self._ix.size)
        
        # Get heads array and confluences dictionary
        heads = self.get_stream_poi("heads", "IND")
        confs = {conf:[] for conf in self.get_stream_poi("confluences", "IND")}
        
        # Sort heads by elevation
        elev = self._zx[ixcix[heads]]
        spos = np.argsort(-elev)
        heads = heads[spos]
        
        # Prepare auxiliary arrays
        gi = np.zeros(self._ncells)
        r2 = np.zeros(self._ncells)
        
        # Taking sequentally all the heads and compute downstream flow
        for head in heads:
            processing = True
            head_cell = head
            mid_cell = self._ixc[ixcix[head_cell]]
            mouth_cell = self._ixc[ixcix[mid_cell]]
        
            if ixcix[mid_cell] == 0:
                # Channel type 2 (mid_cell is an outlet)
                continue
            elif ixcix[mouth_cell]== 0:
                # Channel type 3 (mouth_cell is an outlet)
                processing = False
            # Check si cabecera es una confluencia de solo dos celdas
            elif gi[mid_cell] != 0:
                processing = False
            
            # Obtenemos datos de elevacion y distancias
            win = [mouth_cell, mid_cell, head_cell]
            xi = x_arr[ixcix[win]]
            yi = y_arr[ixcix[win]]
           
            # Calculamos pendiente de celda central por regresion
            poli, SCR = np.polyfit(xi, yi, deg = 1, full = True)[:2]
            g = poli[0]
            if g == 0: 
                g = 0.000001
        
            # Calculamos gradient y R2
            if yi.size * yi.var() == 0:
                R2 = 1 # Puntos colineares
            else:
                R2 = float(1 - SCR/(yi.size * yi.var()))
        
            gi[mid_cell] = g
            r2[mid_cell] = R2
                    
            while processing:
                # Verificamos si estamos al final (en un outlet)  
                if not ixcix[mouth_cell]==0:
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
                
                if len(win) <= 3:
                    processing = False
                    
                # Obtenemos datos de elevacion y distancias para calcular pendientes
                xi = x_arr[ixcix[win]]
                yi = y_arr[ixcix[win]]
                poli, SCR = np.polyfit(xi, yi, deg = 1, full = True)[:2]
                g = poli[0]
                if g == 0: 
                    g = 0.000001
                
                if yi.size * yi.var() == 0:
                    R2 = 1
                else:
                    R2 = float(1 - SCR/(yi.size * yi.var()))
                            
                # Comprobamos si celda central ha sido procesada
                if gi[mid_cell] == 0:
                    gi[mid_cell] = g
                    r2[mid_cell] = R2
                else:
                    # Si ha sido procesada, es una confluencia 
                    processing = False
                    if len(confs[mid_cell]) == 0:
                        confs[mid_cell].append(gi[mid_cell])
                        confs[mid_cell].append(g)
                    else:
                        confs[mid_cell].append(g)
                        
        # Calculamos valores medios en confluencias                
        for cell in confs.keys():
            if len(confs[cell]) > 0:
                gi[cell] = np.mean(np.array(confs[cell]))
        
        # Llenamos array del objeto Network
        if kind == 'slp':
            self._slp = gi[self._ix]
            self._r2slp = r2[self._ix]
            self._slp_npoints = npoints
        elif kind == 'ksn':
            self._ksn = gi[self._ix]
            self._r2ksn = r2[self._ix]
            self._ksn_npoints = npoints
        
    def get_stream_poi(self, kind="heads", coords="CELL"):
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
            # Check input parameters
            if kind not in ['heads', 'confluences', 'outlets']:
                kind = 'heads'
            if coords not in ['CELL', 'XY', 'IND']:
                coords = 'CELL'

            # Get grid channel cells
            w = np.zeros(self._ncells, dtype=np.bool)
            w[self._ix] = True
            w[self._ixc] = True
            
            # Build a sparse array with giver-receivers cells
            aux_vals = np.ones(self._ix.shape, dtype=np.int8)
            sp_arr = csc_matrix((aux_vals, (self._ix, self._ixc)), shape=(self._ncells, self._ncells))
            
            # Get stream POI according the selected type
            if kind == 'heads':
                # Heads will be channel cells marked only as givers (ix) but not as receivers (ixc) 
                sum_arr = np.asarray(np.sum(sp_arr, 0)).ravel()
                out_pos = (sum_arr == 0) & w
            elif kind == 'confluences':
                # Confluences will be channel cells with two or givers
                sum_arr = np.asarray(np.sum(sp_arr, 0)).ravel()
                out_pos = sum_arr > 1
            elif kind == 'outlets':
                # Outlets will be channel cells marked only as receivers (ixc) but not as givers (ix) 
                sum_arr = np.asarray(np.sum(sp_arr, 1)).ravel()
                out_pos = np.logical_and((sum_arr == 0), w)  
                
            out_pos = out_pos.reshape(self._dims)
            row, col = np.where(out_pos)
            
            if coords=="CELL":
                return np.array((row, col)).T
            elif coords=="XY":
                xi, yi = self.cell_2_xy(row, col)
                return np.array((xi, yi)).T
            elif coords=="IND":
                return self.cell_2_ind(row, col)

    def snap_points(self, input_points, kind="channel"):
        """
        Snap input points to channel cells or to stream POI
        
        Parameters:
        ===========
        array : *numpy.ndarray*
          Numpy 2-D ndarray, which first two columns are x and y coordinates [x, y, ...]
        kind : *str* {'channel', 'heads', 'confluences', 'outlets'}  
            Kind of point to snap input points
        
        Returns:
        ===========
        numpy.ndarray
          Numpy ndarray with two columns [xi, yi] with the snap points
        """
        if kind not in  ['channel', 'heads', 'confluences', 'outlets']:
            kind = 'channel'
        
        # Extract a numpy array with the coordinate to snap the points
        if kind in ['heads', 'confluences', 'outlets']:
            poi = self.get_stream_poi(kind, "XY")     
        else:
            row, col = self.ind_2_cell(self._ix)
            x, y = self.cell_2_xy(row, col)
            poi = np.array((x, y)).T
        
        # Get array reshaped for the calculation
        xi = input_points[:, 0].reshape((input_points.shape[0], 1))
        yi = input_points[:, 1].reshape((input_points.shape[0], 1))
        xci = poi[:, 0].reshape((1, poi.shape[0]))
        yci = poi[:, 1].reshape((1, poi.shape[0]))
        
        # Calculate distances and get minimum
        di = np.sqrt((xi - xci)**2 + (yi - yci)**2 )
        pos = np.argmin(di, axis=1)
        
        return poi[pos]
    
    def export_to_points(self, path):
        """
        Export channel points to a semicolon-delimited text file
        
        path : str
          Path for the output text file
        """
        cab = "id;x;y;z;distance;area;chi;slope;ksn;r2_slope;r2_ksn"
        row, col = self.ind_2_cell(self._ix)
        x, y = self.cell_2_xy(row, col)
        
        out_arr = np.array((self._ix, x, y, self._zx, self._dx, self._ax, self._chi, 
                            self._slp, self._ksn, self._r2slp, self._r2ksn)).T
        np.savetxt(path, out_arr, delimiter=";", header=cab, comments="")
    
    def get_streams(self, asgrid=True):
        """
        This function outputs a grid representation of the Network object

        Parameters:
        ===========
        asgrid : *bool*
          Indicates if the network is returned as topopy.Grid (True) or as a numpy.array
        """
        # Get grid channel cells
        w = np.zeros(self._ncells, dtype=np.int8)
        w[self._ix] = 1
        w[self._ixc] = 1
        w = w.reshape(self._dims)
        # Return grid
        if asgrid:
            return self._create_output_grid(w, 0)
        else:
            return w
    
    def get_stream_segments(self, asgrid=True):
        """
        This function extract a drainage network by using a determined area threshold
        and output the numerated stream segments.

        Parameters:
        ===========
        asgrid : *bool*
          Indicates if the network is returned as topopy.Grid (True) or as a numpy.array
        """
        # Get heads and confluences and merge them
        head_ind = self.get_stream_poi("heads", "IND")
        conf_ind = self.get_stream_poi("confluences", "IND")
        all_ind = np.append(head_ind, conf_ind)
        del conf_ind, head_ind # Clean up
        
        # We created a zeros arrays and put in confluences and heads their id
        # Those id will be consecutive numbers starting in one
        seg_arr = np.zeros(self._ncells, dtype=np.int32)
        for n, inds in enumerate(all_ind):
            seg_arr[inds] = n+1
        
        # Move throught channel list. If receiver is 0, give receiver the same id that giver.
        # If a receiver is not 0, that means that we are in a confluence. 
        for n in range(len(self._ix)):
            if seg_arr[self._ixc[n]] == 0:
                seg_arr[self._ixc[n]] = seg_arr[self._ix[n]]
        
        # Reshape and output
        seg_arr = seg_arr.reshape(self._dims)
        if asgrid:
            return self._create_output_grid(seg_arr, 0)
        else:
            return seg_arr
        
    def get_stream_orders(self, kind="strahler", asgrid=True):
        """
        This function extract streams orderded by strahler or shreeve. Cell values
        will have a value acording with the order of the segment they belong
    
        Parameters:
        ===========
        kind : *str* {'strahler', 'shreeve'}
        asgrid : *bool*
          Indicates if the selfwork is returned as topopy.Grid (True) or as a numpy.array
        """
        if kind not in ['strahler', 'shreeve']:
            return
        
        # Get grid channel cells
        str_ord = np.zeros(self._ncells, dtype=np.int8)
        str_ord[self._ix] = 1
        str_ord[self._ixc] = 1
        visited = np.zeros(self._ncells, dtype=np.int8)
    
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
        str_ord = str_ord.reshape(self._dims)
        
        if asgrid:
            return self._create_output_grid(str_ord, nodata_value=0)
        else:
            return str_ord

    def export_to_shp(self, path, con=False):
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
        layer = dataset.CreateLayer("rivers", sp, ogr.wkbLineString)
        layer.CreateField(ogr.FieldDefn("segid", ogr.OFTInteger))
        layer.CreateField(ogr.FieldDefn("order", ogr.OFTInteger))
        layer.CreateField(ogr.FieldDefn("flowto", ogr.OFTInteger))
        
        # Get channel segments and orders
        ch_seg = self.get_stream_segments(False).ravel()
        ch_ord = self.get_stream_orders(asgrid=False).ravel()
        ch_seg = ch_seg[self._ix]
        ch_ord = ch_ord[self._ix]
        
        # Get ixcix auxiliar array
        ixcix = np.zeros(self._ncells, np.int)
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
            order = ch_ord[ixcix[first]]
            if ixcix[last] == 0:
                flowto = idx
            else:
                flowto = ch_seg[ixcix[last]]
            
            # Add feature
            feat = ogr.Feature(layer.GetLayerDefn())
            feat.SetField("segid", int(idx))
            feat.SetField("order", int(order))
            feat.SetField("flowto", int(flowto))
            row, col = self.ind_2_cell(ch_ix)
            xi, yi = self.cell_2_xy(row, col)
            
            geom = ogr.Geometry(ogr.wkbLineString)
            
            for n in range(xi.size):
                geom.AddPoint(xi[n], yi[n])
                
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
        layer = dataset.CreateLayer("rivers", sp, ogr.wkbLineString)
        layer.CreateField(ogr.FieldDefn("order", ogr.OFTInteger))
        
        # Get ixcix auxiliar array
        ixcix = np.zeros(self._ncells, np.int)
        ixcix[self._ix] = np.arange(self._ix.size)
        
        # Get heads, confluences and orders
        heads = self.get_stream_poi("heads", "IND")
        confs = self.get_stream_poi("confluences", "IND")
        ch_ord = self.get_stream_orders(asgrid=False).ravel()
        
        # Get confluences where strahler index increases
        strahler_confs = []
        for conf in confs:
            conf_order = ch_ord[conf]
            givers = self._ix[np.where(self._ixc==conf)]
            giv_orders = ch_ord[givers]
            if giv_orders.max() < conf_order:
                strahler_confs.append(conf)
                
        # Append strahler confluences to heads
        heads = np.append(heads, np.array(strahler_confs))    
        
        # Iterate heads
        for head in heads:
            cell = head
            river_data = [cell]
            processing = True
            while processing:
                next_cell = self._ixc[ixcix[cell]]
                river_data.append(next_cell)
                if next_cell in confs:
                    if ch_ord[next_cell] > ch_ord[cell]:
                        processing = False
                    else:
                        river_data.append(next_cell)
                        cell = next_cell
                elif ixcix[next_cell] == 0:
                    processing = False
                else:
                    river_data.append(next_cell)
                    cell = next_cell    
                            
            # Add feature
            feat = ogr.Feature(layer.GetLayerDefn())
            row, col = self.ind_2_cell(river_data)
            xi, yi = self.cell_2_xy(row, col)
            geom = ogr.Geometry(ogr.wkbLineString)
            for n in range(xi.size):
                geom.AddPoint(xi[n], yi[n])
                
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
        grid.copy_layout(self)
        grid._nodata = nodata_value
        grid._array = array
        grid._tipo = str(array.dtype)
        return grid


class BNetwork(Network):
    
    def __init__(self, net, basin=None, heads=None, bid=1):
        """
        Class to manipulate cells from a drainage network from a single basin network. 
        This class inhereits all methods and properties of the Network class plus some 
        new methods to manipulate the drainage network and the trunk channel. 
        
        Parameters:
        -----------
        net : *topopy.Network* | *str*
          Network instance or path to a previously saved BNetwork file
        basin : *numpy.ndarray* | *topopy.Grid*
          Numpy array or topopy Grid representing the drainage basin. If array or Grid have more than
          one basin, set the basinid properly. 
        heads : *list* or *numpy.ndarray*
          List with [x, y] coordinates for basin the main head or 2 column numpy.ndarray with head(s) 
          coordinate(s). If more than one, the first one is considered the main head (trunk channel). 
          Head(s) will snap to Network heads.
        basinid : *int*
          Id value that identifies the basin cells in case that basin will have more than one basin.        
        """
        
        # If net is a str, load it
        if type(net)== str:
            self._load(net)
            return
        
        #try:
        # If basin is not a numpy array is a topopy.Grid
        if not isinstance(basin, np.ndarray):
            basin = basin.read_array()
        
        # Set basin cells to 1
        basin = np.where(basin==bid, 1, 0)
                
        # Get basin extent
        c1 = basin.max(axis=0).argmax() - 1
        r1 = basin.max(axis=1).argmax() - 1
        c2 = basin.shape[1] - np.fliplr(basin).max(axis=0).argmax() + 1 
        r2 = basin.shape[0] - np.flipud(basin).max(axis=1).argmax() + 1
        
        # Check boundary conditions
        if c1 < 0:
            c1 = 0
        if r1 < 0:
            r1 = 0
        if c2 >= basin.shape[1]:
            c2 = basin.shape[1] - 1
        if r2 >= basin.shape[0]:
            r2 = basin.shape[0] - 1
            
        basin_cl = basin[r1:r2, c1:c2] # clipped basin
        
        # Create the new grid
        self._size = (basin_cl.shape[1], basin_cl.shape[0])
        self._dims = (basin_cl.shape[0], basin_cl.shape[1])
        geot = net._geot
        ULx = geot[0] + geot[1] * c1
        ULy = geot[3] + geot[5] * r1
        self._geot = (ULx, geot[1], 0.0, ULy, 0.0, geot[5])
        self._cellsize = (geot[1], geot[5])
        self._proj = net._proj
        self._ncells = basin_cl.size
        self._threshold = net._threshold
        self._thetaref = net._thetaref
        self._slp_npoints = net._slp_npoints
        self._ksn_npoints = net._ksn_npoints
        
        # Get only points inside the basin
        basin = basin.astype(np.bool).ravel()
        I = basin[net._ix]
        self._ax = net._ax[I]
        self._zx = net._zx[I]
        self._dd = net._dd[I]
        self._dx = net._dx[I]
        self._chi = net._chi[I]
        self._slp = net._slp[I]
        self._ksn = net._ksn[I]
        self._r2slp = net._r2slp[I]
        self._r2ksn = net._r2ksn[I]
        
        # Get new indices for the new grid
        ix = net._ix[I]
        ixc = net._ixc[I]
        rowix, colix = net.ind_2_cell(ix)
        rowixc, colixc = net.ind_2_cell(ixc)
        nrowix = rowix - r1
        ncolix = colix - c1
        nrowixc = rowixc - r1
        ncolixc = colixc - c1
        self._ix = self.cell_2_ind(nrowix, ncolix)
        self._ixc = self.cell_2_ind(nrowixc, ncolixc)
        
        # Set the basin heads and sort them by elevation
        ixcix = np.zeros(self._ncells, np.int)
        ixcix[self._ix] = np.arange(self._ix.size)
        bheads = self.get_stream_poi(kind="heads", coords="IND")
        elev = self._zx[ixcix[bheads]]
        spos = np.argsort(-elev)
        bheads = bheads[spos].tolist()
        
        if heads is None:
            self._nheads = 1
            self._heads = np.array(bheads)
        else:
            # Snap user heads to Network heads
            heads = self.snap_points(heads, kind="heads")
            self._nheads = heads.shape[0]
            row, col = self.xy_2_cell(heads[:, 0], heads[:, 1])
            heads = self.cell_2_ind(row, col)
            heads_in = []
            for head in heads:
                if head in bheads:
                    heads_in.append(head)
                    bheads.remove(head)
                
            self._heads = np.append(np.array(heads_in), np.array(bheads))

    def _load(self, path):
        """
        Loads a BNetwork instance saved in the disk.
        
        Parameter:
        ==========
           Path to the saved BNetwork object (*.net file)
        """
#        try:
        # Call to the parent Network._load() function
        super()._load(path)
        
        # Open again the *.net file to get the heads
        fr = open(path, "r")
        lineas = fr.readlines()
        self._heads = np.array(lineas[3].split(";")).astype(np.int)
        self._nheads = self._heads[0]
        self._heads = self._heads[1:]
        fr.close()
#        except:
#            raise NetworkError("Error loading the BNetwork object")
            
    def save(self, path):
        """
        Saves the BNetwork instance to disk. It will be saved as text file (*.net) plus
        a numpy binary array with the same name (*.npy) with the cell data.
        
        Parameters:
        ===========
        path : *str*
          Path to save the network object, with *.net extension
        """
        # Call to the parent Network.save() function
        super().save(path)
        
        # Open again the *.net file to append the heads
        netfile = open(os.path.splitext(path)[0] + ".net", "a")
        heads = ";".join(self._heads.astype(np.str))
        netfile.write("\n" + str(self._nheads) + ";" + heads)
        netfile.close()

    def chi_plot(self, ax=None):
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
            
        main_ch = self.get_main_channel()
        ax.plot(self._chi, self._zx, color="0.75", ls="None", marker=".", ms=1)
        ax.plot(main_ch[:, 5], main_ch[:, 2], ls="-", c="0.3", lw=1)
        ax.set_xlim(xmin=0)
        ax.set_ylim(ymin=0)
        
        ax.set_xlabel("Chi [m]")
        ax.set_ylabel("Elevation [m]")
        ax.set_title("Chi plot (m/n = {0:.2f})".format(self._thetaref))
   
    def chi_analysis(self, draw=False):
        pass

### ^^^^ UP HERE ALL FUNCTIONS TESTED ^^^^
    def get_main_channel(self):
        chcells = [self._heads[0]]
        ixcix = np.zeros(self._ncells, np.int)
        ixcix[self._ix] = np.arange(self._ix.size)
        nextcell = self._ixc[ixcix[self._heads[0]]]
        
        while ixcix[nextcell] != 0:
            chcells.append(nextcell)
            nextcell = self._ixc[ixcix[nextcell]]
                     
        row, col = self.ind_2_cell(chcells)
        xi, yi = self.cell_2_xy(row, col)
        auxarr = np.zeros(self._ncells).astype(np.bool)
        auxarr[chcells] = True
        I = auxarr[self._ix]
        ai = self._ax[I]
        zi = self._zx[I]
        di = self._dx[I]
        chi = self._chi[I]
        slp = self._slp[I]
        ksn = self._ksn[I]
        r2slp = self._r2slp[I]
        r2ksn = self._r2ksn[I]
#        length = di[0] - di[-1]
#        di -= di[-1]
#        di = length - di

        outarr = np.array([xi, yi, zi, di, ai, chi, slp, ksn, r2slp, r2ksn]).T
        return outarr

    def get_channels(self, nchannels=0):
        aux_arr = np.zeros(self._ncells, np.bool)
        ixcix = np.zeros(self._ncells, np.int)
        ixcix[self._ix] = np.arange(self._ix.size)
        canales = []
        
        if nchannels == 0:
            heads = self._heads
        else:
            heads = self._heads[:nchannels]
        
        for head in heads:
            chcells = [head]
            aux_arr[head] = True
            nextcell = self._ixc[ixcix[head]]
            aux_arr[nextcell] = True
            while ixcix[nextcell] != 0:
                chcells.append(nextcell)
                nextcell = self._ixc[ixcix[nextcell]]
                if aux_arr[nextcell]==False:
                    aux_arr[nextcell] = True
                else:
                    chcells.append(nextcell)
                    break
            
            row, col = self.ind_2_cell(chcells)
            xi, yi = self.cell_2_xy(row, col)
            auxarr = np.zeros(self._ncells).astype(np.bool)
            auxarr[chcells] = True
            I = auxarr[self._ix]
            ai = self._ax[I]
            zi = self._zx[I]
            di = self._dx[I]
            chi = self._chi[I]
            slp = self._slp[I]
            ksn = self._ksn[I]
            r2slp = self._r2slp[I]
            r2ksn = self._r2ksn[I]
            chandata = np.array([xi, yi, zi, di, ai, chi, slp, ksn, r2slp, r2ksn]).T
            canales.append(chandata)
        
        return canales
     
class NetworkError(Exception):
    pass