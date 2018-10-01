# -*- coding: utf-8 -*-

# network.py
# Jose Vicente Perez Pena
# Dpto. Geodinamica-Universidad de Granada
# 18071 Granada, Spain
# vperez@ugr.es // geolovic@gmail.com
#
# MIT License (see LICENSE file)
# Version: 1.0
# December 26, 2017
#
# Last modified 31 May, 2018

import numpy as np
import gdal
from scipy.sparse import csc_matrix
from . import Grid, PRaster, Flow, DEM
  

class Network(PRaster):

    def __init__(self, dem, flow, threshold):
        """
        Class to manipulate cells from a Network, which is defined by applying
        a threshold to a flow accumulation raster derived from a topological 
        sorted Flow object.
        
        Parameters:
        -----------
        dem : *topopy.DEM* object
          Digital Elevation Model instance
        flow : *topopy.Flow* object
          Flow direccion instance
        threshold : *int*
          Number the cells to initiate a channel
        """
        # Set PRaster properties
        self._size = flow.get_size()
        self._dims = flow.get_dims()
        self._geot = flow.get_geotransform()
        self._cellsize = flow.get_cellsize()
        self._proj = flow.get_projection()
        self._ncells = flow.get_ncells()
      
        # Get sort Nodes for channel cells
        fac = flow.get_flow_accumulation(nodata=False, asgrid=False)
        w = fac > threshold
        w = w.ravel()
        I   = w[flow._ix]
        self._ix  = flow._ix[I]
        self._ixc = flow._ixc[I]
        self._ax = fac.ravel()[self._ix] * self._cellsize**2 # Area in map units
        self._zx = dem.read_array().ravel()[self._ix]
        
        ## TODO
        self._slp = np.zeros(self._ix.shape)
        self._chi = np.zeros(self._ix.shape)
    
    def calculate_di(self):
        di = np.zeros(self._ncells, np.float)
        
        nix = self._ix.size
        for n in np.arange(nix)[::-1]:
            grow, gcol = self.ind_2_cell(self._ix[n])
            gx, gy = self.cell_2_xy(grow, gcol)
            rrow, rcol = self.ind_2_cell(self._ixc[n])
            rx, ry = self.cell_2_xy(rrow, rcol)
            d_gr = np.sqrt((gx - rx)**2 + (gy - ry)**2) 
            di[self._ix[n]] = di[self._ixc[n]] + d_gr
    
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
        
        # Get array with channel cells
        ix  = self._ix
        ixc = self._ixc
        w = np.zeros(self._ncells, dtype=np.bool)
        w[self._ix] = True
        w[self._ixc] = True
        
        # Build a sparse array with giver-receivers cells
        aux_vals = np.ones(ix.shape, dtype=np.int8)
        sp_arr = csc_matrix((aux_vals, (ix, ixc)), shape=(self._ncells, self._ncells))
        
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

    def get_streams(self, asgrid=True):
        """
        This function extract a drainage network by using a determined area threshold

        Parameters:
        ===========
        asgrid : *bool*
          Indicates if the network is returned as topopy.Grid (True) or as a numpy.array
        """
        # Get grid channel cells
        w = np.zeros(self._ncells, dtype=np.int8)
        w[self._chcells] = 1
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
        # Get heads
        head_ind = self.get_stream_poi("heads")
        head_ind = self.cell_2_ind(head_ind[0], head_ind[1])
        
        # Get confluences
        conf_ind = self.get_stream_poi("confluences")
        conf_ind = self.cell_2_ind(conf_ind[0], conf_ind[1])
                
        # Merge all heads and confluences
        all_ind = np.append(head_ind, conf_ind)
        del conf_ind, head_ind # Clean up
        
        # We created a zeros arrays and put in confuences and heads their id
        # Those id will be consecutive numbers starting in one
        seg_arr = np.zeros(self._ncells, dtype=np.int32)
        for n, inds in enumerate(all_ind):
            seg_arr[inds] = n+1
        
        # Move throught giver list. If receiver is 0, give the same id that giver.
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

    def get_stream_order(self, kind="strahler", asgrid=True):
        """
        This function extract streams orderded by strahler or shreeve
    
        Parameters:
        ===========
        kind : *str* {'strahler', 'shreeve'}
        asgrid : *bool*
          Indicates if the network is returned as topopy.Grid (True) or as a numpy.array
        """
        if kind not in ['strahler', 'shreeve']:
            return
        
        # Get grid channel cells
        str_ord = np.zeros(self._ncells, dtype=np.int8)
        str_ord[self._chcells] = 1
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
    
    def get_main_channels(self, heads=None):
        # Aceptar salida a vector
        pass

    def calculate_chi(self, thetaref=0.45, a0=1.0):
        tmpix = self._ix[::-1]
        tmpixc = self._ixc[::-1]
        tmpai = np.ones(self._ncells, dtype=np.int)
        tmpai[self._ix] = self._ax
        tmpai = self._ax[::-1]**thetaref
        chi = np.zeros(self._ncells, dtype=np.float)
        try:
            for n in range(tmpix.size):
                gr, gc = self.ind_2_cell(tmpix[n])
                rr, rc = self.ind_2_cell(tmpixc[n])
                gx, gy = self.cell_2_xy(gr, gc)
                rx, ry = self.cell_2_xy(rr, rc)
                dx = np.sqrt((rx-gx)**2 + (ry-gy)**2)
                chi[tmpixc[n]] = chi[tmpix[n]] + (a0 * dx / tmpai[tmpixc[n]])
        except:
            print(n)
            print(chi.size)
            print(tmpix[n])
            print(tmpixc[n])
                
    
    def calculate_slope(self, dist=None):
        pass
    
    def snap_points(self, row, col, kind="heads"):
        """
        Snap points to network points of interest {heads, confluences, or outlets}
        
        Parameters:
        -----------
        row : row indexes (number, list, or numpy.ndarray) of input points
        col : column indexes (number, list, or numpy.ndarray) of input points
        
        Return:
        =======
        Tuple with (row, col) indices of snap points as np.ndarrays
        """
        n_points = len(row)
        snap_row = []
        snap_col = []
        # Get stream poi
        poi = self.get_stream_poi(kind = kind)
        # Calculate distances
        for n in range(n_points):
            dist = np.sqrt((row[n]-poi[0])**2 + (col[n]-poi[1])**2)
            min_pos = np.argmin(dist)
            snap_row.append(poi[0][min_pos])
            snap_col.append(poi[1][min_pos])
            
        return np.array(snap_row), np.array(snap_col)
    
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
    
#    def get_stream_poi(self, threshold, kind="heads"):
#        """
#        This function finds points of interest of the drainage network. These points of interest
#        can be 'heads', 'confluences' or 'outlets'.
#        
#        Parameters:
#        -----------
#        threshold : *int* 
#          Area threshold to initiate a channels (in cells)
#        kind : *str* {'heads', 'confluences', 'outlets'}
#          Kind of point of interest to return. 
#          
#        Returns:
#        -----------
#        (row, col) : *tuple*
#          Tuple of numpy nd arrays with the location of the points of interest
#          
#        References:
#        -----------
#        The algoritms to extract the point of interest have been adapted to Python 
#        from Topotoolbox matlab codes developed by Wolfgang Schwanghart (version of 17. 
#        August, 2017). These smart algoritms use sparse arrays with giver-receiver indexes, to 
#        derive point of interest in a really efficient way. Cite:
#                
#        Schwanghart, W., Scherler, D., 2014. Short Communication: TopoToolbox 2 - 
#        MATLAB-based software for topographic analysis and modeling in Earth 
#        surface sciences. Earth Surf. Dyn. 2, 1–7. https://doi.org/10.5194/esurf-2-1-2014
#        """
#        # Check input parameters
#        if kind not in ['heads', 'confluences', 'outlets']:
#            return np.array([]), np.array([])
#        
#        # Get the Flow Accumulation and select cells that meet the threholds
#        fac = self.flow_accumulation(nodata = False).read_array()
#        w = fac > threshold
#        del fac
#        
#        # Build a sparse array with giver-receivers cells
#        w = w.ravel()
#        I   = w[self._ix]
#        ix  = self._ix[I]
#        ixc = self._ixc[I]
#        aux_vals = np.ones(ix.shape, dtype=np.int8)
#    
#        sp_arr = csc_matrix((aux_vals, (ix, ixc)), shape=(self._ncells, self._ncells))
#        del I, ix, ixc, aux_vals # Clean up (Don't want to leave stuff in memory)
#        
#        if kind == 'heads':
#            # Heads will be channel cells marked only as givers (ix) but not as receivers (ixc) 
#            sum_arr = np.asarray(np.sum(sp_arr, 0)).ravel()
#            out_pos = (sum_arr == 0) & w
#        elif kind == 'confluences':
#            # Confluences will be channel cells with two or givers
#            sum_arr = np.asarray(np.sum(sp_arr, 0)).ravel()
#            out_pos = sum_arr > 1
#        elif kind == 'outlets':
#            # Outlets will be channel cells marked only as receivers (ix) but not as givers (ixc) 
#            sum_arr = np.asarray(np.sum(sp_arr, 1)).ravel()
#            out_pos = (sum_arr == 0) & w  
#            
#        out_pos = out_pos.reshape(self._dims)
#        row, col = np.where(out_pos)
#        
#        return row, col    
#
#    def get_stream_order(self, threshold, kind="strahler", asgrid=True):
#        """
#        This function extract stream orders by using a determined area threshold
#    
#        Parameters:
#        ===========
#        threshold : *int* 
#          Area threshold to initiate a channels (in cells)
#        kind : *str* {'strahler', 'shreeve'}
#        asgrid : *bool*
#          Indicates if the network is returned as topopy.Grid (True) or as a numpy.array
#        """
#        if kind not in ['strahler', 'shreeve']:
#            return
#        
#        fac = self.get_flow_accumulation(nodata=False, asgrid=False)
#        w = fac > threshold
#        w = w.ravel()
#        I   = w[self._ix]
#        ix  = self._ix[I]
#        ixc = self._ixc[I]
#    
#        str_ord = np.copy(w).astype(np.int8)
#        visited = np.zeros(self._ncells, dtype=np.int8)
#    
#        if kind == 'strahler':
#            for n in range(len(ix)):
#                if (str_ord[ixc[n]] == str_ord[ix[n]]) & visited[ixc[n]]:
#                    str_ord[ixc[n]] = str_ord[ixc[n]] + 1
#                else:
#                    str_ord[ixc[n]] = max(str_ord[ix[n]], str_ord[ixc[n]])
#                    visited[ixc[n]] = True
#        elif kind == 'shreeve':
#            for n in range(len(ix)):
#                if visited[ixc[n]]:
#                    str_ord[ixc[n]] = str_ord[ixc[n]] + str_ord[ix[n]]
#                else:
#                    str_ord[ixc[n]] = max(str_ord[ix[n]], str_ord[ixc[n]])
#                    visited[ixc[n]] = True
#        str_ord = str_ord.reshape(self._dims)
#        
#        if asgrid:
#            return self._create_output_grid(str_ord, nodata_value=0)
#        else:
#            return str_ord
#
#    def get_network(self, threshold, asgrid=True):
#        """
#        This function extract a drainage network by using a determined area threshold
#
#        Parameters:
#        ===========
#        threshold : *int* 
#          Area threshold to initiate a channels (in cells)
#        asgrid : *bool*
#          Indicates if the network is returned as topopy.Grid (True) or as a numpy.array
#        """
#        # Get the Flow Accumulation and select cells that meet the threholds
#        fac = self.get_flow_accumulation(nodata = False, asgrid=False)
#        w = fac > threshold
#        w = w.astype(np.int8)
#        if asgrid:
#            return self._create_output_grid(w, 0)
#        else:
#            return w
#        
   
#
#    def _get_presills(self, flats, sills, dem):
#        """
#        This functions extracts the presill pixel locations (i.e. pixelimmediately 
#        upstream to sill pixels)- Adapted from TopoToolbox matlab codes.
#        """
#        zvals = dem.read_array()
#        flats_arr = flats.read_array().astype("bool")
#        dims = zvals.shape
#        row, col = sills.find()
#        
#        rowadd = np.array([-1, -1, 0, 1, 1,  1,  0, -1])
#        coladd = np.array([ 0,  1, 1, 1, 0, -1, -1, -1])
#        ps_rows = np.array([], dtype=np.int32)
#        ps_cols = np.array([], dtype=np.int32)
#        
#        for n in range(8):
#            rowp = row + rowadd[n]
#            colp = col + coladd[n]
#            # Avoid neighbors outside array (remove cells and their neighbors)
#            valid_rc = (rowp >= 0) & (colp >= 0) & (rowp < dims[0]) & (colp < dims[1])
#            rowp = rowp[valid_rc]
#            colp = colp[valid_rc]
#            # Discard cells (row-col pairs) that do not fullfill both conditions
#            cond01 = zvals[row[valid_rc], col[valid_rc]] == zvals[rowp, colp]
#            cond02 = flats_arr[rowp, colp]
#            valid_pix = np.logical_and(cond01, cond02)
#            ps_rows = np.append(ps_rows, rowp[valid_pix])
#            ps_cols = np.append(ps_cols, colp[valid_pix])
#      
#        ps_pos = list(zip(ps_rows, ps_cols))
#        return ps_pos
#       
#    def _get_topodiff(self, topodiff, flats):
#        """
#        This function calculate an auxiliar topography to sort the flats areas
#        """
#        tweight = 2
#        carvemin = 0.1                       
#        struct = np.ones((3, 3), dtype="int8")
#        lbl_arr, nlbl = ndimage.label(flats.read_array(), structure=struct)
#        lbls = np.arange(1, nlbl + 1)
#        
#        for lbl in lbls:
#            topodiff[lbl_arr == lbl] = (topodiff[lbl_arr==lbl].max() - topodiff[lbl_arr == lbl])**tweight + carvemin
#                    
#        return topodiff
#    
#    def _get_weights(self, flats, topodiff, ps_pos):
#        """
#        This function calculate weights in the flats areas by doing a cost-distance analysis.
#        It uses presill positions as seed locations, and an auxiliar topography as friction
#        surface.
#        """
#        flats_arr = flats.read_array().astype("bool")
#        flats_arr = np.invert(flats_arr)
#        topodiff[flats_arr] = 99999
#        if len(ps_pos) > 0:
#            lg = graph.MCP_Geometric(topodiff)
#            topodiff = lg.find_costs(starts=ps_pos)[0] + 1
#        topodiff[flats_arr] = -99999
#        
#        return topodiff
#    
#    def _sort_pixels(self, dem, weights):
#        """
#        Sort the cells of a DEM in descending order. It uses weights for the flats areas
#        """
#        # Sort the flat areas
#        rdem = dem.read_array().ravel()
#        rweights = weights.ravel()
#        ix_flats = np.argsort(-rweights, kind='quicksort')
#        
#        # Sort the rest of the pixels from the DEM
#        ndx = np.arange(self._ncells, dtype=np.int)
#        ndx = ndx[ix_flats]
#        ix = ndx[np.argsort(-rdem[ndx], kind='mergesort')]
#        
#        return ix
#    
#    def _get_receivers(self, ix, dem):
#        
#        dem = dem.read_array()
#        rdem = dem.ravel()
#        
#        pp = np.zeros(self._dims, dtype=np.int32)
#        IX = np.arange(self._ncells, dtype=np.int32)
#        pp = pp.ravel()
#        pp[ix] = IX
#        pp = pp.reshape(self._dims)
#                
#        # Get cardinal neighbors
#        footprint= np.array([[0, 1, 0],
#                             [1, 1, 1],
#                             [0, 1, 0]], dtype=np.int)
#        IXC1 = ndimage.morphology.grey_dilation(pp, footprint=footprint)
#        xxx1 = np.copy(IXC1)
#        IX = IXC1.ravel()[ix]
#        IXC1 = ix[IX]
#        G1   = (rdem[ix]-rdem[IXC1])/(self._cellsize)
#        
#        # Get diagonal neighbors
#        footprint= np.array([[1, 0, 1],
#                             [0, 1, 0],
#                             [1, 0, 1]], dtype=np.int)
#        IXC2 = ndimage.morphology.grey_dilation(pp, footprint=footprint)
#        xxx2 = np.copy(IXC2)
#        IX = IXC2.ravel()[ix]
#        IXC2 = ix[IX]
#        G2   = (rdem[ix]-rdem[IXC2])/(self._cellsize * np.sqrt(2))
#        
#        del IX
#        
#        # Get the steepest one
#        I  = (G1<=G2) & (xxx2.ravel()[ix]>xxx1.ravel()[ix])
#        ixc = IXC1
#        ixc[I] = IXC2[I]
#        
#        return ixc
class FlowError(Exception):
    pass