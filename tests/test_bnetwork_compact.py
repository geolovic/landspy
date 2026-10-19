"""Basin extraction and channel compatibility against the original traversal.

Last modified: 2026-10-19 11:56:48 CEST (Europe/Madrid).
"""
import unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
from landspy import Flow, Network, BNetwork, Basin, DEM, Grid
from landspy.network import Channel, NetworkError

DATA = Path(__file__).resolve().parent / "data" / "in"

def reference_channels(self, nchannels=None, min_length=0):
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

class CompactBNetworkTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.flow = Flow(str(DATA / "small25_fd.tif"))
        cls.net = Network(cls.flow)
        outlets = np.loadtxt(DATA / "small25_bnet_outlets.txt", delimiter=";")
        cls.grid = cls.flow.drainageBasins(cls.net.snapPoints(outlets))
        cls.heads = np.loadtxt(DATA / "small25_bnet_heads.txt", delimiter=";")
        cls.dem = DEM(str(DATA / "small25.tif"))

    def assert_channel(self, actual, expected):
        for name in expected.__dict__:
            np.testing.assert_equal(getattr(actual, name), getattr(expected, name),
                                    err_msg=name)

    def test_constructor_and_channels(self):
        for bid in (1, 2, 3):
            mask = self.grid.readArray() == bid
            rows, cols = np.where(mask)
            r0, c0 = rows.min(), cols.min()
            selected = np.flatnonzero(mask.ravel()[self.net._ix])[:-1]
            for grid in (self.grid, Basin(self.dem, self.grid, bid)):
                for heads in (None, self.heads):
                    with self.subTest(bid=bid, kind=type(grid).__name__, heads=heads is not None):
                        b = BNetwork(self.net, grid, heads, bid)
                        self.assertEqual(b.getDims(), (rows.max()-r0+1, cols.max()-c0+1))
                        for name in ("_ax", "_zx", "_dd", "_dx", "_chi", "_slp", "_ksn", "_r2slp", "_r2ksn"):
                            np.testing.assert_array_equal(getattr(b, name), getattr(self.net, name)[selected])
                            self.assertFalse(np.shares_memory(getattr(b, name), getattr(self.net, name)))
                        for name in ("_ix", "_ixc"):
                            rr, cc = self.net.indToCell(getattr(self.net, name)[selected])
                            np.testing.assert_array_equal(getattr(b, name), b.cellToInd(rr-r0, cc-c0))
                        for count in (None, 0, 1, 4, "ALL"):
                            for length in (0, 1000, 100000, np.nan):
                                expected = reference_channels(b, count, length)
                                actual = b.getChannels(count, length)
                                self.assertEqual(len(actual), len(expected))
                                for a, e in zip(actual, expected):
                                    self.assert_channel(a, e)

    def test_single_channel_only_materializes_requested_object(self):
        b = BNetwork(self.net, self.grid)
        expected = reference_channels(b, "ALL")
        for idx in (0, 1, len(expected)-1, -1, -len(expected)):
            with self.subTest(idx=idx), patch("landspy.network.Channel", wraps=Channel) as constructor:
                self.assert_channel(b.getChannel(idx), expected[idx])
                self.assertEqual(constructor.call_count, 1)
        with self.assertRaises(NetworkError):
            b.getChannel(len(expected))
        with self.assertRaises(IndexError):
            b.getChannel(-len(expected)-1)
        with self.assertRaises(TypeError):
            b.getChannel(0.5)

    def test_bounding_box_across_blocks(self):
        net = Network()
        net._size = (2000, 1200)
        net._ix = np.arange(520, 540, dtype=np.int64) * 2000 + 42
        net._ixc = net._ix + 2000
        for name in ("_ax", "_zx", "_dd", "_dx", "_chi", "_slp", "_ksn", "_r2slp", "_r2ksn"):
            setattr(net, name, np.arange(20, 0, -1, dtype=np.float64))
        grid = Grid()
        grid.setArray(np.zeros(net.getDims(), dtype=np.uint8))
        grid._geot = net._geot
        grid.readArray()[520:541, 42] = 7
        basin = BNetwork(net, grid, bid=7)
        self.assertEqual(basin.getDims(), (21, 1))
        np.testing.assert_array_equal(basin._ix, np.arange(19))
        np.testing.assert_array_equal(basin._ixc, np.arange(1, 20))

    def test_inputs_unchanged_and_empty_basin(self):
        values = self.grid.readArray().copy()
        BNetwork(self.net, self.grid, self.heads, 2)
        np.testing.assert_array_equal(values, self.grid.readArray())
        b = BNetwork(self.net, self.grid, bid=999)
        self.assertEqual(b.getSize(), (1, 1))


if __name__ == "__main__":
    unittest.main()
