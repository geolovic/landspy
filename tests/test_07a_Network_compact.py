"""Network traversals on sparse layouts, independent of generated test data."""

import unittest
from unittest.mock import patch

import numpy as np

from landspy import DEM, Flow, Network
from landspy._network import compact_nodes, accumulate_downstream


class CompactNetworkTest(unittest.TestCase):
    def test_reverse_topology_with_shared_receivers_and_multiple_outlets(self):
        ix = np.array([100, 200, 300, 400], dtype=np.int64)
        ixc = np.array([300, 300, 500, 600], dtype=np.int64)
        nodes, givers, receivers = compact_nodes(ix, ixc)
        result = accumulate_downstream(givers, receivers,
                                       np.array([1., 2., 3., 4.]), nodes.size)
        np.testing.assert_array_equal(result, [4., 5., 3., 4.])

    def test_large_cell_ids_do_not_allocate_full_raster(self):
        net = Network()
        net._size = (100000, 100000)
        net._ix = np.array([2**32 + 1, 2**32 + 2], dtype=np.int64)
        net._ixc = np.array([2**32 + 2, 2**32 + 3], dtype=np.int64)
        net._dd = np.array([2., 3.])
        net._ax = np.array([1., 4.])
        net._zx = np.array([20., 10.])
        net._dx = np.array([5., 3.])
        zeros = np.zeros

        def bounded_zeros(shape, *args, **kwargs):
            self.assertLessEqual(int(np.prod(shape)), 32)
            return zeros(shape, *args, **kwargs)

        with patch('numpy.zeros', side_effect=bounded_zeros):
            net.calculateChi(thetaref=.5, a0=2.)
            net.calculateGradients(5, 'slp')
            net.calculateGradients(5, 'ksn')
            for kind, expected in [('heads', [2**32 + 1]),
                                   ('outlets', [2**32 + 3]),
                                   ('confluences', [])]:
                np.testing.assert_array_equal(net.streamPoi(kind, 'IND'), expected)
        np.testing.assert_array_equal(net._chi, [7., 3.])
        np.testing.assert_allclose(net._slp, [0., 5.], rtol=1e-14)
        np.testing.assert_allclose(net._ksn, [0., 2.5], rtol=1e-14)
        np.testing.assert_array_equal(net._r2slp, [0., 1.])
        np.testing.assert_array_equal(net._r2ksn, [0., 1.])

    def test_distances_preserve_coordinate_rounding(self):
        dem = DEM()
        dem.setArray(np.array([[10, 9, 8], [7, 6, 5], [4, 3, 2]], np.int16))
        dem._geot = (1e12, .13, 0., 2e12, 0., -.27)
        flow = Flow(dem)
        net = Network(flow, threshold=1)
        expected_dd = np.zeros(net._ix.size)
        expected_dx = np.zeros(net.getNCells())
        for n in range(net._ix.size - 1, -1, -1):
            gx, gy = net.cellToXY(*net.indToCell(net._ix[n]))
            rx, ry = net.cellToXY(*net.indToCell(net._ixc[n]))
            expected_dd[n] = np.sqrt((gx-rx)**2 + (gy-ry)**2)
            expected_dx[net._ix[n]] = expected_dx[net._ixc[n]] + expected_dd[n]
        np.testing.assert_array_equal(net._dd, expected_dd)
        np.testing.assert_array_equal(net._dx, expected_dx[net._ix])

    def test_threshold_above_accumulation_has_no_network_nodes(self):
        dem = DEM()
        dem.setArray(np.array([[3, 2], [2, 1]], np.int16))
        net = Network(Flow(dem), threshold=5, gradients=True)
        for name in ['_ix', '_ixc', '_dd', '_dx', '_chi', '_slp', '_ksn']:
            self.assertEqual(getattr(net, name).size, 0)
        for kind in ['heads', 'outlets', 'confluences']:
            self.assertEqual(net.streamPoi(kind, 'XY').shape, (0, 2))


if __name__ == '__main__':
    unittest.main()
