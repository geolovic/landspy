"""Numerical regressions for the compiled Network construction traversals."""

# Last modified: 2026-10-16 11:19:50 CEST (Europe/Madrid).

import unittest

import numpy as np

from landspy import Flow, Grid, Network
from landspy._network_kernels import receiver_positions, downstream_sums


def branching_flow(dtype='uint32'):
    flow = Flow()
    flow._size = (7, 6)
    flow._geot = (359123.7105, 27.006999386539967, 0.,
                  599062.1513, 0., -23.006514936323367)
    # Two outlets, merging tributaries, diagonal and cardinal edges. Givers
    # are topologically ordered but deliberately not sorted by cell index.
    flow._ix = np.array([31, 0, 8, 9, 16, 17, 24, 25, 5], dtype=dtype)
    flow._ixc = np.array([24, 8, 16, 16, 24, 25, 32, 32, 6], dtype=dtype)
    flow._zx = np.arange(20., 11., -1)
    flow._nodata_pos = np.setdiff1d(np.arange(42), np.union1d(flow._ix, flow._ixc))
    return flow


def reference_accumulation(flow, initial=None):
    values = (np.ones(flow.getNCells(), dtype=np.uint32) if initial is None
              else np.array(initial, dtype=np.float64).ravel().copy())
    for giver, receiver in zip(flow._ix, flow._ixc):
        values[receiver] += values[giver]
    values[flow._nodata_pos] = 0
    return values


def reference_metrics(flow, threshold, theta, a0=1.):
    accumulation = reference_accumulation(flow)
    selected = accumulation[flow._ix] >= threshold
    ix, ixc = flow._ix[selected], flow._ixc[selected]
    area = accumulation[ix]
    dd = np.zeros(ix.size)
    distance = np.zeros(flow.getNCells())
    chi = np.zeros(flow.getNCells())
    for n in range(ix.size - 1, -1, -1):
        gx, gy = flow.cellToXY(*flow.indToCell(ix[n]))
        rx, ry = flow.cellToXY(*flow.indToCell(ixc[n]))
        dd[n] = np.sqrt((gx - rx)**2 + (gy - ry)**2)
        distance[ix[n]] = distance[ixc[n]] + dd[n]
        chi[ix[n]] = chi[ixc[n]] + a0 * dd[n] / area[n]**theta
    return dict(_ix=ix, _ixc=ixc, _ax=area, _zx=flow._zx[selected],
                _dd=dd, _dx=distance[ix], _chi=chi[ix])


class NetworkConstructionTest(unittest.TestCase):
    def test_constructor_matches_raster_reference(self):
        for dtype in ('uint32', 'int32', 'int64'):
            flow = branching_flow(dtype)
            for threshold in (0, 1, 2.5, 3, 5, 100):
                for theta in (0., .45, 1.):
                    with self.subTest(dtype=dtype, threshold=threshold, theta=theta):
                        net = Network(flow, threshold=threshold, thetaref=theta)
                        effective = threshold or flow.getNCells() * .0025
                        expected = reference_metrics(flow, effective, theta)
                        for name, values in expected.items():
                            np.testing.assert_allclose(getattr(net, name), values,
                                                       rtol=1e-14, atol=1e-14)
                        np.testing.assert_array_equal(net._slp, np.zeros(net._ix.size))
                        np.testing.assert_array_equal(net._ksn, np.zeros(net._ix.size))

    def test_recalculate_chi_with_reference_area(self):
        flow = branching_flow()
        net = Network(flow, threshold=1)
        for theta, a0 in ((.3, 2.), (.7, .5), (0., 1.)):
            net.calculateChi(theta, a0)
            expected = reference_metrics(flow, 1, theta, a0)
            np.testing.assert_allclose(net._chi, expected['_chi'], rtol=1e-14)
            self.assertEqual(net._thetaref, theta)

    def test_accumulation_weighted_and_nodata(self):
        flow = branching_flow()
        weights = Grid()
        weights.setArray(np.linspace(-1, 2, 42).reshape(flow.getDims()))
        weights._geot = flow._geot
        for source in (None, weights):
            initial = None if source is None else source.readArray()
            expected = reference_accumulation(flow, initial)
            actual = flow.flowAccumulation(weights=source, nodata=False, asgrid=False)
            np.testing.assert_array_equal(actual.ravel(), expected)
            grid = flow.flowAccumulation(weights=source, nodata=True)
            expected[flow._nodata_pos] = np.iinfo(np.uint32).max
            np.testing.assert_array_equal(grid.readArray().ravel(), expected)
        np.testing.assert_array_equal(weights.readArray(), np.linspace(-1, 2, 42).reshape(flow.getDims()))

    def test_compact_lookup_with_large_cell_ids(self):
        # These cell IDs would require a 32 GiB float64 lookup raster.
        ix = np.array([2**32 - 2, 7, 2**32 - 1], dtype=np.uint32)
        ixc = np.array([2**32 - 1, 2**32 - 1, 8], dtype=np.uint32)
        receivers = receiver_positions(ix, ixc)
        np.testing.assert_array_equal(receivers, [2, 2, -1])
        np.testing.assert_array_equal(downstream_sums(receivers, np.array([2., 3., 5.])),
                                      [7., 8., 5.])

    def test_empty_flow_and_single_edge(self):
        flow = Flow()
        net = Network(flow, threshold=1)
        self.assertEqual(net._ix.size, 0)
        self.assertEqual(net._chi.size, 0)
        flow = branching_flow()
        net = Network(flow, threshold=5)
        self.assertEqual(net._ix.size, 1)
        np.testing.assert_array_equal(net._dx, net._dd)


if __name__ == '__main__':
    unittest.main()
