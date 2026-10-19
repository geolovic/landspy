"""Regression checks for compact, compiled Network gradients.

Last modified: 2026-10-19 10:46:34 CEST (Europe/Madrid).
"""
import unittest
from unittest.mock import patch
from pathlib import Path

import numpy as np

from landspy import Flow, Network
from landspy._gradients import centered_fit, compact_topology, reference_fit
from landspy._network_kernels import edge_lengths, receiver_positions, downstream_sums


def dense_reference(net, npoints, kind='slp', fit=reference_fit):
    """Frozen pre-optimization traversal, including discarded regressions."""
    x = net._chi if kind == 'ksn' else net._dx
    y = net._zx
    lookup = np.zeros(net.getNCells(), np.int64)
    lookup[net._ix] = np.arange(net._ix.size)
    heads = net.streamPoi('heads', 'IND')
    heads = heads[np.argsort(-y[lookup[heads]])]
    gradient, r2 = np.zeros(net.getNCells()), np.zeros(net.getNCells())
    for head in heads:
        mid = net._ixc[lookup[head]]
        mouth = net._ixc[lookup[mid]]
        win = [mouth, mid, head]
        processing = True
        if lookup[mid] == 0 or gradient[mid] != 0:
            continue
        elif lookup[mouth] == 0:
            processing = False
            win = [mid, mid, head]
        gradient[mid], r2[mid] = fit(x[lookup[win]], y[lookup[win]])
        in_outlet = False
        while processing:
            if not in_outlet:
                next_cell = net._ixc[lookup[mouth]]
                if lookup[next_cell] == 0:
                    win.insert(0, next_cell)
                    win.pop()
                elif len(win) < 2 * npoints + 1:
                    win.insert(0, next_cell)
                    win.insert(0, net._ixc[lookup[next_cell]])
                else:
                    win.insert(0, next_cell)
                    win.pop()
            else:
                win.pop()
                win.pop()
            mid = net._ixc[lookup[mid]]
            mouth = win[0]
            if lookup[mouth] == 0:
                in_outlet = True
                win[0] = win[1]
            if len(win) <= 3:
                processing = False
            g, score = fit(x[lookup[win]], y[lookup[win]])
            if gradient[mid] == 0:
                gradient[mid], r2[mid] = g, score
            else:
                processing = False
    return gradient[net._ix], r2[net._ix]


def branching_network():
    net = Network()
    net._size = (8, 8)
    net._geot = (100., 27., 0., 200., 0., -23.)
    net._ix = np.array([1, 9, 2, 10, 17, 18, 25, 26, 33, 34, 41, 42], dtype='uint32')
    net._ixc = np.array([10, 17, 10, 18, 25, 26, 33, 34, 41, 42, 49, 50], dtype='uint32')
    # Equal-elevation heads exercise raster order as the tie-breaking input.
    net._zx = np.array([100., 100., 100., 90., 90., 80., 70., 70., 60., 60., 50., 50.])
    net._ax = np.arange(1, 13, dtype='uint32')
    net._dd = edge_lengths(net._ix, net._ixc, 8, net._geot)
    net._dx = downstream_sums(receiver_positions(net._ix, net._ixc), net._dd)
    net.calculateChi()
    return net


class LinearFitTest(unittest.TestCase):
    def test_centered_fit_matches_polyfit(self):
        net = Network()
        rng = np.random.default_rng(913)
        for count in (3, 5, 11, 31):
            for offset in (0., 1000., 100000.):
                for slope in (-.2, 0., .001, .5, 10.):
                    with self.subTest(count=count, offset=offset, slope=slope):
                        x = offset + np.arange(count) * 27.
                        y = 300. + slope * (x - offset) + rng.normal(0., .01, count)
                        np.testing.assert_allclose(net.polynomial_fit(x, y), reference_fit(x, y),
                                                   rtol=1e-9, atol=1e-11)

    def test_reference_fallback_and_exception_behavior(self):
        net = Network()
        cases = [(np.arange(3.), np.array([1., np.nan, 2.])),
                 (np.ones(3), np.ones(3)),
                 (np.ones(3), np.arange(3.)),
                 (np.arange(3.), np.full(3, .1)),
                 (np.arange(3.), np.full(3, 1e16)),
                 (np.arange(2.), np.ones(2)),
                 (np.arange(2.), np.arange(2.)),
                 (1e12 + np.arange(5) * .001, np.arange(5.)),
                 (np.arange(5.), 1e12 + np.arange(5) * .001)]
        for x, y in cases:
            with self.subTest(x=x, y=y):
                self.assertFalse(centered_fit(x, y)[2])
                try:
                    expected = reference_fit(x, y)
                except Exception as error:
                    with self.assertRaises(type(error)):
                        net.polynomial_fit(x, y)
                else:
                    np.testing.assert_allclose(net.polynomial_fit(x, y), expected,
                                               rtol=0, atol=0, equal_nan=True)

    def test_floor_constant_elevations_and_readonly_inputs(self):
        net = Network()
        x = np.arange(11.)
        y = np.full(11, 100.)
        x.flags.writeable = y.flags.writeable = False
        self.assertEqual(net.polynomial_fit(x, y), (.001, 1.))
        np.testing.assert_allclose(net.polynomial_fit(x[::-1], np.arange(11.)),
                                   reference_fit(x[::-1], np.arange(11.)), atol=1e-12)


class GradientsTraversalTest(unittest.TestCase):
    def assert_reference(self, net, npoints, r2_atol=1e-11):
        for kind, gname, rname in (('slp', '_slp', '_r2slp'), ('ksn', '_ksn', '_r2ksn')):
            expected, expected_r2 = dense_reference(net, npoints, kind)
            self.assertIsNone(net.calculateGradients(npoints, kind))
            np.testing.assert_allclose(getattr(net, gname), expected, rtol=1e-9, atol=1e-11)
            np.testing.assert_allclose(getattr(net, rname), expected_r2, rtol=1e-9, atol=r2_atol)
            np.testing.assert_array_equal(getattr(net, gname) == 0, expected == 0)
            self.assertEqual(getattr(net, '_ksn_np' if kind == 'ksn' else '_slp_np'), npoints)

    def test_confluences_tied_heads_and_two_outlets(self):
        for npoints in (0, 1, 2, 5, 10, 1000000):
            with self.subTest(npoints=npoints):
                self.assert_reference(branching_network(), npoints)

    def test_real_flow_float64_and_legacy_float32(self):
        path = Path(__file__).resolve().parent / 'data/in/small25_fd.tif'
        for dtype in ('float32', 'float64'):
            for npoints in (0, 2, 5, 10):
                with self.subTest(dtype=dtype, npoints=npoints):
                    net = Network(Flow(str(path)), threshold=10)
                    net._zx = net._zx.astype(dtype)
                    # The original used float32 y.var() on serialized Flow data.
                    self.assert_reference(net, npoints, 1e-7 if dtype == 'float32' else 1e-11)

    def test_fused_equals_separate_and_constructor(self):
        net = branching_network()
        net.calculateGradients(5, 'slp')
        net.calculateGradients(5, 'ksn')
        expected = [getattr(net, name).copy() for name in ('_slp', '_r2slp', '_ksn', '_r2ksn')]
        net._calculate_gradients(5, 'slp', both=True)
        for name, values in zip(('_slp', '_r2slp', '_ksn', '_r2ksn'), expected):
            np.testing.assert_array_equal(getattr(net, name), values)
        path = Path(__file__).resolve().parent / 'data/in/small25_fd.tif'
        flow = Flow(str(path))
        constructed = Network(flow, threshold=100, gradients=True)
        separate = Network(flow, threshold=100)
        separate.calculateGradients(5, 'slp')
        separate.calculateGradients(5, 'ksn')
        for name in ('_slp', '_r2slp', '_ksn', '_r2ksn'):
            np.testing.assert_array_equal(getattr(constructed, name), getattr(separate, name))

    def test_no_full_dem_lookup_and_unmodified_inputs(self):
        net = branching_network()
        expected = dense_reference(net, 5)
        before = {name: getattr(net, name).copy() for name in ('_ix', '_ixc', '_zx', '_dx', '_chi')}
        net._size = (1000000, 1000000)
        net.calculateGradients(5)
        np.testing.assert_allclose(net._slp, expected[0], rtol=1e-9, atol=1e-11)
        np.testing.assert_allclose(net._r2slp, expected[1], rtol=1e-9, atol=1e-11)
        for name, values in before.items():
            np.testing.assert_array_equal(getattr(net, name), values)
        receivers, heads = compact_topology(net._ix, net._ixc, net._zx)
        self.assertEqual(receivers.size, net._ix.size + 2)

    def test_empty_and_short_networks(self):
        for count in (0, 1, 2, 3):
            net = Network()
            net._size = (max(count + 1, 1), 1)
            net._ix = np.arange(count, dtype='uint32')
            net._ixc = net._ix + 1
            net._zx = np.arange(count, 0, -1, dtype='float64')
            net._dx = net._zx.copy()
            net._chi = net._zx.copy()
            self.assert_reference(net, 5)
        self.assert_reference(Network(), 0)

    def test_custom_fitter_is_honored_with_independent_windows(self):
        net = branching_network()
        seen = []
        expected_seen = []
        def custom(x, y):
            seen.append((x, y))
            return float(len(seen)), .75
        def reference_custom(x, y):
            expected_seen.append((x, y))
            return float(len(expected_seen)), .75
        expected = dense_reference(net, 5, fit=reference_custom)
        net.polynomial_fit = custom
        net.calculateGradients(5)
        np.testing.assert_array_equal(net._slp, expected[0])
        np.testing.assert_array_equal(net._r2slp, expected[1])
        self.assertEqual(len(seen), len(expected_seen))
        for (x, y), (rx, ry) in zip(seen, expected_seen):
            np.testing.assert_array_equal(x, rx)
            np.testing.assert_array_equal(y, ry)
        for n in range(len(seen) - 1):
            self.assertFalse(np.shares_memory(seen[n][0], seen[n + 1][0]))

    def test_class_overrides_and_constructor_dispatch(self):
        path = Path(__file__).resolve().parent / 'data/in/small25_fd.tif'
        flow = Flow(str(path))
        seen = []
        def custom_fit(self, x, y):
            seen.append(len(x))
            return 2., .75
        with patch.object(Network, 'polynomial_fit', custom_fit):
            net = Network(flow, threshold=100, gradients=True)
        self.assertTrue(seen)
        np.testing.assert_array_equal(net._slp[net._slp != 0], 2.)
        np.testing.assert_array_equal(net._ksn[net._ksn != 0], 2.)
        calls = []
        def custom_gradients(self, npoints, kind='slp'):
            calls.append((npoints, kind))
        with patch.object(Network, 'calculateGradients', custom_gradients):
            Network(flow, threshold=100, gradients=True)
        self.assertEqual(calls, [(5, 'slp'), (5, 'ksn')])


if __name__ == '__main__':
    unittest.main()
