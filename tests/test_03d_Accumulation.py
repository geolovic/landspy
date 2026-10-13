"""Accumulation arithmetic and output contracts, independent of saved fixtures."""

import unittest
from unittest.mock import patch

import numpy as np

from landspy import Flow, Grid
from landspy._accumulation import accumulate_flow


class AccumulationTest(unittest.TestCase):
    def make_flow(self, dtype=np.uint32):
        flow = Flow()
        flow._size = (4, 2)
        flow._ix = np.array([0, 1, 2, 3, 4], dtype=dtype)
        flow._ixc = np.array([2, 2, 3, 5, 6], dtype=dtype)
        flow._nodata_pos = np.array([7], dtype=np.int64)
        return flow

    def test_shared_receivers_outlets_and_nodata(self):
        for index_dtype in [np.uint32, np.int32, np.int64]:
            for nodata in [False, True]:
                for asgrid in [False, True]:
                    with self.subTest(dtype=index_dtype, nodata=nodata, asgrid=asgrid):
                        flow = self.make_flow(index_dtype)
                        out = flow.flowAccumulation(nodata=nodata, asgrid=asgrid)
                        array = out.readArray() if asgrid else out
                        sentinel = np.iinfo(np.uint32).max if nodata else 0
                        expected = np.array([[1, 1, 3, 4], [1, 5, 2, sentinel]], np.uint32)
                        np.testing.assert_array_equal(array, expected)
                        self.assertEqual(array.dtype, np.uint32)
                        if asgrid:
                            self.assertEqual(out.getNodata(), sentinel if nodata else None)
                            self.assertEqual(out.getGeot(), flow.getGeot())
                            self.assertEqual(out.getCRS(), flow.getCRS())

    def test_weighted_order_precision_and_input_preservation(self):
        flow = self.make_flow()
        for dtype in [np.int16, np.float32, np.float64]:
            for nodata in [False, True]:
                with self.subTest(dtype=dtype, nodata=nodata):
                    values = np.array([[3, -2, 5, 7], [0, -1, 4, 99]], dtype=dtype)
                    if dtype == np.float64:
                        values[0] = [1e16, -1e16, 1, .25]
                    weights = Grid()
                    weights.setArray(values.copy())
                    weights._geot = flow._geot
                    expected = values.ravel().astype(np.float64)
                    for giver, receiver in zip(flow._ix, flow._ixc):
                        expected[receiver] += expected[giver]
                    expected[7] = np.iinfo(np.uint32).max if nodata else 0
                    out = flow.flowAccumulation(weights, nodata=nodata, asgrid=False)
                    np.testing.assert_array_equal(out, expected.reshape(2, 4))
                    self.assertEqual(out.dtype, np.float64)
                    np.testing.assert_array_equal(weights.readArray(), values)

    def test_unsigned_overflow(self):
        values = np.array([2**32-1, 1, 1], np.uint32)
        accumulate_flow(np.array([0, 1], np.uint32), np.array([1, 2], np.uint32), values)
        np.testing.assert_array_equal(values, np.array([2**32-1, 0, 1], np.uint32))

    def test_random_topologies_preserve_sums_exactly(self):
        rng = np.random.default_rng(38)
        givers = np.arange(511, dtype=np.int64)
        receivers = np.array([rng.integers(i+1, 512) for i in givers])
        for dtype in [np.uint32, np.float64]:
            with self.subTest(dtype=dtype):
                initial = (rng.integers(0, 2**32, 512, dtype=np.uint32)
                           if dtype == np.uint32 else rng.normal(size=512))
                expected = initial.copy()
                with np.errstate(over='ignore'):
                    for giver, receiver in zip(givers, receivers):
                        expected[receiver] += expected[giver]
                actual = initial.copy()
                accumulate_flow(givers, receivers, actual)
                np.testing.assert_array_equal(actual, expected)

    def test_nonfinite_weights_keep_numpy_behavior(self):
        flow = self.make_flow()
        weights = Grid()
        values = np.array([[np.inf, -np.inf, 0, 1], [np.nan, 1, 1, 2]])
        weights.setArray(values)
        weights._geot = flow._geot
        expected = values.ravel().copy()
        with np.errstate(invalid='ignore'):
            for giver, receiver in zip(flow._ix, flow._ixc):
                expected[receiver] += expected[giver]
        expected[7] = 0
        np.testing.assert_array_equal(flow.flowAccumulation(weights, nodata=False, asgrid=False), expected.reshape(2, 4))

    def test_empty_and_disconnected(self):
        np.testing.assert_array_equal(Flow().flowAccumulation(asgrid=False), [[1]])
        flow = self.make_flow()
        flow._ix = flow._ix[:0]
        flow._ixc = flow._ixc[:0]
        flow._nodata_pos = np.arange(8)
        np.testing.assert_array_equal(flow.flowAccumulation(nodata=False, asgrid=False), np.zeros((2, 4), np.uint32))

    def test_nodata_does_not_expand_to_coordinate_arrays(self):
        flow = self.make_flow()
        with patch('landspy.flow.np.unravel_index', side_effect=AssertionError('unnecessary coordinate arrays')):
            out = flow.flowAccumulation(nodata=False, asgrid=False)
        self.assertEqual(out[1, 3], 0)


if __name__ == '__main__':
    unittest.main()
