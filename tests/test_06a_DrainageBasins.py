"""Basin-label equivalence, outlet ordering and seed-ID regressions."""

# Last modified: 2026-10-16 11:38:32 CEST (Europe/Madrid).

import unittest

import numpy as np

from landspy import Flow
from landspy.flow import FlowError
from landspy._basins import outlet_mask


def branching_flow(dtype='uint32'):
    flow = Flow()
    flow._size = (5, 5)
    flow._geot = (100., 2., 0., 200., 0., -3.)
    flow._proj = 'test projection'
    flow._ix = np.array([0, 1, 5, 6, 7, 12, 16, 17, 23], dtype=dtype)
    flow._ixc = np.array([6, 6, 6, 12, 12, 17, 17, 22, 24], dtype=dtype)
    flow._zx = np.arange(9., 0., -1)
    flow._nodata_pos = np.setdiff1d(np.arange(25), np.union1d(flow._ix, flow._ixc))
    return flow


def reference_labels(flow, seeds=(), ids=None):
    """Original Python traversal, including single-seed and empty fallback rules."""
    labels = np.zeros(flow.getNCells(), dtype=np.int64)
    if len(seeds):
        if len(seeds) == 1:
            labels[seeds] = 1
        else:
            for n, seed in enumerate(seeds):
                labels[seed] = ids[n] if ids is not None else n + 1
        for giver, receiver in zip(flow._ix[::-1], flow._ixc[::-1]):
            if labels[receiver] != 0 and labels[giver] == 0:
                labels[giver] = labels[receiver]
    else:
        count = 0
        for giver, receiver in zip(flow._ix[::-1], flow._ixc[::-1]):
            if labels[receiver] == 0:
                count += 1
                labels[receiver] = count
            labels[giver] = labels[receiver]
    return labels.reshape(flow.getDims())


def coordinates(flow, indices, ids=None):
    x, y = flow.cellToXY(*flow.indToCell(np.asarray(indices, dtype=np.int64)))
    return np.column_stack((x, y) if ids is None else (x, y, ids))


class DrainageBasinsRegressionTest(unittest.TestCase):
    def test_all_basins_exact_ids_and_metadata(self):
        for dtype in ('uint32', 'int32', 'int64'):
            with self.subTest(dtype=dtype):
                flow = branching_flow(dtype)
                original = flow._ix.copy(), flow._ixc.copy()
                grid = flow.drainageBasins(min_area=0)
                np.testing.assert_array_equal(grid.readArray(), reference_labels(flow))
                self.assertEqual(grid.readArray().dtype, np.dtype('int64'))
                self.assertEqual(grid.getNodata(), 0)
                self.assertEqual(grid.getGeot(), flow.getGeot())
                self.assertEqual(grid.getSize(), flow.getSize())
                self.assertEqual(grid.getCRS(), flow.getCRS())
                np.testing.assert_array_equal(flow._ix, original[0])
                np.testing.assert_array_equal(flow._ixc, original[1])

    def test_minimum_area_and_strict_threshold(self):
        for dtype in ('uint32', 'int32', 'int64'):
            flow = branching_flow(dtype)
            accumulation = flow.flowAccumulation(nodata=False, asgrid=False).ravel()
            for min_area in (.001, .04, .08, .2, 1., 2.):
                with self.subTest(dtype=dtype, min_area=min_area):
                    threshold = int(flow.getNCells() * min_area)
                    seeds = flow.streamPoi(threshold, kind='outlets', coords='IND')
                    actual_seeds = np.flatnonzero(outlet_mask(flow._ix, flow._ixc, accumulation, threshold))
                    np.testing.assert_array_equal(actual_seeds, seeds)
                    actual = flow.drainageBasins(min_area=min_area, asgrid=False)
                    np.testing.assert_array_equal(actual, reference_labels(flow, seeds))

    def test_custom_duplicate_nested_and_isolated_outlets(self):
        flow = branching_flow()
        cases = [([6, 22], [31, 42]), ([22, 6], [42, 31]),
                 ([6, 6, 22], [1, 2, 3]), ([6, 22], [0, -7]),
                 ([6, 22], [3.9, 2**40]), ([4, 24], [70, 80]),
                 ([22], [999]), ([], [])]
        for seeds, ids in cases:
            with self.subTest(seeds=seeds, ids=ids):
                outlets = coordinates(flow, seeds, ids)
                before = outlets.copy()
                result = flow.drainageBasins(outlets, asgrid=False)
                np.testing.assert_array_equal(result, reference_labels(flow, seeds, ids))
                np.testing.assert_array_equal(outlets, before)

    def test_coordinate_forms_and_default_ids(self):
        flow = branching_flow()
        single = coordinates(flow, [22])[0]
        for outlets in (single.tolist(), tuple(single), single.reshape(1, 2)):
            np.testing.assert_array_equal(flow.drainageBasins(outlets, asgrid=False),
                                          reference_labels(flow, [22]))
        outlets = coordinates(flow, [24, 22])
        np.testing.assert_array_equal(flow.drainageBasins(outlets, asgrid=False),
                                      reference_labels(flow, [24, 22]))

    def test_invalid_coordinates_and_empty_flow(self):
        flow = branching_flow()
        for outlets in ([95., 199.], np.array([[95., 199.]])):
            with self.assertRaises(FlowError):
                flow.drainageBasins(outlets)
        for min_area in (0, .005):
            result = Flow().drainageBasins(min_area=min_area, asgrid=False)
            np.testing.assert_array_equal(result, np.zeros((1, 1), dtype='int64'))


if __name__ == '__main__':
    unittest.main()
