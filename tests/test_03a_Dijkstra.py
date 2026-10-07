"""Compare compiled geometric distances against the previous MCP solver."""

import unittest

import numpy as np
from skimage.graph import MCP_Geometric

from landspy._dijkstra import cost_distances, flat_distances, _distances


class DijkstraTest(unittest.TestCase):
    def test_unit_mask_float32(self):
        mask = np.ones((33, 35), dtype=bool)
        mask[8:28, 15] = False
        costs = np.where(mask, 1.0, np.inf)
        expected = MCP_Geometric(costs).find_costs([(0, 0), (32, 34)])[0]
        actual = flat_distances(mask, [(0, 0), (32, 34), (0, 0)])
        self.assertEqual(actual.dtype, np.dtype('float32'))
        np.testing.assert_allclose(actual, expected, rtol=2e-6, atol=2e-6)

    def test_unit_mask_barriers_and_no_seeds(self):
        mask = np.array([[True, False, True], [True, False, True]])
        distances = flat_distances(mask, [(0, 0)])
        np.testing.assert_array_equal(distances[:, 0], [0, 1])
        self.assertTrue(np.isinf(distances[:, 1:]).all())
        self.assertTrue(np.isinf(flat_distances(mask, [])).all())

    def test_mcp_reference(self):
        rng = np.random.default_rng(72)
        for shape in ((1, 1), (1, 35), (35, 1), (19, 23)):
            for uniform in (False, True):
                costs = (np.ones(shape) if uniform else
                         rng.uniform(0.1, 10, shape))
                costs[rng.random(shape) < 0.2] = np.inf
                costs[0, 0] = costs[-1, -1] = 1
                starts = [(0, 0), (shape[0] - 1, shape[1] - 1), (0, 0)]
                for order in ('C', 'F'):
                    with self.subTest(shape=shape, uniform=uniform, order=order):
                        surface = np.array(costs, order=order)
                        expected = MCP_Geometric(surface).find_costs(starts)[0]
                        actual = cost_distances(surface, starts)
                        np.testing.assert_allclose(actual, expected, rtol=1e-14,
                                                   atol=1e-14)
                        self.assertEqual(actual.dtype, np.dtype('float64'))

    def test_growing_heap_and_wide_indices(self):
        costs = np.ones((40, 40))
        expected = MCP_Geometric(costs).find_costs([(20, 20)])[0]
        # A tiny initial frontier forces repeated growth; int64 exercises the
        # large-raster specialization without allocating a huge raster.
        actual = _distances(costs, np.array([[20, 20]], dtype=np.int64),
                            np.full(costs.size, -1, dtype=np.int64),
                            np.empty(1, dtype=np.int64))
        np.testing.assert_allclose(actual, expected, rtol=1e-14, atol=1e-14)

    def test_unreachable_and_empty_seeds(self):
        costs = np.array([[1, np.inf, 1], [1, np.inf, 1]])
        result = cost_distances(costs, [(0, 0)])
        np.testing.assert_array_equal(result[:, 0], [0, 1])
        self.assertTrue(np.isinf(result[:, 1:]).all())
        self.assertTrue(np.isinf(cost_distances(costs, [])).all())

    def test_invalid_seed(self):
        with self.assertRaises(ValueError):
            cost_distances(np.ones((2, 2)), [(2, 0)])


if __name__ == '__main__':
    unittest.main()
