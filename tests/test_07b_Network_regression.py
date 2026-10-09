"""Centered regressions against the original NumPy SVD calculation."""

import unittest
from unittest.mock import patch
import numpy as np
from landspy import Network
from landspy._network import linear_fit


def reference(x, y):
    coefficients, residuals = np.polyfit(x, y, 1, full=True)[:2]
    r2 = 1. if y.var() == 0. else float(1. - residuals[0] / (y.size * y.var()))
    return max(coefficients[0], .001), r2


class NetworkRegressionTest(unittest.TestCase):
    def test_random_windows_against_svd(self):
        rng = np.random.default_rng(418)
        network = Network()
        for size in [3, 5, 11, 21, 41]:
            for _ in range(30):
                x = np.cumsum(rng.uniform(.1, 20., size))
                y = rng.uniform(-100., 100., size)
                np.testing.assert_allclose(network.polynomial_fit(x, y), reference(x, y),
                                           rtol=1e-10, atol=1e-12)

    def test_flat_negative_and_repeated_points(self):
        for x, y in [(np.arange(5.), np.ones(5)),
                     (np.arange(5.), -np.arange(5.)),
                     (np.array([0., 0., 1.]), np.array([0., 0., 2.]))]:
            np.testing.assert_allclose(Network().polynomial_fit(x, y), reference(x, y),
                                       rtol=1e-10, atol=1e-12)

    def test_loaded_float32_elevations_keep_r2_denominator(self):
        rng = np.random.default_rng(67)
        for size in [3, 5, 11, 21]:
            x = np.cumsum(rng.uniform(.1, 20., size))
            y = rng.uniform(-100., 100., size).astype(np.float32)
            np.testing.assert_allclose(Network().polynomial_fit(x, y), reference(x, y),
                                       rtol=1e-10, atol=1e-12)

    def test_ill_conditioned_windows_keep_polyfit(self):
        for x, y in [(1e12 + np.arange(5.), np.arange(5.)),
                     (np.arange(5.), 1e12 + np.arange(5.)),
                     (np.ones(5), np.ones(5)),
                     (np.array([0., 1., np.nan]), np.ones(3))]:
            self.assertFalse(linear_fit(x, y)[2])
        x = 1e12 + np.arange(5.)
        y = np.arange(5.)
        with patch('numpy.polyfit', wraps=np.polyfit) as original:
            result = Network().polynomial_fit(x, y)
            original.assert_called_once()
        np.testing.assert_array_equal(result, reference(x, y))

    def test_short_windows_keep_original_exception(self):
        with self.assertRaises(IndexError):
            Network().polynomial_fit(np.array([0., 1.]), np.array([0., 1.]))


if __name__ == '__main__':
    unittest.main()
