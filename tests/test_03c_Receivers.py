"""Regression checks against the original dilation-based receiver rules."""

import unittest

import numpy as np
from scipy import ndimage

from landspy.flow import get_receivers


def dilation_reference(ix, elevations, cellsize, order):
    ranks = np.zeros(elevations.shape, dtype=np.int32)
    flat = ranks.ravel(order=order).copy()
    flat[ix] = np.arange(ix.size, dtype=np.int32)
    ranks = flat.reshape(elevations.shape, order=order)
    cardinal = ndimage.grey_dilation(
        ranks, footprint=np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]]))
    diagonal = ndimage.grey_dilation(
        ranks, footprint=np.array([[1, 0, 1], [0, 1, 0], [1, 0, 1]]))
    first_rank, second_rank = cardinal.ravel(order=order)[ix], diagonal.ravel(order=order)[ix]
    first, second = ix[first_rank], ix[second_rank]
    values = elevations.ravel(order=order)
    g1 = (values[ix] - values[first]) / cellsize
    g2 = (values[ix] - values[second]) / (cellsize * np.sqrt(2))
    return np.where((g1 <= g2) & (second_rank > first_rank), second, first).astype('uint32')


class ReceiversTest(unittest.TestCase):
    def test_reference_types_orders_and_edges(self):
        rng = np.random.default_rng(27)
        for shape in ((1, 1), (1, 31), (37, 1), (17, 29)):
            for dtype in ('int8', 'uint8', 'int16', 'uint16', 'int32',
                          'uint32', 'int64', 'uint64', 'float16', 'float32', 'float64', 'longdouble'):
                dt = np.dtype(dtype)
                if dt.kind in 'iu':
                    limits = np.iinfo(dt)
                    elevations = rng.choice(np.array([limits.min, limits.max, 0, 1], dtype=dt), shape)
                else:
                    elevations = rng.normal(0, 10, shape).astype(dt)
                ix = rng.permutation(elevations.size).astype('uint32')
                before = elevations.copy()
                for order in ('C', 'F'):
                    for cellsize in (27.006757161431667, np.float32(3.71), 1.0):
                        with self.subTest(shape=shape, dtype=dtype, order=order, cellsize=cellsize):
                            with np.errstate(over='ignore', invalid='ignore'):
                                expected = dilation_reference(ix, elevations, cellsize, order)
                                actual = get_receivers(ix, elevations, cellsize, order)
                            np.testing.assert_array_equal(actual, expected)
                            self.assertEqual(actual.dtype, expected.dtype)
                np.testing.assert_array_equal(elevations, before)

    def test_flat_ties_and_noncontiguous_input(self):
        rng = np.random.default_rng(45)
        for elevations in (np.ones((9, 11)), np.arange(99).reshape(9, 11)[:, ::-1]):
            ix = rng.permutation(elevations.size).astype('uint32')
            for order in ('C', 'F'):
                np.testing.assert_array_equal(
                    get_receivers(ix, elevations, 2.3, order),
                    dilation_reference(ix, elevations, 2.3, order))


if __name__ == '__main__':
    unittest.main()
