"""Receiver selection without full-raster dilations or gradient arrays."""

import numpy as np
from numba import njit


@njit(cache=True)
def _build_ranks(ix, ranks):
    for position in range(ix.size):
        ranks[ix[position]] = position


@njit(cache=True)
def _candidate_ranks(node, ranks, rows, cols, fortran):
    if fortran:
        row, col = node % rows, node // rows
        row_stride, col_stride = 1, rows
    else:
        row, col = node // cols, node % cols
        row_stride, col_stride = cols, 1
    cardinal = ranks[node]
    diagonal = cardinal
    for dr in range(-1, 2):
        nr = min(max(row + dr, 0), rows - 1)
        for dc in range(-1, 2):
            if dr == 0 and dc == 0:
                continue
            # With one-cell offsets, clipping reproduces scipy's reflect
            # boundary mode, including diagonal-to-cardinal edge reflection.
            nc = min(max(col + dc, 0), cols - 1)
            rank = ranks[nr * row_stride + nc * col_stride]
            if dr == 0 or dc == 0:
                cardinal = max(cardinal, rank)
            else:
                diagonal = max(diagonal, rank)
    return cardinal, diagonal


@njit(cache=True, error_model='numpy')
def _select(ix, values, ranks, output, rows, cols, fortran,
            cardinal_length, diagonal_length, differences, g1, g2):
    for position in range(ix.size):
        node = ix[position]
        cardinal, diagonal = _candidate_ranks(node, ranks, rows, cols, fortran)
        first, second = ix[cardinal], ix[diagonal]
        # NumPy subtracts in the input dtype, including integer overflow.
        # These tiny buffers preserve those casts and ufunc division dtypes
        # without allocating full-raster differences or gradients.
        differences[0] = values[node] - values[first]
        differences[1] = values[node] - values[second]
        g1[0] = differences[0]
        g2[0] = differences[1]
        g1[0] = g1[0] / cardinal_length
        g2[0] = g2[0] / diagonal_length
        output[position] = second if g1[0] <= g2[0] and diagonal > cardinal else first


@njit(cache=True)
def _block_candidates(ix, ranks, first, second, rows, cols, fortran):
    for position in range(ix.size):
        cardinal, diagonal = _candidate_ranks(ix[position], ranks, rows, cols, fortran)
        first[position] = cardinal
        second[position] = diagonal


def receiver_indices(ix, elevations, cellsize, order='C'):
    """Preserve rank/gradient receiver rules, including reflected boundaries.

    For contiguous input in the requested order, only ranks and uint32 output
    scale with raster size. Other layouts may need a ravel copy. Float16 and
    extended floats, unsupported by Numba arithmetic, use bounded NumPy blocks
    with compiled neighbours.
    """
    if order not in ('C', 'F'):
        raise ValueError('Receiver order must be C or F')
    elevations = np.asarray(elevations)
    ix = np.asarray(ix)
    if elevations.dtype.kind == 'b':
        raise TypeError('Boolean elevation subtraction is not supported')
    if elevations.ndim != 2 or ix.ndim != 1 or ix.size != elevations.size:
        raise ValueError('Receivers require a 2D DEM and one index per cell')
    if ix.dtype.kind not in 'iu':
        raise TypeError('Receiver indices must be integers')
    if ix.size and (ix.min() < 0 or ix.max() >= ix.size):
        raise ValueError('Receiver index outside the DEM')
    rows, cols = elevations.shape
    values = elevations.ravel(order=order)
    rank_dtype = np.int32 if ix.size <= np.iinfo(np.int32).max else np.int64
    ranks = np.zeros(ix.size, dtype=rank_dtype)
    _build_ranks(ix, ranks)
    output = np.empty(ix.size, dtype=np.uint32)
    diagonal_length = cellsize * np.sqrt(2)
    if values.dtype.kind == 'f' and values.dtype.itemsize not in (4, 8):
        # Keep unsupported floating-point precision without a full conversion.
        for start in range(0, ix.size, 65536):
            block = ix[start:start + 65536]
            first = np.empty(block.size, dtype=rank_dtype)
            second = np.empty(block.size, dtype=rank_dtype)
            _block_candidates(block, ranks, first, second, rows, cols, order == 'F')
            a, b = ix[first], ix[second]
            g1 = (values[block] - values[a]) / cellsize
            g2 = (values[block] - values[b]) / diagonal_length
            output[start:start + block.size] = np.where((g1 <= g2) & (second > first), b, a)
    else:
        # Ask NumPy for the actual ufunc promotion rules on this version.
        probe = np.zeros(1, dtype=values.dtype)
        g1_dtype = (probe / cellsize).dtype
        g2_dtype = (probe / diagonal_length).dtype
        _select(ix, values, ranks, output, rows, cols, order == 'F',
                g1_dtype.type(cellsize), g2_dtype.type(diagonal_length),
                np.empty(2, dtype=values.dtype), np.empty(1, dtype=g1_dtype),
                np.empty(1, dtype=g2_dtype))
    return output
