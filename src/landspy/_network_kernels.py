"""Compiled topological traversals with channel-sized distance/chi storage."""

# Last modified: 2026-10-19 11:56:48 CEST (Europe/Madrid).

import numpy as np
from numba import njit


@njit(cache=True)
def accumulate(ix, ixc, values):
    """Accumulate in place, retaining the caller's integer/float dtype."""
    for n in range(ix.size):
        values[ixc[n]] += values[ix[n]]


@njit(cache=True)
def channel_mask(ix, accumulation, threshold):
    """Threshold givers without gathering a full-sized accumulation copy."""
    selected = np.empty(ix.size, dtype=np.bool_)
    for n in range(ix.size):
        selected[n] = accumulation[ix[n]] >= threshold
    return selected


@njit(cache=True)
def receiver_positions(ix, ixc):
    """Map receiver cell IDs to channel positions; outlets have position -1.

    The topologically ordered givers need not be numerically sorted. A sorted
    index and binary searches avoid a lookup raster spanning the entire DEM.
    """
    ordered = np.argsort(ix)
    positions = np.full(ix.size, -1, dtype=np.int64)
    for n in range(ix.size):
        target = ixc[n]
        low, high = 0, ix.size
        while low < high:
            middle = (low + high) // 2
            if ix[ordered[middle]] < target:
                low = middle + 1
            else:
                high = middle
        if low < ix.size and ix[ordered[low]] == target:
            positions[n] = ordered[low]
    return positions


@njit(cache=True)
def edge_lengths(ix, ixc, ncols, geot):
    """Use the same pixel-center coordinates as PRaster.cellToXY."""
    lengths = np.empty(ix.size, dtype=np.float64)
    for n in range(ix.size):
        grow, gcol = ix[n] // ncols, ix[n] % ncols
        rrow, rcol = ixc[n] // ncols, ixc[n] % ncols
        gx = geot[0] + geot[1] * gcol + geot[1] / 2
        gy = geot[3] + geot[5] * grow + geot[5] / 2
        rx = geot[0] + geot[1] * rcol + geot[1] / 2
        ry = geot[3] + geot[5] * rrow + geot[5] / 2
        lengths[n] = np.sqrt((gx - rx)**2 + (gy - ry)**2)
    return lengths


@njit(cache=True)
def downstream_sums(receivers, increments):
    """Sum downstream in reverse topological order, with zero at outlets."""
    values = np.zeros(receivers.size, dtype=np.float64)
    for n in range(receivers.size - 1, -1, -1):
        receiver = receivers[n]
        downstream = values[receiver] if receiver >= 0 else 0.0
        values[n] = downstream + increments[n]
    return values


@njit(cache=True)
def basin_channel_path(head, receivers, giver_positions, owners, channel_id):
    """Trace one legacy basin channel using only network-sized arrays."""
    path = [head]
    owners[head] = 1
    next_node = receivers[giver_positions[head]]
    owners[next_node] = 1
    flowto = -1
    while giver_positions[next_node] > 0:
        path.append(next_node)
        next_node = receivers[giver_positions[next_node]]
        if owners[next_node] == -1:
            owners[next_node] = channel_id
        else:
            flowto = owners[next_node]
            path.append(next_node)
            break
    return np.asarray(path, dtype=np.int64), flowto
