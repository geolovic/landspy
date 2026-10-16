"""Compiled drainage-basin labeling and thresholded outlet selection."""

# Last modified: 2026-10-16 11:38:32 CEST (Europe/Madrid).

import numpy as np
from numba import njit


@njit(cache=True)
def outlet_mask(ix, ixc, accumulation, threshold):
    """Return receivers without outgoing selected edges, in raster space.

    Match streamPoi's strict accumulation > threshold rule without building
    its sparse adjacency matrix, degree array, or selected-edge copies.
    """
    outlets = np.zeros(accumulation.size, dtype=np.bool_)
    # Givers precede their receivers. Mark a candidate when it receives flow;
    # its own later edge clears it again unless it is a selected-network sink.
    for n in range(ix.size):
        if accumulation[ix[n]] > threshold:
            outlets[ix[n]] = False
            outlets[ixc[n]] = True
    return outlets


@njit(cache=True)
def label_basins(ix, ixc, basins, all_basins):
    """Propagate int64 basin labels upstream in place, preserving seed IDs."""
    nbasins = 0
    for n in range(ix.size - 1, -1, -1):
        giver, receiver = ix[n], ixc[n]
        if all_basins:
            if basins[receiver] == 0:
                nbasins += 1
                basins[receiver] = nbasins
            basins[giver] = basins[receiver]
        elif basins[receiver] != 0 and basins[giver] == 0:
            basins[giver] = basins[receiver]
