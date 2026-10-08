"""Compact topological traversals for drainage networks."""

import numpy as np
from numba import njit


def compact_nodes(ix, ixc):
    """Return sorted cell IDs and edge endpoints in compact node coordinates."""
    nodes, inverse = np.unique(np.concatenate((ix, ixc)), return_inverse=True)
    return nodes, inverse[:ix.size], inverse[ix.size:]


@njit(cache=True)
def accumulate_downstream(givers, receivers, increments, node_count):
    """Preserve the reverse topological addition order, without a DEM-sized array."""
    values = np.zeros(node_count, np.float64)
    for n in range(givers.size - 1, -1, -1):
        values[givers[n]] = values[receivers[n]] + increments[n]
    return values[givers]
