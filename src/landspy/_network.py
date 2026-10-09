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


@njit(cache=True)
def linear_fit(x, y, normalization=-1.):
    """Centered degree-one regression; flag exceptional windows for polyfit.

    No fastmath: the variance and residual sums retain their operation order.
    Returns (gradient, R2, valid), with gradient floored at 0.001 and R2
    computed from the unfloored fit. A nonnegative normalization replaces
    the y sum-of-squares denominator; a negative value uses the centred sum.
    Ill-conditioned, non-finite or short inputs return (0, 0, False); the
    caller, not this helper, must perform any NumPy SVD fallback.
    """
    count = x.size
    if count < 3 or y.size != count:
        return 0., 0., False
    sx = 0.
    sy = 0.
    largest_x = 0.
    largest_y = 0.
    lo = x[0]
    hi = x[0]
    lo_y = y[0]
    hi_y = y[0]
    for i in range(count):
        if not np.isfinite(x[i]) or not np.isfinite(y[i]):
            return 0., 0., False
        sx += x[i]
        sy += y[i]
        largest_x = max(largest_x, abs(x[i]))
        largest_y = max(largest_y, abs(y[i]))
        lo = min(lo, x[i])
        hi = max(hi, x[i])
        lo_y = min(lo_y, y[i])
        hi_y = max(hi_y, y[i])
    # Large offsets relative to spread can change polyfit's effective rank
    # and amplify its rounding errors. Preserve that behavior via the fallback.
    if hi - lo <= largest_x * 1e-4:
        return 0., 0., False
    if hi_y != lo_y and hi_y - lo_y <= largest_y * 1e-4:
        return 0., 0., False
    mx = sx / count
    my = sy / count
    xx = 0.
    xy = 0.
    yy = 0.
    for i in range(count):
        dx = x[i] - mx
        dy = y[i] - my
        xx += dx * dx
        xy += dx * dy
        yy += dy * dy
    if xx == 0. or not np.isfinite(xx) or not np.isfinite(xy) or not np.isfinite(yy):
        return 0., 0., False
    gradient = xy / xx
    residual = 0.
    for i in range(count):
        difference = (y[i] - my) - gradient * (x[i] - mx)
        residual += difference * difference
    denominator = yy if normalization < 0. else normalization
    r2 = 1. if denominator == 0. else 1. - residual / denominator
    if not np.isfinite(gradient) or not np.isfinite(r2):
        return 0., 0., False
    return max(gradient, 0.001), r2, True
