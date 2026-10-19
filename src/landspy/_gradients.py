"""Compact channel traversal and centered float64 gradient regressions."""

# Last modified: 2026-10-19 10:46:34 CEST (Europe/Madrid).

import numpy as np
from numba import njit, objmode

from ._network_kernels import receiver_positions


def reference_fit(x, y):
    """Original polyfit behavior for degenerate or ill-conditioned windows."""
    coefficients, residuals = np.polyfit(x, y, deg=1, full=True)[:2]
    gradient = coefficients[0]
    if gradient < 0.001:
        gradient = 0.001
    variance = y.size * y.var()
    r2 = 1 if variance == 0 else float(1 - residuals[0] / variance)
    return gradient, r2


@njit(cache=True)
def centered_fit(x, y):
    """Return slope, R2 and whether the centered fit is well conditioned."""
    count = x.size
    if count < 3:
        return 0., 0., False
    x0, y0 = x[0], y[0]
    mean_x, mean_y = 0., 0.
    max_x, max_y = 0., 0.
    for i in range(count):
        if not np.isfinite(x[i]) or not np.isfinite(y[i]):
            return 0., 0., False
        mean_x += x[i] - x0
        mean_y += y[i] - y0
        max_x = max(max_x, abs(x[i]))
        max_y = max(max_y, abs(y[i]))
    mean_x /= count
    mean_y /= count
    xx, xy, yy, spread_x, spread_y = 0., 0., 0., 0., 0.
    for i in range(count):
        dx = (x[i] - x0) - mean_x
        dy = (y[i] - y0) - mean_y
        xx += dx * dx
        xy += dx * dy
        yy += dy * dy
        spread_x = max(spread_x, abs(dx))
        spread_y = max(spread_y, abs(dy))
    tolerance = count * np.finfo(np.float64).eps * 32
    if (xx == 0 or not np.isfinite(xx + xy + yy)
            or spread_x <= max_x * tolerance
            or (yy > 0 and spread_y <= max_y * tolerance)):
        return 0., 0., False
    # NumPy's variance of a repeated noninteger can be nonzero through
    # rounding. Keep its historical policy via the reference implementation.
    if yy == 0 and (y0 != np.rint(y0) or abs(y0) > 1e12):
        return 0., 0., False
    slope = xy / xx
    if not np.isfinite(slope):
        return 0., 0., False
    residual = 0.
    for i in range(count):
        error = ((y[i] - y0) - mean_y) - slope * ((x[i] - x0) - mean_x)
        residual += error * error
    if not np.isfinite(residual):
        return 0., 0., False
    r2 = 1. if yy == 0 else 1. - residual / yy
    return max(slope, 0.001), r2, True


@njit(cache=True)
def fit_window(x, y):
    """Use the fast fit normally, preserving polyfit's exceptional cases."""
    gradient, r2, regular = centered_fit(x, y)
    if not regular:
        with objmode(gradient='float64', r2='float64'):
            gradient, r2 = reference_fit(x, y)
    return gradient, r2


def compact_topology(ix, ixc, elevations):
    """Index only channel givers and terminal receivers, never the DEM grid."""
    count = ix.size
    if not count:
        return np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64)
    receivers = receiver_positions(ix, ixc)
    missing = receivers < 0
    outlets, inverse = np.unique(ixc[missing], return_inverse=True)
    receivers[missing] = count + inverse
    if outlets.size:
        # Missing cells mapped to position zero in the legacy dense lookup.
        # Preserve its next-cell placeholder without a full-raster array.
        receivers = np.concatenate((receivers, np.full(outlets.size, receivers[0], dtype=np.int64)))
    incoming = np.zeros(receivers.size, dtype=np.bool_)
    incoming[receivers[:count]] = True
    heads = np.flatnonzero(~incoming[:count])
    # streamPoi returned heads in raster-index order before elevation sorting.
    heads = heads[np.argsort(ix[heads])]
    heads = heads[np.argsort(-elevations[heads])]
    return receivers, heads


@njit(cache=True)
def gradient_traversal(receivers, heads, x, y, npoints, fit, second_x, both,
                       skip_discarded):
    """Preserve legacy windows/branch priorities using compact node indices."""
    count = y.size
    slopes = np.zeros(receivers.size, dtype=np.float64)
    r2 = np.zeros(receivers.size, dtype=np.float64)
    other = np.zeros(receivers.size if both else 0, dtype=np.float64)
    other_r2 = np.zeros(receivers.size if both else 0, dtype=np.float64)
    capacity = max(3, min(int(2 * npoints + 1) + 2, count + 3))
    window = np.empty(capacity, dtype=np.int64)
    xbuf = np.empty(capacity, dtype=np.float64)
    ybuf = np.empty(capacity, dtype=np.float64)
    second_buf = np.empty(capacity if both else 0, dtype=np.float64)
    for head in heads:
        mid = receivers[head]
        mouth = receivers[mid]
        if mid == 0 or mid >= count or slopes[mid] != 0:
            continue
        processing = not (mouth == 0 or mouth >= count)
        window[0] = mouth if processing else mid
        window[1], window[2] = mid, head
        length = 3
        for k in range(length):
            sample = window[k] if window[k] < count else 0
            xbuf[k], ybuf[k] = x[sample], y[sample]
            if both:
                second_buf[k] = second_x[sample]
        slopes[mid], r2[mid] = fit(xbuf[:length], ybuf[:length])
        if both:
            other[mid], other_r2[mid] = fit(second_buf[:length], ybuf[:length])
        in_outlet = False
        while processing:
            if not in_outlet:
                next_cell = receivers[mouth]
                if next_cell == 0 or next_cell >= count or length >= 2 * npoints + 1:
                    for k in range(length - 1, 0, -1):
                        window[k] = window[k - 1]
                    window[0] = next_cell
                else:
                    for k in range(length - 1, -1, -1):
                        window[k + 2] = window[k]
                    window[0], window[1] = receivers[next_cell], next_cell
                    length += 2
            else:
                length -= 2
            mid = receivers[mid]
            mouth = window[0]
            if mouth == 0 or mouth >= count:
                in_outlet = True
                window[0] = window[1]
            if length <= 3:
                processing = False
            # The reference computed and discarded a regression here when a
            # tributary reached a cell already processed by a higher head.
            already_processed = slopes[mid] != 0
            if already_processed and skip_discarded:
                break
            for k in range(length):
                sample = window[k] if window[k] < count else 0
                xbuf[k], ybuf[k] = x[sample], y[sample]
                if both:
                    second_buf[k] = second_x[sample]
            gradient, score = fit(xbuf[:length], ybuf[:length])
            if already_processed:
                break
            slopes[mid], r2[mid] = gradient, score
            if both:
                other[mid], other_r2[mid] = fit(second_buf[:length], ybuf[:length])
    return slopes[:count], r2[:count], other[:count], other_r2[:count]
