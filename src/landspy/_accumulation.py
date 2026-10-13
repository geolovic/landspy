"""Sequential accumulation over topologically ordered drainage edges."""

from numba import njit


@njit(cache=True)
def accumulate_flow(givers, receivers, values):
    """Accumulate in place, preserving edge order and the buffer's dtype.

    Shared receivers and downstream dependencies require a sequential loop.
    No fastmath: weighted float64 sums keep the original addition order;
    uint32 assignments retain the original unsigned overflow behavior.
    """
    for n in range(givers.size):
        values[receivers[n]] += values[givers[n]]
