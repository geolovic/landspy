"""Compiled eight-neighbour cost distances without traceback or edge maps."""

import numpy as np
from numba import njit


@njit(cache=True)
def _less(a, b, distance):
    return distance[a] < distance[b] or (distance[a] == distance[b] and a < b)


@njit(cache=True)
def _sift_up(heap, positions, distance, slot):
    node = heap[slot]
    while slot:
        parent = (slot - 1) // 2
        if not _less(node, heap[parent], distance):
            break
        heap[slot] = heap[parent]
        positions[heap[slot]] = slot
        slot = parent
    heap[slot] = node
    positions[node] = slot


@njit(cache=True)
def _distance_kernel(costs, starts, positions, heap, distance, uniform):
    rows, cols = costs.shape
    count = 0
    for k in range(starts.shape[0]):
        node = starts[k, 0] * cols + starts[k, 1]
        if positions[node] != -1:
            continue
        distance[node] = 0.0
        heap[count] = node
        positions[node] = count
        _sift_up(heap, positions, distance, count)
        count += 1

    diagonal = np.sqrt(2.0)
    while count:
        node = heap[0]
        count -= 1
        positions[node] = -2  # Finalized: never enter the heap again.
        if count:
            last = heap[count]
            slot = 0
            while 2 * slot + 1 < count:
                child = 2 * slot + 1
                if child + 1 < count and _less(heap[child + 1], heap[child], distance):
                    child += 1
                if not _less(heap[child], last, distance):
                    break
                heap[slot] = heap[child]
                positions[heap[slot]] = slot
                slot = child
            heap[slot] = last
            positions[last] = slot

        row, col = node // cols, node % cols
        if uniform and not costs[row, col]:
            continue
        old_cost = 1.0 if uniform else costs[row, col]
        for dr in range(-1, 2):
            nr = row + dr
            if nr < 0 or nr >= rows:
                continue
            for dc in range(-1, 2):
                nc = col + dc
                if nc < 0 or nc >= cols or (dr == 0 and dc == 0):
                    continue
                neighbour = nr * cols + nc
                if positions[neighbour] == -2:
                    continue
                length = diagonal if dr != 0 and dc != 0 else 1.0
                if uniform:
                    if not costs[nr, nc]:
                        continue
                    # Round both the step and accumulated distance in float32,
                    # before comparing/updating heap priorities.
                    candidate = np.float32(distance[node] + np.float32(length))
                else:
                    new_cost = costs[nr, nc]
                    if not np.isfinite(new_cost) or new_cost < 0:
                        continue
                    # Preserve MCP_Geometric's floating-point operation order.
                    candidate = distance[node] + length * 0.5 * (old_cost + new_cost)
                if candidate >= distance[neighbour] or not np.isfinite(candidate):
                    continue
                distance[neighbour] = candidate
                slot = positions[neighbour]
                if slot == -1:
                    if count == heap.size:
                        grown = np.empty(min(costs.size, 2 * heap.size), dtype=heap.dtype)
                        grown[:count] = heap[:count]
                        heap = grown
                    slot = count
                    heap[slot] = neighbour
                    positions[neighbour] = slot
                    count += 1
                _sift_up(heap, positions, distance, slot)
    return distance.reshape((rows, cols))


@njit(cache=True)
def _distances(costs, starts, positions, heap):
    distance = np.full(costs.size, np.inf, dtype=np.float64)
    return _distance_kernel(costs, starts, positions, heap, distance, False)


def cost_distances(costs, starts):
    """Minimum geometric costs from zero-cost seeds on a 2D friction surface.

    Eight neighbours use mean endpoint friction times length (1 or sqrt(2)).
    Infinite and negative costs block entry. Distances are float64; native
    heap and position indices use int32 when the raster fits, else int64.
    The heap has one entry per frontier cell and grows only as needed.
    """
    costs = np.asarray(costs, dtype=np.float64)
    if costs.ndim != 2 or costs.size == 0:
        raise ValueError('Dijkstra requires a nonempty two-dimensional surface')
    starts = np.asarray(starts, dtype=np.int64).reshape((-1, 2))
    if (np.any(starts < 0) or np.any(starts[:, 0] >= costs.shape[0])
            or np.any(starts[:, 1] >= costs.shape[1])):
        raise ValueError('Dijkstra seed outside the friction surface')
    index_dtype = np.int32 if costs.size <= np.iinfo(np.int32).max else np.int64
    positions = np.full(costs.size, -1, dtype=index_dtype)
    heap = np.empty(min(costs.size, max(16, len(starts))), dtype=index_dtype)
    return _distances(costs, starts, positions, heap)


def flat_distances(flats, starts):
    """Float32 unit-friction distances using only a traversability mask."""
    flats = np.asarray(flats, dtype=np.bool_)
    if flats.ndim != 2 or flats.size == 0:
        raise ValueError('Dijkstra requires a nonempty two-dimensional mask')
    starts = np.asarray(starts, dtype=np.int64).reshape((-1, 2))
    if (np.any(starts < 0) or np.any(starts[:, 0] >= flats.shape[0])
            or np.any(starts[:, 1] >= flats.shape[1])):
        raise ValueError('Dijkstra seed outside the flat mask')
    index_dtype = np.int32 if flats.size <= np.iinfo(np.int32).max else np.int64
    positions = np.full(flats.size, -1, dtype=index_dtype)
    heap = np.empty(min(flats.size, max(16, len(starts))), dtype=index_dtype)
    distance = np.full(flats.size, np.inf, dtype=np.float32)
    return _distance_kernel(flats, starts, positions, heap, distance, True)
