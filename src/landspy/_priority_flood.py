"""Eight-neighbour depression filling with a compiled Priority-Flood loop."""

import heapq

import numpy as np
from numba import njit


@njit(cache=True)
def _flood(filled, visited, pit):
    rows, cols = filled.shape
    values = filled.reshape(-1)
    # Establish the heap's native tuple type without Python objects per cell.
    heap = [(values[0], 0)]
    heapq.heappop(heap)

    # Only the raster boundary is an outlet.
    for row in range(rows):
        for col in (0, cols - 1):
            index = row * cols + col
            if not visited[index]:
                visited[index] = True
                heapq.heappush(heap, (values[index], index))
    for col in range(cols):
        for row in (0, rows - 1):
            index = row * cols + col
            if not visited[index]:
                visited[index] = True
                heapq.heappush(heap, (values[index], index))

    head = 0
    tail = 0
    while heap or head < tail:
        # Cells within a depression share its spill elevation. Process these
        # with a FIFO instead of paying for a heap operation on every cell.
        if head < tail:
            index = int(pit[head])
            head += 1
            level = values[index]
        else:
            level, index = heapq.heappop(heap)

        row = index // cols
        col = index % cols
        for dr in range(-1, 2):
            nr = row + dr
            if nr < 0 or nr >= rows:
                continue
            for dc in range(-1, 2):
                nc = col + dc
                if nc < 0 or nc >= cols:
                    continue
                neighbour = nr * cols + nc
                if visited[neighbour]:
                    continue
                visited[neighbour] = True
                if values[neighbour] <= level:
                    values[neighbour] = level
                    pit[tail] = neighbour
                    tail += 1
                else:
                    heapq.heappush(heap, (values[neighbour], neighbour))


def priority_flood(array):
    """Return a filled copy, preserving dtype and leaving input untouched.

    All cells are processed as elevations. Auxiliary arrays use
    one byte per cell for visitation and four bytes per cell for the FIFO
    (eight bytes when the raster exceeds the signed 32-bit index range).
    The heap size depends on the terrain. This is an in-memory algorithm.
    """
    if array.ndim != 2 or array.size == 0:
        raise ValueError('Priority-Flood requires a nonempty two-dimensional array')
    if array.dtype.kind not in 'biuf':
        raise TypeError('Priority-Flood requires real numeric elevations')
    if array.dtype.kind == 'f' and np.isnan(array).any():
        raise ValueError('Priority-Flood requires elevations without NaN values')
    # Numba does not support float16 arrays; float32 represents them exactly.
    dtype = np.float32 if array.dtype == np.dtype('float16') else array.dtype
    filled = np.array(array, dtype=dtype, order='C', copy=True)
    visited = np.zeros(array.size, dtype=np.bool_)
    index_dtype = np.int32 if array.size <= np.iinfo(np.int32).max else np.int64
    pit = np.empty(array.size, dtype=index_dtype)
    _flood(filled, visited, pit)
    return filled.astype(array.dtype, copy=False)
