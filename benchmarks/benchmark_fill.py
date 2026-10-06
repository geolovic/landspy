"""Measure Priority-Flood in a separate process (Linux/macOS).

Example: python benchmarks/benchmark_fill.py --size 4096
Use --size 14480 for an approximately 800 MiB DEM.
Peak RSS includes interpreter, imports, warmup, input and fill work arrays.
"""

import argparse
import hashlib
import json
import os
import resource
import subprocess
import sys
import time


def worker(args):
    import numpy as np
    from landspy import DEM

    # Warm up separately from the timed fill. Report the cost rather than
    # silently counting JIT compilation as terrain processing time.
    warmup = DEM()
    warmup.setArray(np.ones((3, 3), dtype='float32'))
    start = time.perf_counter()
    warmup.fill()
    warmup_seconds = time.perf_counter() - start

    dem = DEM()
    dem._array = np.empty((args.size, args.size), dtype='float32')
    dem._size = (args.size, args.size)
    dem._tipo = 'float32'
    if args.terrain == 'bowl':
        dem._array.fill(0)
        dem._array[0, :] = dem._array[-1, :] = 100
        dem._array[:, 0] = dem._array[:, -1] = 100
    else:
        rng = np.random.default_rng(7)
        for row in range(0, args.size, 64):
            block = dem._array[row:row + 64]
            block[:] = rng.integers(0, 2000, block.shape, dtype=np.int16)
    start = time.perf_counter()
    result = dem.fill(inplace=True)
    elapsed = time.perf_counter() - start
    peak_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_mib = peak_rss / (1024**2 if sys.platform == 'darwin' else 1024)

    digest = hashlib.sha256()
    for values in result.readArray():
        digest.update(memoryview(values))
        if args.terrain == 'bowl':
            assert np.all(values == 100)
    return {'method': 'priority-flood', 'size': args.size, 'dtype': 'float32',
            'terrain': args.terrain,
            'input_MiB': dem.readArray().nbytes / 1024**2,
            'warmup_seconds': round(warmup_seconds, 3),
            'fill_seconds': round(elapsed, 3), 'peak_rss_MiB': round(peak_mib, 1),
            'sha256': digest.hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--size', type=int, default=4096)
    parser.add_argument('--terrain', choices=('random', 'bowl'), default='random')
    parser.add_argument('--worker', action='store_true', help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.size < 1:
        parser.error('size must be positive')
    if args.worker:
        print(json.dumps(worker(args)))
        return
    output = subprocess.check_output([
        sys.executable, os.path.abspath(__file__), '--worker',
        '--size', str(args.size), '--terrain', args.terrain], universal_newlines=True)
    print(json.dumps(json.loads(output)))


if __name__ == '__main__':
    main()
