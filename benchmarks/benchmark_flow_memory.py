"""Compare complete Flow construction with a previous Git revision.

python benchmarks/benchmark_flow_memory.py --size 14480 --baseline 3390f70
Each revision runs sequentially in its own process. Requires local Git history.
"""

import argparse
import hashlib
import json
import os
import resource
import subprocess
import sys
import time
import types


def worker(args):
    import numpy as np
    from landspy import DEM, Flow

    source = subprocess.check_output(
        ['git', 'show', args.baseline + ':src/landspy/flow.py'],
        cwd=os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        universal_newlines=True)
    baseline = types.ModuleType('landspy._memory_baseline')
    baseline.__package__ = 'landspy'
    exec(compile(source, '<baseline-flow>', 'exec'), baseline.__dict__)
    constructor = baseline.Flow if args.worker == 'before' else Flow

    # Warm the compiled fill kernel separately from the timed construction.
    warmup = DEM()
    warmup.setArray(np.ones((3, 3), dtype='float32'))
    warmup.fill()

    dem = DEM()
    dem._array = np.random.default_rng(38).integers(
        0, 2000, (args.size, args.size), dtype='int16').astype('float32')
    dem._size = (args.size, args.size)
    dem._tipo = 'float32'
    start = time.perf_counter()
    flow = constructor(dem)
    elapsed = time.perf_counter() - start
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_mib = peak / (1024**2 if sys.platform == 'darwin' else 1024)
    arrays = {name: getattr(flow, name) for name in ('_ix', '_ixc', '_zx')}
    return {'mode': args.worker, 'baseline': args.baseline,
            'input_MiB': dem.readArray().nbytes / 1024**2,
            'seconds': round(elapsed, 3), 'peak_MiB': round(peak_mib, 1),
            'checksums': {name: hashlib.sha256(memoryview(array)).hexdigest()
                          for name, array in arrays.items()},
            'dtypes': {name: str(array.dtype) for name, array in arrays.items()}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--size', type=int, default=2048)
    parser.add_argument('--baseline', default='3390f70')
    parser.add_argument('--worker', choices=('before', 'after'), help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.size < 1:
        parser.error('size must be positive')
    if args.worker:
        print(json.dumps(worker(args)))
        return
    results = []
    for mode in ('before', 'after'):
        output = subprocess.check_output([
            sys.executable, os.path.abspath(__file__), '--worker', mode,
            '--size', str(args.size), '--baseline', args.baseline],
            universal_newlines=True)
        result = json.loads(output)
        results.append(result)
        print(json.dumps(result), flush=True)
    assert results[0]['checksums'] == results[1]['checksums'], 'Flow arrays differ'
    assert results[0]['dtypes'] == results[1]['dtypes'], 'Flow dtypes differ'


if __name__ == '__main__':
    main()
