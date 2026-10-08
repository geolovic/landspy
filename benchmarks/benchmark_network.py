"""Compare Network construction in separate processes against a Git revision.

Use --flow for a saved real Flow, or the default sparse chain in a large layout.
RSS includes imports and Flow loading; construction time excludes both and JIT
warmup. The synthetic chain is a memory stress case, not a realistic DEM.
"""

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import resource
import subprocess
import sys
import tempfile
import time

import numpy as np
from landspy import Flow, Network
from landspy._network import accumulate_downstream


def worker(args):
    cls = Network
    if args.revision:
        repo = Path(__file__).resolve().parents[1]
        source = subprocess.check_output(
            ['git', 'show', f'{args.revision}:src/landspy/network.py'], cwd=repo)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'network.py'
            path.write_bytes(source)
            spec = importlib.util.spec_from_file_location('landspy._benchmark_network', path)
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            cls = module.Network
    if args.flow:
        flow = Flow(args.flow)
        threshold = args.threshold
    else:
        flow = Flow()
        flow._size = (args.size, args.size)
        flow._ix = np.arange(args.size - 1, dtype=np.uint32)
        flow._ixc = flow._ix + 1
        flow._zx = np.arange(args.size, 1, -1, dtype=np.float64)
        flow._nodata_pos = np.array([], dtype=np.int64)
        threshold = 1 if args.threshold == 0 else args.threshold
    accumulate_downstream(np.array([0]), np.array([1]), np.array([1.]), 2)
    start = time.perf_counter()
    network = cls(flow, threshold=threshold, gradients=args.gradients)
    seconds = time.perf_counter() - start
    checksum = hashlib.sha256()
    for name in ['_ix', '_ixc', '_ax', '_zx', '_dd', '_dx', '_chi',
                 '_slp', '_ksn', '_r2slp', '_r2ksn']:
        array = getattr(network, name)
        checksum.update(name.encode())
        checksum.update(str(array.dtype).encode())
        checksum.update(array.tobytes())
    return {'revision': args.revision or 'working-tree', 'seconds': seconds,
            'peak_rss_mib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024,
            'dem_cells': flow.getNCells(), 'network_cells': network._ix.size,
            'threshold': threshold, 'gradients': args.gradients,
            'sha256': checksum.hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', default='6c351a43f46ac8c73ec98b0b1eab500c6c53318d')
    parser.add_argument('--flow')
    parser.add_argument('--size', type=int, default=4096)
    parser.add_argument('--threshold', type=int, default=0)
    parser.add_argument('--gradients', action='store_true')
    parser.add_argument('--runs', type=int, default=3)
    parser.add_argument('--worker', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--revision', help=argparse.SUPPRESS)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.worker:
        print(json.dumps(worker(args)))
        return
    if args.size < 2 or args.runs < 1:
        parser.error('size must be >= 2 and runs must be >= 1')
    common = [sys.executable, str(Path(__file__).resolve()), '--worker',
              '--size', str(args.size), '--threshold', str(args.threshold)]
    if args.flow:
        common += ['--flow', str(Path(args.flow).resolve())]
    if args.gradients:
        common += ['--gradients']
    results = []
    for _ in range(args.runs):
        for revision in [args.baseline, None]:
            command = common + (['--revision', revision] if revision else [])
            results.append(json.loads(subprocess.check_output(command, text=True)))
    if len({result['sha256'] for result in results}) != 1:
        raise RuntimeError('Network outputs differ')
    output = json.dumps({'runs': results}, indent=2) + '\n'
    if args.output:
        args.output.write_text(output)
    print(output, end='')


if __name__ == '__main__':
    main()
